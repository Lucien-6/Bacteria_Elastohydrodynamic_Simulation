function dZ = ode_coupled_rhs(t, Z, params)
%ODE_COUPLED_RHS  Core ODE right-hand side for coupled body-pili system
%
%  Description:
%    Computes dZ/dt for the coupled bacterial body + elastic pili system.
%    Evaluation steps:
%      1-7. Build and solve the coupled linear system (delegated to
%           build_coupled_lin_sys, which handles geometry reconstruction,
%           global mobility/resistance, Q/B matrix assembly, and
%           clamped boundary conditions)
%      8.   Compute elastic restoring moments for each pilus
%      8b.  Subtract per-pili active internal moments (if enabled)
%      8c.  Subtract background flow contribution (if enabled)
%      9.   Solve  lin_sys * z_dot = R_global
%     10.   Pack state rates into dZ
%
%  Syntax:
%    dZ = ode_coupled_rhs(t, Z, params)
%
%  Inputs:
%    t      - Current dimensionless time
%    Z      - [n x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    dZ - [n x 1] Time derivative of state vector
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06

%% Extract frequently used parameters
N_segs = params.N_segments;
sigma  = params.sigma;
Eg     = params.Eg;
g_hat  = params.g_hat(:);
M_pili = params.M_pili;

% Determine root coupling type for each pilus (default: clamped)
rc = cell(M_pili, 1);
for k = 1:M_pili
    if isfield(params, 'root_coupling') && numel(params.root_coupling) >= k
        rc{k} = params.root_coupling{k};
    else
        rc{k} = 'clamped';
    end
end

% Cross-call state for motor speed-limiting (free_torsion + constant motor).
%
% MATLAB structs are passed by value, so any field assignment to `params`
% inside this function would be invisible to subsequent ODE calls (the
% solver re-passes the original `params` every time).  We therefore use a
% persistent variable to remember the most recent dpsi_1 rate.
%
% NOTE on numerical semantics:
%   ode15s is an adaptive multi-step implicit solver that may evaluate the
%   RHS multiple times per accepted step (Newton iterations, rejected
%   steps, error estimation, dense-output interpolation).  DPSI1_PREV
%   therefore caches dpsi_1 from the LAST RHS EVALUATION, which is in
%   general NOT the previous physically-accepted time step.  The
%   speed-limiting model
%       tau_eff = tau_motor * (1 - |omega| / omega_max)
%   is consequently applied as an EXPLICIT LAGGED APPROXIMATION rather
%   than an implicit nonlinear coupling.  This is intentional — it keeps
%   the linear system unchanged and avoids per-call Newton iterations —
%   and is accurate when |omega| stays well below omega_max.  Configurations
%   that drive the motor close to stall should reduce AbsTol/RelTol or
%   increase max_omega_Hz to widen the safety margin.
%
% A run-token (`params.run_token`, set once per simulation in
% init_coupled_params) lets us reset the cache between independent runs in
% the same MATLAB session.  `clear ode_coupled_rhs` (issued in
% run_simulation.m before each ODE integration) is also performed as a
% defensive belt-and-braces reset.
persistent DPSI1_PREV LAST_RUN_TOKEN;
current_token = [];
if isfield(params, 'run_token')
    current_token = params.run_token;
end
if isempty(DPSI1_PREV) || ~isequal(LAST_RUN_TOKEN, current_token) ...
        || numel(DPSI1_PREV) ~= M_pili
    DPSI1_PREV    = cell(M_pili, 1);   % all empty -> first-call zero speed
    LAST_RUN_TOKEN = current_token;
end

%% Steps 1-7: Build coupled linear system (geometry + hydrodynamics + root BC)
has_bg_flow = ~strcmp(params.bg_flow.type, 'off');
if has_bg_flow
    [lin_sys, geom, hydro] = build_coupled_lin_sys(Z, params);
else
    [lin_sys, geom] = build_coupled_lin_sys(Z, params);
end

%% Step 8: Compute elastic restoring moments (RHS)
n_state  = 6 + 3 * sum(N_segs);
R_global = zeros(n_state, 1);
R_global(1:3) = -Eg * g_hat;

row_offset = 6;
for k = 1:M_pili
    N_k = N_segs(k);
    [kappa1, kappa2, kappa3] = compute_pili_elastic_moments( ...
        geom.d1_basis{k}, geom.d2_basis{k}, geom.d3_basis{k}, ...
        N_k, sigma, params.kb0{k}, rc{k});
    idx_k = row_offset + (1:3*N_k)';
    R_global(idx_k) = reshape([kappa1'; kappa2'; kappa3'], [], 1);

    %% Step 8b: Subtract active internal moments (per-pili, time-varying)
    if ~strcmp(params.am_cfg{k}.type, 'off')
        d1_k  = geom.d1_basis{k};
        d2_k  = geom.d2_basis{k};
        d3_k  = geom.d3_basis{k};
        [d1N_k, d2N_k, d3N_k] = compute_node_directors(d1_k, d2_k, d3_k, N_k);
        M_act = coupled_active_moment(t, k, N_k, d1_k, d2_k, d3_k, ...
                                      d1N_k, d2N_k, d3N_k, params.am_cfg{k}, rc{k});
        R_global(idx_k) = R_global(idx_k) - M_act(:);
    end

    %% Step 8d: Apply motor torque for free-torsion pili
    if strcmp(rc{k}, 'free_torsion') && isfield(params, 'motor_cfg') && ...
            ~strcmp(params.motor_cfg{k}.type, 'off')
        mc_k = params.motor_cfg{k};
        switch mc_k.type
            case 'constant'
                tau_motor = mc_k.torque_dimless;
                % Speed-limiting: read dpsi_1 from the most recent RHS
                % evaluation (persistent cache, see header note).  This
                % is an EXPLICIT LAGGED feedback — the cached value is
                % NOT guaranteed to be the previous accepted time step
                % under ode15s' adaptive multi-step behaviour, but the
                % linear stall-curve model tolerates this when the motor
                % operates well below omega_max.  On the first ODE call
                % DPSI1_PREV{k} is empty -> zero starting speed.
                omega_motor = 0;
                if k <= numel(DPSI1_PREV) && ~isempty(DPSI1_PREV{k})
                    omega_motor = DPSI1_PREV{k};
                end
                % Linear torque reduction as |omega| approaches max_omega.
                max_om = mc_k.max_omega_dimless;
                if abs(omega_motor) < max_om
                    tau_eff = tau_motor * (1 - abs(omega_motor) / max_om);
                else
                    tau_eff = 0;
                end
                % Apply motor torque to pili root (d3/torsional component)
                r3 = row_offset + 3;
                R_global(r3) = R_global(r3) + tau_eff;
                % Apply reaction torque to body (Newton's 3rd law).
                %
                % The motor exerts +tau along d3_lab on segment 1 of pilus k;
                % the body therefore receives -tau * d3_lab in the lab frame.
                % Rows 4-6 of R_global are the body lab-frame angular-momentum
                % equation (cf. B_global rows 4-6 in build_coupled_lin_sys.m,
                % which integrate r x f about the centre of mass), so the
                % physically correct reaction is the lab-frame d3 of segment 1.
                %
                % Earlier versions (<= v7.0.x) used J_euler(3,:)' here, which
                % is the kinematic mapping  omega -> dot_psi  (i.e. M^{-1} of
                % the ZYZ Euler-angle parameterisation), not the d3 vector.
                % That choice violated angular-momentum balance whenever
                % J_euler(3,:)' deviated from d3_lab (it always has a zero
                % z-component, while d3_lab generally does not).
                d3_lab_seg1 = geom.d3{k}(1, :)';
                R_global(4:6) = R_global(4:6) - tau_eff * d3_lab_seg1;
        end
    end

    row_offset = row_offset + 3 * N_k;
end

%% Step 8c: Subtract background flow contribution
if has_bg_flow
    u_bg = coupled_background_flow(geom.X_all, geom, params);
    back = apply_Ainv_u(u_bg, hydro);
    R_global = R_global - back;
end

%% Apply root BC RHS zeros
%  Indexes the rows of R_global that correspond to the (theta_1, phi_1,
%  psi_1) BC equations replaced in build_coupled_lin_sys Step 7.  For
%  free_torsion, only theta_1/phi_1 are kinematic constraints (RHS = 0);
%  psi_1's row is the natural moment-balance equation and must keep its
%  RHS (elastic kappa_3 + motor torque) intact.
%
%  The offset variable starts at 6 (skipping body V/Omega rows) and
%  advances by 3*N_k per pilus, mirroring the row_offset bookkeeping in
%  build_coupled_lin_sys Step 7.
row_offset_bc = 6;
for k = 1:M_pili
    N_k = N_segs(k);
    r1 = row_offset_bc + 1;
    R_global(r1)   = 0;      % theta_1 always zeroed
    R_global(r1+1) = 0;      % phi_1 always zeroed
    if strcmp(rc{k}, 'clamped')
        R_global(r1+2) = 0;  % psi_1 zeroed for clamped
    end
    % free_torsion: r1+2 preserved (elastic moment + motor torque)
    row_offset_bc = row_offset_bc + 3 * N_k;
end

%% Step 9: Solve
z_dot = lin_sys \ R_global;

% Store dpsi_1 for the lagged motor speed-limiting feedback on the next
% RHS call (see header note on adaptive-solver semantics).  Only
% free_torsion pili with a constant motor read the cache; other
% configurations still write it (cheap, no branching) but never consume.
for k = 1:M_pili
    N_k = N_segs(k);
    psi1_col = 6 + 3*sum(N_segs(1:k-1)) + 2*N_k + 1;
    DPSI1_PREV{k} = z_dot(psi1_col);
end

%% Step 10: Pack state rates
V     = z_dot(1:3);
Omega = z_dot(4:6);

dqdt       = quat_kinematics(geom.q, Omega);
q_curr     = geom.q;
dqdt       = dqdt - dot(dqdt, q_curr) * q_curr;

pili_rates = cell(params.M_pili, 1);
col_offset = 6;
for k = 1:params.M_pili
    N_k = N_segs(k);
    pili_rates{k}.dtheta = z_dot(col_offset + 1         : col_offset + N_k);
    pili_rates{k}.dphi   = z_dot(col_offset + N_k + 1   : col_offset + 2*N_k);
    pili_rates{k}.dpsi   = z_dot(col_offset + 2*N_k + 1 : col_offset + 3*N_k);
    col_offset = col_offset + 3 * N_k;
end

dZ = pack_state_rates(V, dqdt, pili_rates, N_segs);

end


%% ==================== Local Helper Functions ====================

function back = apply_Ainv_u(u_bg, hydro)
%APPLY_AINV_U  Compute B_global * Ainv * u_bg for background flow coupling
%
%  Prefers the reusable LU decomposition (hydro.dM_weighted) built by
%  build_coupled_lin_sys to skip a second full factorisation of the
%  mobility matrix; falls back to the raw matrix for backward
%  compatibility if the field is absent.
    n_t = hydro.n_trans;
    u_t = u_bg(1:n_t);
    u_w = u_bg(n_t+1:end);

    if isfield(hydro, 'dM_weighted') && ~isempty(hydro.dM_weighted)
        Ainv = hydro.dM_weighted;
    else
        Ainv = hydro.M_weighted;
    end

    if ~isempty(hydro.A_twist_total)
        Aw = hydro.A_twist_total * u_w;
        if ~isempty(hydro.M_rotlet)
            u_t = u_t - hydro.M_rotlet * Aw;
        end
        Au_t = Ainv \ u_t;
        Au   = [Au_t; Aw];
    else
        Au = Ainv \ u_t;
    end

    back = hydro.B_global * Au;
end
