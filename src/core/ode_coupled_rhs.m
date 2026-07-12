function dZ = ode_coupled_rhs(t, Z, params)
%ODE_COUPLED_RHS  Core ODE right-hand side for coupled body-pili system
%
%  Description:
%    Computes dZ/dt for the coupled bacterial body + elastic pili system.
%    The ODE state has length 7 + 3*sum(N_segments):
%      Z = [X_cm(3); q(4); theta; phi; psi].
%    The internal generalized-rate vector has length 6 + 3*sum(N_segments):
%      z_dot = [V(3); Omega(3); dtheta; dphi; dpsi],
%    which is also the row and column dimension of lin_sys.
%
%    Evaluation steps:
%      1-4, 6-7. Build the coupled linear system (delegated to
%           build_coupled_lin_sys.m, which handles Step 1: geometry
%           reconstruction, Step 2: global mobility/resistance, Step 3-4:
%           Q/B matrix assembly, Step 6: linear system assembly,
%           Step 7: root boundary conditions)
%      5.     (Not applicable — elastic restoring moments are a RHS term,
%           computed here as Step 8, not part of the LHS linear system)
%      8.     Compute elastic restoring moments for each pilus
%      8b.    Subtract per-pili active internal moments (if enabled)
%      8d.    Record free-torsion motor rows and torque-law data
%      8c.    Subtract the background-flow contribution (if enabled);
%           Step 8d intentionally precedes Step 8c in the implementation
%      BC RHS. Set theta_1/phi_1 rows to zero for every pilus and also
%           psi_1 for clamped roots; preserve the free-torsion psi_1
%           moment-balance RHS
%      9.     Solve lin_sys * z_dot = R_global.  Motor speed limiting is
%           resolved exactly and self-consistently within this RHS call
%           from the motor-free solution and unit-torque sensitivities;
%           no lagged motor state is used
%     10.     Pack state rates into dZ
%
%  Syntax:
%    dZ = ode_coupled_rhs(t, Z, params)
%
%  Inputs:
%    t      - Current dimensionless time
%    Z      - [(7 + 3*sum(N_segments)) x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    dZ - [(7 + 3*sum(N_segments)) x 1] Time derivative of state vector
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-07-12

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

% Motor speed-limiting (free_torsion + constant motor) is handled by an
% EXACT in-call self-consistent solve (see Step 8d / Step 9 below), NOT by a
% lagged persistent cache.  Because the assembled system is linear in the
% motor torque, the stall-curve coupling
%       tau_eff = tau_motor * (1 - |omega| / omega_max),   omega = dpsi_1
% is resolved within a single RHS evaluation using one extra back-substitution
% per motor against the already-factored lin_sys (multi-column right-hand
% side, so lin_sys is factorised only once).  This keeps the RHS a PURE
% function of (t, Z) -- a prerequisite for ode15s' internal Jacobian and
% error control -- and is exact even near stall (defect-2 remediation).
% The former persistent DPSI1_PREV / run-token machinery has been removed.

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

% Collectors for motorised pili (free_torsion + constant motor).  The motor
% torque is NOT added to R_global in this loop; R_global is assembled as the
% motor-free baseline RHS and the stall-curve torque is resolved exactly in
% Step 9 (defect-2 remediation).
motor_rows    = [];   % moment-balance row index (segment-1 psi_1 row) per motor
motor_psicols = [];   % state-rate column index of dpsi_1 per motor
motor_tau     = [];   % constant motor torque (dimensionless)
motor_ommax   = [];   % stall speed (dimensionless)

row_offset = 6;
col_state  = 6;       % running offset to locate each pilus's dpsi_1 column
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

    %% Step 8d: Record motor torque for free-torsion pili (resolved in Step 9)
    if strcmp(rc{k}, 'free_torsion') && isfield(params, 'motor_cfg') && ...
            ~strcmp(params.motor_cfg{k}.type, 'off')
        mc_k = params.motor_cfg{k};
        switch mc_k.type
            case 'constant'
                % The motor torque enters the segment-1 d3 (psi_1) moment-
                % balance row, which under free_torsion is the natural
                % equation driving dpsi_1.  We only RECORD the row, the
                % corresponding dpsi_1 column, and the stall-curve parameters
                % here; the self-consistent torque is solved EXACTLY in Step 9
                % (the system is linear in tau, so this needs no iteration on
                % the full linear system).
                %
                % NOTE: NO term is added to R_global(4:6).  Rows 4-6 are the
                % FULL SYSTEM angular-momentum balance.  The motor torque
                % between body and pilus is INTERNAL: the action (+tau) on the
                % pilus and the reaction (-tau) on the body cancel when summed
                % over the whole body+pilus system (Newton's third law), so no
                % contribution belongs in the system moment equation.  The
                % body's reaction motion is propagated implicitly through the
                % global mobility solve and the root kinematic constraints
                % (theta_1, phi_1 locked to body Omega).
                motor_rows(end+1)    = row_offset + 3;          %#ok<AGROW>
                motor_psicols(end+1) = col_state + 2*N_k + 1;   %#ok<AGROW>
                motor_tau(end+1)     = mc_k.torque_dimless;     %#ok<AGROW>
                motor_ommax(end+1)   = mc_k.max_omega_dimless;  %#ok<AGROW>
            otherwise
                error('ode_coupled_rhs:UnknownMotorType', ...
                    ['Pili %d: motor_cfg.type = ''%s'' is not implemented in RHS. ', ...
                     'Supported: constant.'], k, mc_k.type);
        end
    end

    row_offset = row_offset + 3 * N_k;
    col_state  = col_state  + 3 * N_k;
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

%% Step 9: Solve (with exact self-consistent motor speed-limiting)
%  The assembled system is linear in the motor torques, so
%      z_dot = zdot_base + sum_i tau_i * S(:, i),
%  where zdot_base solves the motor-free baseline RHS and S(:, i) is the
%  response to a unit torque on motor i.  Both come from a SINGLE
%  factorisation of lin_sys via a multi-column right-hand side.  The
%  stall-curve law
%      tau_i = tau_motor_i * (1 - |omega_i| / omega_max_i),  omega_i = dpsi_1,i
%  is then solved self-consistently in the reduced motor space (closed form
%  for a single motor; damped fixed point for several).  This replaces the
%  former lagged persistent cache and makes the RHS a pure function of
%  (t, Z) (defect-2 remediation).
n_motor = numel(motor_rows);
if n_motor == 0
    z_dot = lin_sys \ R_global;
else
    E_motor = zeros(n_state, n_motor);
    for i = 1:n_motor
        E_motor(motor_rows(i), i) = 1;
    end

    sol_aug   = lin_sys \ [R_global, E_motor];
    zdot_base = sol_aug(:, 1);
    S_motor   = sol_aug(:, 2:end);

    a_omega = zdot_base(motor_psicols);     % [n_motor x 1] baseline dpsi_1
    B_omega = S_motor(motor_psicols, :);    % [n_motor x n_motor] sensitivity

    tau_eff = solve_motor_self_consistent(a_omega(:), B_omega, ...
        motor_tau(:), motor_ommax(:), get_rhs_logf(params));

    z_dot = zdot_base + S_motor * tau_eff;
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

    if isfield(hydro, 'M_twist') && ~isempty(hydro.M_twist)
        Minv_u = Ainv \ u_t;
        rhs_w  = u_w - hydro.C_twist' * Minv_u;

        if isfield(hydro, 'dS_twist') && ~isempty(hydro.dS_twist)
            m_twist = hydro.dS_twist \ rhs_w;
        else
            m_twist = hydro.S_twist \ rhs_w;
        end

        if isfield(hydro, 'Minv_C_twist') && ~isempty(hydro.Minv_C_twist)
            Au_t = Minv_u - hydro.Minv_C_twist * m_twist;
        else
            Au_t = Ainv \ (u_t - hydro.C_twist * m_twist);
        end
        Au = [Au_t; hydro.D_twist * m_twist];
    else
        Au = Ainv \ u_t;
    end

    back = hydro.B_global * Au;
end


function tau = solve_motor_self_consistent(a, B, tau_motor, om_max, logf_h)
%SOLVE_MOTOR_SELF_CONSISTENT  Exact stall-curve torque(s) for free-torsion motors
%
%  Solves the reduced self-consistency problem
%      omega = a + B * tau,
%      tau_i = motor_law(omega_i),
%  where  motor_law(om) = tau_motor*(1 - |om|/om_max)  if |om| < om_max,
%  else 0.  Here a (= baseline dpsi_1) and B (= sensitivity of dpsi_1 to a
%  unit motor torque) are precomputed from a single factorisation of the
%  coupled linear system, so this routine performs NO further linear solves.
%
%  For a single motor the solution is closed form; for several motors a
%  damped fixed-point iteration is used (the map is a contraction whenever
%  the feedback gain |B*tau_motor/om_max| < 1, i.e. away from the critical
%  stall-feedback limit).
    n = numel(a);
    if n == 1
        tau = motor_law_root_scalar(a, B(1,1), tau_motor(1), om_max(1));
        return;
    end

    tau   = zeros(n, 1);
    relax = 0.7;
    for it = 1:500 %#ok<NASGU>
        omega   = a + B * tau;
        tau_t   = motor_law_vec(omega, tau_motor, om_max);
        tau_new = (1 - relax) * tau + relax * tau_t;
        if norm(tau_new - tau, inf) <= 1e-13 * (1 + norm(tau_new, inf))
            tau = tau_new;
            return;
        end
        tau = tau_new;
    end
    if ~isempty(logf_h)
        logf_h('[ODE] WARNING: motor self-consistent solve did not fully converge (n=%d motors).\n', n);
    end
end


function tau = motor_law_root_scalar(a, b, tm, om_max)
%MOTOR_LAW_ROOT_SCALAR  Closed-form self-consistent torque for one motor
%  f(omega) = a + b*tau(omega) - omega is continuous, piecewise linear and
%  strictly decreasing (for subcritical feedback), so it has a unique root.
%  Exactly one of the three branches below satisfies its own consistency
%  test; the saturated branch (tau = 0) is the continuous fallback.
    if ~(om_max > 0)
        tau = tm;   % no speed limiting if om_max is non-positive
        return;
    end

    % Active branch, omega >= 0:  tau = tm*(1 - omega/om_max)
    om_p = (a + b*tm) / (1 + b*tm/om_max);
    if om_p >= 0 && abs(om_p) < om_max
        tau = tm * (1 - abs(om_p)/om_max);
        return;
    end

    % Active branch, omega < 0:   tau = tm*(1 + omega/om_max)
    om_n = (a + b*tm) / (1 - b*tm/om_max);
    if om_n < 0 && abs(om_n) < om_max
        tau = tm * (1 - abs(om_n)/om_max);
        return;
    end

    % Saturated branch: |omega| >= om_max  =>  tau = 0
    tau = 0;
end


function tau = motor_law_vec(omega, tau_motor, om_max)
%MOTOR_LAW_VEC  Elementwise stall-curve torque from current angular speeds
    n   = numel(omega);
    tau = zeros(n, 1);
    for i = 1:n
        if om_max(i) > 0 && abs(omega(i)) < om_max(i)
            tau(i) = tau_motor(i) * (1 - abs(omega(i)) / om_max(i));
        else
            tau(i) = 0;
        end
    end
end


function logf_h = get_rhs_logf(params)
%GET_RHS_LOGF  Retrieve logf handle from params (empty if unavailable)
    if isstruct(params) && isfield(params, 'logf') && ~isempty(params.logf)
        logf_h = params.logf;
    else
        logf_h = [];
    end
end
