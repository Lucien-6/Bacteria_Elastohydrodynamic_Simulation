function [dZ, flow_hydro] = ode_coupled_rhs(t, Z, params)
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
%      8a.    Subtract the Brownian thermal generalised body force held
%           constant over the current Brownian interval (if enabled)
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
%    [dZ, flow_hydro] = ode_coupled_rhs(t, Z, params)
%
%  The second output is optional.  It reuses the mobility matrices already
%  built for dZ and is intended for flow-field reconstruction.  Callers that
%  request only dZ do not copy those matrices.  flow_hydro contains
%  M_weighted, M_twist, M_rotlet, D_twist, C_twist, and C_twist_reverse.
%  It does not include LU factors.
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
%  Updated: 2026-07-13
%  Updated: 2026-09-04 — motor law: signed (asymmetric) linear torque-speed
%           relation in the rotor speed relative to the body, using the full
%           axial rate dpsi_1 + cos(theta_1) dphi_1 - Omega_body . d3_1
%           v8.4.0 — free_torsion psi_1 row carries only the boundary (motor)
%           torque; the interior twist moment is no longer double counted.
%           v8.5.0 — Step 8a: Brownian thermal generalised force realised
%           through the coupled RHS (params.brown_force_dimless).

%% Extract frequently used parameters
if ~isscalar(t) || ~isfinite(t) || any(~isfinite(Z(:)))
    error('numerical_diagnostics:NonFiniteState', ...
        'ODE time and state must be finite.');
end
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
% motor torque, the linear torque-speed law
%       tau_eff = tau_motor * (1 - omega_rel / omega_max),
%       omega_rel = sign(tau_motor) * [ (dpsi_1 + cos(theta_1) dphi_1)
%                                       - Omega_body . d3_1 ]
% is resolved within a single RHS evaluation using one extra back-substitution
% per motor against the already-factored lin_sys (multi-column right-hand
% side, so lin_sys is factorised only once).  omega_rel is the rotor speed
% RELATIVE to the cell body, measured positive in the motor's driving
% direction: it uses the full axial angular velocity of segment 1
% (dpsi_1 + cos(theta_1) dphi_1, not dpsi_1 alone) minus the body's spin
% about the same axis.  The law is intentionally asymmetric (no |omega|):
% rotation against the motor (omega_rel < 0) raises the torque above the
% stall value and rotation faster than omega_max makes the motor a brake,
% following the measured linear torque-speed relation of the bacterial
% flagellar motor (Berg & Turner 1993; Berg 2003).  This keeps the RHS a
% PURE function of (t, Z) -- a prerequisite for ode15s' internal Jacobian
% and error control -- and is exact even near stall.
% The former persistent DPSI1_PREV / run-token machinery has been removed.

%% Steps 1-7: Build coupled linear system (geometry + hydrodynamics + root BC)
has_bg_flow = ~strcmp(params.bg_flow.type, 'off');
want_flow_hydro = nargout >= 2;
if has_bg_flow || want_flow_hydro
    [lin_sys, geom, hydro] = build_coupled_lin_sys(Z, params);
else
    [lin_sys, geom] = build_coupled_lin_sys(Z, params);
end

%% Step 8: Compute elastic restoring moments (RHS)
n_state  = 6 + 3 * sum(N_segs);
R_global = zeros(n_state, 1);
R_global(1:3) = -Eg * g_hat;

%% Step 8a: Brownian thermal generalised force on the body (v8.5.0)
%  Body-only Brownian dynamics is realised as a constant generalised
%  force/torque held over one Brownian interval (see
%  sample_brownian_body_force.m).  It enters with the same sign convention
%  as gravity: an applied dimensionless force f enters R as -f.  The field
%  is zero/empty for deterministic runs and between intervals, so the RHS
%  remains a pure function of (t, Z) within any single ode15s segment.
if isfield(params, 'brown_force_dimless') && ~isempty(params.brown_force_dimless)
    R_global(1:6) = R_global(1:6) - params.brown_force_dimless(:);
end

% Collectors for motorised pili (free_torsion + constant motor).  The motor
% torque is NOT added to R_global in this loop; R_global is assembled as the
% motor-free baseline RHS and the stall-curve torque is resolved exactly in
% Step 9 (defect-2 remediation).
motor_rows    = [];   % moment-balance row index (segment-1 psi_1 row) per motor
motor_obs     = zeros(n_state, 0);  % omega_rel = motor_obs(:,i)' * z_dot per motor
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
                %
                % Rotor speed observation vector: omega_rel = obs' * z_dot.
                % Axial angular velocity of segment 1 about its own d3 is
                %   omega_1 = dpsi_1 + cos(theta_1) * dphi_1
                % (Walker Euler-angle kinematics, same expression used for
                % the twist-rate rows in build_coupled_lin_sys).  The motor
                % sits between rotor (segment 1) and cell body, so the speed
                % that enters the torque-speed law is the RELATIVE spin
                %   omega_rel = omega_1 - Omega_body . d3_1(lab).
                % d3 is stored in the pilus basis; rot_k' maps it to the lab.
                theta_1 = geom.pili_angles{k}.theta(1);
                d3_1_lab = params.rot{k}' * geom.d3_basis{k}(1, :)';
                obs = zeros(n_state, 1);
                obs(4:6)                    = -d3_1_lab;
                obs(col_state + N_k + 1)    = cos(theta_1);   % dphi_1
                obs(col_state + 2*N_k + 1)  = 1;              % dpsi_1
                motor_rows(end+1)    = row_offset + 3;          %#ok<AGROW>
                motor_obs(:, end+1)  = obs;                     %#ok<AGROW>
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
%  psi_1's row is the whole-rod d3 moment balance about the root.  Its RHS
%  is the body->rod boundary torque only: compute_pili_elastic_moments and
%  coupled_active_moment already return zero for the root row (v8.4.0), so
%  after this loop the row carries the background-flow term (if any) and
%  receives the motor torque in Step 9.  A free bearing without a motor
%  therefore enforces zero net axial hydrodynamic torque on the filament.
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
    % free_torsion: r1+2 preserved (boundary torque row; motor added in Step 9)
    row_offset_bc = row_offset_bc + 3 * N_k;
end

%% Step 9: Solve (with exact self-consistent motor speed-limiting)
%  The assembled system is linear in the motor torques, so
%      z_dot = zdot_base + sum_i tau_i * S(:, i),
%  where zdot_base solves the motor-free baseline RHS and S(:, i) is the
%  response to a unit torque on motor i.  Both come from a SINGLE
%  factorisation of lin_sys via a multi-column right-hand side.  The
%  linear torque-speed law
%      tau_i = tau_motor_i * (1 - omega_rel,i / omega_max_i),
%      omega_rel,i = sign(tau_motor_i) * motor_obs(:,i)' * z_dot
%  is then solved self-consistently in the reduced motor space.  Because
%  the law is linear in omega_rel, the reduced problem is a small linear
%  system with a unique solution (see solve_motor_self_consistent).  This
%  replaces the former lagged persistent cache and makes the RHS a pure
%  function of (t, Z) (defect-2 remediation).
n_motor = numel(motor_rows);
if any(~isfinite(R_global))
    error('numerical_diagnostics:NonFiniteRhs', ...
        'RHS vector R_global must be finite.');
end
if n_motor == 0
    z_dot = lin_sys \ R_global;
    check_linear_solve(lin_sys, z_dot, R_global, ...
        'ode_coupled_rhs base solve');
else
    % Sign convention of the motor torque.  In this system R_global holds
    % the NEGATIVE of the applied generalised force (gravity enters as
    % R(1:3) = -Eg*g_hat, and a positive elastic kappa_3 in a psi row
    % untwists the rod, i.e. rotates it negatively about d3).  A motor torque
    % +tau that drives the rotor in the positive right-handed sense about
    % the root tangent d3 (counter-clockwise when viewed from the filament
    % tip looking toward the cell body -- the E. coli "run" sense for a
    % left-handed filament) must therefore enter its psi_1 row as -tau.
    E_motor = zeros(n_state, n_motor);
    for i = 1:n_motor
        E_motor(motor_rows(i), i) = -1;
    end

    sol_aug   = lin_sys \ [R_global, E_motor];
    check_linear_solve(lin_sys, sol_aug, [R_global, E_motor], ...
        'ode_coupled_rhs augmented motor solve');
    zdot_base = sol_aug(:, 1);
    S_motor   = sol_aug(:, 2:end);

    % Reduced motor problem in the RELATIVE rotor speed (lab-independent):
    %   omega_rel = a_omega + B_omega * tau   (before the sign convention
    %   applied inside solve_motor_self_consistent)
    a_omega = motor_obs' * zdot_base;       % [n_motor x 1] baseline omega_rel
    B_omega = motor_obs' * S_motor;         % [n_motor x n_motor] sensitivity

    tau_eff = solve_motor_self_consistent(a_omega(:), B_omega, ...
        motor_tau(:), motor_ommax(:), get_rhs_logf(params));

    z_dot = zdot_base + S_motor * tau_eff;
    check_linear_solve(lin_sys, z_dot, R_global + E_motor * tau_eff, ...
        'ode_coupled_rhs motor-coupled solution');
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

if want_flow_hydro
    flow_hydro = flow_hydro_for_coefficients(hydro);
else
    flow_hydro = [];
end

end


%% ==================== Local Helper Functions ====================

function flow_hydro = flow_hydro_for_coefficients(hydro)
%FLOW_HYDRO_FOR_COEFFICIENTS  Mobility blocks for the flow-field force solve.
%  LU factors and the other coupled-system intermediates stay local to the
%  right-hand-side evaluation.
flow_hydro = struct();
flow_hydro.M_weighted = hydro.M_weighted;
flow_hydro.M_twist = hydro.M_twist;
flow_hydro.M_rotlet = hydro.M_rotlet;
flow_hydro.D_twist = hydro.D_twist;
flow_hydro.C_twist = hydro.C_twist;
flow_hydro.C_twist_reverse = hydro.C_twist_reverse;
end

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
        check_linear_solve(hydro.M_weighted, Minv_u, u_t, ...
            'ode_coupled_rhs background translational solve');
        rhs_w  = u_w - hydro.C_twist_reverse * Minv_u;

        if isfield(hydro, 'dS_twist') && ~isempty(hydro.dS_twist)
            m_twist = hydro.dS_twist \ rhs_w;
        else
            m_twist = hydro.S_twist \ rhs_w;
        end
        check_linear_solve(hydro.S_twist, m_twist, rhs_w, ...
            'ode_coupled_rhs background Schur solve');

        if isfield(hydro, 'Minv_C_twist') && ~isempty(hydro.Minv_C_twist)
            Au_t = Minv_u - hydro.Minv_C_twist * m_twist;
        else
            translational_rhs = u_t - hydro.C_twist * m_twist;
            Au_t = Ainv \ translational_rhs;
            check_linear_solve(hydro.M_weighted, Au_t, ...
                translational_rhs, ...
                'ode_coupled_rhs background coupled translational solve');
        end
        Au = [Au_t; hydro.D_twist * m_twist];
    else
        Au = Ainv \ u_t;
        check_linear_solve(hydro.M_weighted, Au, u_t, ...
            'ode_coupled_rhs background solve');
    end

    back = hydro.B_global * Au;
end


function logf_h = get_rhs_logf(params)
%GET_RHS_LOGF  Retrieve logf handle from params (empty if unavailable)
    if isstruct(params) && isfield(params, 'logf') && ~isempty(params.logf)
        logf_h = params.logf;
    else
        logf_h = [];
    end
end
