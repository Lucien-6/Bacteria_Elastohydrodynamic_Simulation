function J = ode_coupled_jacobian(t, Z, params)
%ODE_COUPLED_JACOBIAN  Fast frozen-hydrodynamics Jacobian for ode15s
%
%  Description:
%    Computes an accurate Newton-iteration Jacobian J = d(dZ/dt)/dZ for the
%    coupled body-pili elastohydrodynamic ODE system:
%      lin_sys(Z) * z_dot = R(t, Z).
%
%    In stiff filament elastohydrodynamics, the large O(N^4) Jacobian
%    eigenvalues arise from the Kirchhoff elastic restoring moments (and
%    local active moments) dR/dZ, whereas the hydrodynamic resistance
%    operator lin_sys(Z) varies smoothly at O(1).  Evaluating lin_sys(Z)
%    ONCE at the current state and solving lin_sys \ (dR/dZ) via a single
%    LU factorisation with 3-colour banded perturbations of R reduces the
%    Jacobian cost from O(n_ode) full hydrodynamic assemblies to O(1) —
%    typically 50x-400x faster per Jacobian evaluation — without altering
%    the ODE right-hand side or ode15s local error control.
%
%    Exact terms included:
%      - Full elastic + active moment derivative dR/d(theta, phi, psi)
%        (root segment handled individually for root-frame-coupled active
%        moments; interior segments 3-coloured for local stencils, or
%        column-wise for non-local traveling_wave moments).
%      - Exact linearised motor torque-speed feedback for free_torsion
%        motorised pili.
%      - Exact quaternion kinematics block d(dq/dt)/d(q, Omega).
%      - Exact full finite-difference Jacobian for bare-body cases (M_pili == 0).
%
%  Syntax:
%    J = ode_coupled_jacobian(t, Z, params)
%
%  Inputs:
%    t      - Current dimensionless time
%    Z      - [(7 + 3*sum(N_segments)) x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    J      - [n_ode x n_ode] Jacobian matrix for ode15s ('Jacobian' option)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

n_ode  = numel(Z);
M_pili = params.M_pili;
N_segs = params.N_segments;

% ---- Bare body (M_pili == 0): 7-column exact FD takes < 1 ms ----
if M_pili == 0
    f0 = ode_coupled_rhs(t, Z, params);
    J  = zeros(n_ode, n_ode);
    h  = 1e-7;
    has_bg = ~strcmp(params.bg_flow.type, 'off');
    j_start = 4;
    if has_bg
        j_start = 1;
    end
    for j = j_start:n_ode
        dy = h * max(abs(Z(j)), 1);
        Zj = Z;
        Zj(j) = Z(j) + dy;
        J(:, j) = (ode_coupled_rhs(t, Zj, params) - f0) / dy;
    end
    return;
end

n_state = 6 + 3 * sum(N_segs);

%% Step 1: Assemble coupled linear system ONCE at current Z
has_bg_flow      = ~strcmp(params.bg_flow.type, 'off');
has_pili_thermal = isfield(params, 'brown_pili_force') && ~isempty(params.brown_pili_force);
if has_bg_flow || has_pili_thermal
    [lin_sys, geom, hydro] = build_coupled_lin_sys(Z, params, 'rhs_light');
else
    [lin_sys, geom] = build_coupled_lin_sys(Z, params);
    hydro = [];
end

%% Step 2: Baseline RHS R0 (full, for z_dot0) and R0_local (elastic+active only, for dR/dZ)
[R0, R0_local, motor_info, rc] = assemble_rhs_vector( ...
    t, geom, hydro, params, has_bg_flow, has_pili_thermal);

%% Step 3: Factor lin_sys once and compute baseline z_dot0 (and motor sensitivities)
dL = decomposition(lin_sys, 'lu', 'CheckCondition', false);

if motor_info.n_motor == 0
    z_dot0  = dL \ R0;
    S_motor = [];
    B_omega = [];
else
    E_motor = zeros(n_state, motor_info.n_motor);
    for im = 1:motor_info.n_motor
        E_motor(motor_info.rows(im), im) = -1;
    end
    sol_aug   = dL \ [R0, E_motor];
    zdot_base = sol_aug(:, 1);
    S_motor   = sol_aug(:, 2:end);
    a_omega   = motor_info.obs' * zdot_base;
    B_omega   = motor_info.obs' * S_motor;
    tau_eff   = solve_motor_self_consistent(a_omega(:), B_omega, ...
        motor_info.tau(:), motor_info.ommax(:), []);
    z_dot0    = zdot_base + S_motor * tau_eff;
end

%% Step 4: Sparse / 3-coloured dR/d(angles) (zero hydrodynamic recomputation)
dR_dZ = zeros(n_state, n_ode);
h     = 1e-7;
row_off = 6;
z_off   = 7;

for k = 1:M_pili
    N_k    = N_segs(k);
    rc_k   = rc{k};
    th0    = geom.pili_angles{k}.theta;
    ph0    = geom.pili_angles{k}.phi;
    ps0    = geom.pili_angles{k}.psi;
    r_idx  = row_off + (1 : 3*N_k);
    rk0    = R0_local(r_idx);
    has_am = ~strcmp(params.am_cfg{k}.type, 'off');
    % traveling_wave integrates from tip to root (non-local cumulative sum),
    % so perturb each segment individually.  Custom moments default to the
    % exact j=1 + 3-colour local path (e.g. Chlamydomonas tabulated drive);
    % set am_cfg.is_nonlocal = true for non-local custom moment laws.
    is_nonlocal_am = has_am && (strcmp(params.am_cfg{k}.type, 'traveling_wave') || ...
        (isfield(params.am_cfg{k}, 'is_nonlocal') && logical(params.am_cfg{k}.is_nonlocal)));

    for comp = 1:3
        % Root segment (j = 1) individually to capture any root-frame coupling
        th = th0;  ph = ph0;  ps = ps0;
        if comp == 1
            dy = h * max(abs(th(1)), 1);  th(1) = th(1) + dy;  col = z_off + 1;
        elseif comp == 2
            dy = h * max(abs(ph(1)), 1);  ph(1) = ph(1) + dy;  col = z_off + N_k + 1;
        else
            dy = h * max(abs(ps(1)), 1);  ps(1) = ps(1) + dy;  col = z_off + 2*N_k + 1;
        end
        rk = eval_pilus_rhs(t, k, th, ph, ps, N_k, params, rc_k, has_am);
        dR_dZ(r_idx, col) = (rk - rk0) / dy;

        if ~is_nonlocal_am && N_k >= 2
            % Segment j (j >= 2) only enters nodes j and j+1 in Kirchhoff
            % curvature and node-local active moments -> stride-3 colouring
            % is exact in 3 calls regardless of N_k.
            for color = 0:2
                js = (2 + color) : 3 : N_k;
                if isempty(js)
                    continue;
                end
                th = th0;  ph = ph0;  ps = ps0;
                if comp == 1
                    dys = h * max(abs(th(js)), 1);  th(js) = th(js) + dys;  cols = z_off + js;
                elseif comp == 2
                    dys = h * max(abs(ph(js)), 1);  ph(js) = ph(js) + dys;  cols = z_off + N_k + js;
                else
                    dys = h * max(abs(ps(js)), 1);  ps(js) = ps(js) + dys;  cols = z_off + 2*N_k + js;
                end
                rk = eval_pilus_rhs(t, k, th, ph, ps, N_k, params, rc_k, has_am);
                dr = rk - rk0;
                for ii = 1:numel(js)
                    jj   = js(ii);
                    n_lo = max(2, jj);
                    n_hi = min(N_k, jj + 1);
                    rr   = (3*(n_lo - 1) + 1) : (3 * n_hi);
                    dR_dZ(row_off + rr, cols(ii)) = dr(rr) / dys(ii);
                end
            end
        else
            for j = 2:N_k
                th = th0;  ph = ph0;  ps = ps0;
                if comp == 1
                    dy = h * max(abs(th(j)), 1);  th(j) = th(j) + dy;  col = z_off + j;
                elseif comp == 2
                    dy = h * max(abs(ph(j)), 1);  ph(j) = ph(j) + dy;  col = z_off + N_k + j;
                else
                    dy = h * max(abs(ps(j)), 1);  ps(j) = ps(j) + dy;  col = z_off + 2*N_k + j;
                end
                rk = eval_pilus_rhs(t, k, th, ph, ps, N_k, params, rc_k, has_am);
                dR_dZ(r_idx, col) = (rk - rk0) / dy;
            end
        end
    end

    row_off = row_off + 3 * N_k;
    z_off   = z_off   + 3 * N_k;
end

%% Step 5: Single multi-RHS back-substitution for dz_dot / d(angles)
dzdot_dZ = zeros(n_state, n_ode);
dzdot_dZ(:, 8:end) = dL \ dR_dZ(:, 8:end);

% Linearised motor torque-speed feedback (matching solve_motor_self_consistent):
%   tau_i = tau_0,i * (1 - s_i * obs_i' * z_dot / om_max,i)  for limited motors
%   => (I + diag(slope) * B_omega) * dtau = -diag(slope) * obs' * dzdot_base
if motor_info.n_motor > 0
    om_max  = motor_info.ommax(:);
    limited = isfinite(om_max) & (om_max > 0);
    if any(limited)
        c_vec = zeros(motor_info.n_motor, 1);
        c_vec(limited) = abs(motor_info.tau(limited)) ./ om_max(limited);
        M_red = eye(motor_info.n_motor) + c_vec .* B_omega;
        dtau  = -(M_red \ (c_vec .* (motor_info.obs' * dzdot_dZ(:, 8:end))));
        dzdot_dZ(:, 8:end) = dzdot_dZ(:, 8:end) + S_motor * dtau;
    end
end

%% Step 6: Map generalised-rate Jacobian to ODE state Jacobian (incl. quaternion kinematics)
J = zeros(n_ode, n_ode);
J(1:3, :)   = dzdot_dZ(1:3, :);
J(8:end, :) = dzdot_dZ(7:end, :);

% Exact Galilean translation Jacobian for linear shear background flow:
%   u_bg(X + dX) = u_bg(X) + (U_max / H) * (s_hat' * dX) * d_hat,
%   a uniform shift adds directly to V_cm without bending or spinning the cell.
if has_bg_flow && strcmp(params.bg_flow.type, 'shear_z')
    shear_rate = params.bg_flow.U_max_dimless / params.bg_flow.shear_height_dimless;
    J(1:3, 1:3) = shear_rate * (params.bg_flow.direction(:) * params.bg_flow.shear_axis(:)');
end

q   = geom.q;
Om0 = z_dot0(4:6);
G_q = 0.5 * [-q(2), -q(3), -q(4); ...
              q(1),  q(4), -q(3); ...
             -q(4),  q(1),  q(2); ...
              q(3), -q(2),  q(1)];
P_q = eye(4) - q * q';
J(4:7, 8:end) = (P_q * G_q) * dzdot_dZ(4:6, 8:end);

% Quaternion-to-quaternion block d(dq/dt)/dq at frozen Omega0
dq0 = quat_kinematics(q, Om0);
dq0 = dq0 - dot(dq0, q) * q;
for jq = 1:4
    q_p = Z(4:7);
    dy  = h * max(abs(q_p(jq)), 1);
    q_p(jq) = q_p(jq) + dy;
    qn  = q_p / norm(q_p);
    dqp = quat_kinematics(qn, Om0);
    dqp = dqp - dot(dqp, qn) * qn;
    J(4:7, 3 + jq) = (dqp - dq0) / dy;
end

end


%% ==================== Local Helper Functions ====================

function rk = eval_pilus_rhs(t, k, th, ph, ps, N_k, params, rc_k, has_am)
%EVAL_PILUS_RHS  Elastic + active moment rows for pilus k (no hydrodynamics)
    [d1b, d2b, d3b] = pili_directors(th, ph, ps);
    [d1N, d2N, d3N] = compute_node_directors(d1b, d2b, d3b, N_k);
    [k1, k2, k3]    = compute_pili_elastic_moments( ...
        d1b, d2b, d3b, N_k, params.sigma, params.kb0{k}, rc_k, d1N, d2N, d3N);
    rk = reshape([k1'; k2'; k3'], [], 1);
    if has_am
        M_act = coupled_active_moment(t, k, N_k, d1b, d2b, d3b, ...
                                      d1N, d2N, d3N, params.am_cfg{k}, rc_k);
        rk = rk - M_act(:);
    end
    rk(1) = 0;
    rk(2) = 0;
    if strcmp(rc_k, 'clamped')
        rk(3) = 0;
    end
end


function [R_global, R_local, mi, rc] = assemble_rhs_vector(t, geom, hydro, params, has_bg_flow, has_pili_thermal)
%ASSEMBLE_RHS_VECTOR  Baseline R_global (full) and R_local (elastic+active only)
    N_segs = params.N_segments;
    M_pili = params.M_pili;
    n_state = 6 + 3 * sum(N_segs);

    R_global = zeros(n_state, 1);
    R_global(1:3) = -params.Eg * params.g_hat(:);
    if isfield(params, 'brown_force_dimless') && ~isempty(params.brown_force_dimless)
        R_global(1:6) = R_global(1:6) - params.brown_force_dimless(:);
    end

    rc = cell(M_pili, 1);
    mi.n_motor = 0;
    mi.rows    = zeros(M_pili, 1);
    mi.obs     = zeros(n_state, M_pili);
    mi.tau     = zeros(M_pili, 1);
    mi.ommax   = zeros(M_pili, 1);

    row_off = 6;
    col_st  = 6;
    has_cached_dN = isfield(geom, 'd1N_basis') && numel(geom.d1N_basis) >= M_pili;
    for k = 1:M_pili
        N_k = N_segs(k);
        if isfield(params, 'root_coupling') && numel(params.root_coupling) >= k
            rc{k} = params.root_coupling{k};
        else
            rc{k} = 'clamped';
        end
        d1_k = geom.d1_basis{k};
        d2_k = geom.d2_basis{k};
        d3_k = geom.d3_basis{k};
        if has_cached_dN && ~isempty(geom.d1N_basis{k})
            d1N_k = geom.d1N_basis{k};
            d2N_k = geom.d2N_basis{k};
            d3N_k = geom.d3N_basis{k};
        else
            [d1N_k, d2N_k, d3N_k] = compute_node_directors(d1_k, d2_k, d3_k, N_k);
        end
        [k1, k2, k3] = compute_pili_elastic_moments( ...
            d1_k, d2_k, d3_k, N_k, params.sigma, params.kb0{k}, rc{k}, d1N_k, d2N_k, d3N_k);
        idx_k = row_off + (1 : 3*N_k)';
        R_global(idx_k) = reshape([k1'; k2'; k3'], [], 1);

        if ~strcmp(params.am_cfg{k}.type, 'off')
            M_act = coupled_active_moment(t, k, N_k, d1_k, d2_k, d3_k, ...
                                          d1N_k, d2N_k, d3N_k, params.am_cfg{k}, rc{k});
            R_global(idx_k) = R_global(idx_k) - M_act(:);
        end

        if strcmp(rc{k}, 'free_torsion') && isfield(params, 'motor_cfg') && ...
                ~strcmp(params.motor_cfg{k}.type, 'off')
            mc = params.motor_cfg{k};
            switch mc.type
                case 'constant'
                    th1 = geom.pili_angles{k}.theta(1);
                    d3_1_lab = params.rot{k}' * d3_k(1, :)';
                    obs = zeros(n_state, 1);
                    obs(4:6)                = -d3_1_lab;
                    obs(col_st + N_k + 1)   = cos(th1);
                    obs(col_st + 2*N_k + 1) = 1;
                    mi.n_motor = mi.n_motor + 1;
                    mi.rows(mi.n_motor)    = row_off + 3;
                    mi.obs(:, mi.n_motor)  = obs;
                    mi.tau(mi.n_motor)     = mc.torque_dimless;
                    mi.ommax(mi.n_motor)   = mc.max_omega_dimless;
                otherwise
                    error('ode_coupled_jacobian:UnknownMotorType', ...
                        ['Pili %d: motor_cfg.type = ''%s'' is not implemented in Jacobian. ', ...
                         'Supported: constant.'], k, mc.type);
            end
        end

        row_off = row_off + 3 * N_k;
        col_st  = col_st  + 3 * N_k;
    end

    mi.rows  = mi.rows(1:mi.n_motor);
    mi.obs   = mi.obs(:, 1:mi.n_motor);
    mi.tau   = mi.tau(1:mi.n_motor);
    mi.ommax = mi.ommax(1:mi.n_motor);

    % Snapshot elastic + active moment rows BEFORE subtracting background flow
    % or held pilus thermal forces (with matching root BC zeros), so finite
    % differences (eval_pilus_rhs - rk0) cancel constant forcing terms cleanly.
    R_local = R_global;
    row_off_bc = 6;
    for k = 1:M_pili
        N_k = N_segs(k);
        R_local(row_off_bc + 1) = 0;
        R_local(row_off_bc + 2) = 0;
        if strcmp(rc{k}, 'clamped')
            R_local(row_off_bc + 3) = 0;
        end
        row_off_bc = row_off_bc + 3 * N_k;
    end

    if has_pili_thermal
        pili_tau_in = [];
        if isfield(params, 'brown_pili_tau') && ~isempty(params.brown_pili_tau)
            pili_tau_in = params.brown_pili_tau;
        end
        R_global = R_global - pili_thermal_rhs(params.brown_pili_force, ...
            geom.R_body, hydro.B_global, params, pili_tau_in, geom.d3);
    end

    if has_bg_flow
        u_bg = coupled_background_flow(geom.X_all, geom, params);
        R_global = R_global - apply_Ainv_u_local(u_bg, hydro);
    end

    row_off_bc = 6;
    for k = 1:M_pili
        N_k = N_segs(k);
        R_global(row_off_bc + 1) = 0;
        R_global(row_off_bc + 2) = 0;
        if strcmp(rc{k}, 'clamped')
            R_global(row_off_bc + 3) = 0;
        end
        row_off_bc = row_off_bc + 3 * N_k;
    end
end


function back = apply_Ainv_u_local(u_bg, hydro)
%APPLY_AINV_U_LOCAL  Apply B_global * Ainv * u_bg using Block Schur or LU
    n_t = hydro.n_trans;
    u_t = u_bg(1:n_t);
    u_w = u_bg(n_t+1:end);

    if isfield(hydro, 'block_schur') && ~isempty(hydro.block_schur)
        bs = hydro.block_schur;
        R  = bs.R_body;
        u_t_bf = reshape(R' * reshape(u_t, 3, []), n_t, 1);
        if bs.nP == 0
            f_bf = bs.dM_BB \ u_t_bf;
            back = hydro.B_global * reshape(R * reshape(f_bf, 3, []), n_t, 1);
            return;
        end
        uB = u_t_bf(1:bs.nB);
        uP = u_t_bf(bs.nB+1:end);
        vB = bs.dM_BB \ uB;
        fP = bs.dS_PP \ (uP - bs.M_PB * vB);
        fB = vB - bs.U_BP * fP;
        Minv_u_bf = [fB; fP];
        rhs_w   = u_w - hydro.C_twist_reverse * Minv_u_bf;
        m_twist = hydro.dS_twist \ rhs_w;
        Au_t_bf = Minv_u_bf - hydro.Minv_C_twist * m_twist;
        Au_t    = reshape(R * reshape(Au_t_bf, 3, []), n_t, 1);
        Ns = numel(m_twist);
        Au_tw = zeros(3 * Ns, 1);
        Au_tw(1:3:end) = hydro.d3_all_lab(:, 1) .* m_twist;
        Au_tw(2:3:end) = hydro.d3_all_lab(:, 2) .* m_twist;
        Au_tw(3:3:end) = hydro.d3_all_lab(:, 3) .* m_twist;
        back = hydro.B_global * [Au_t; Au_tw];
        return;
    end

    if isfield(hydro, 'dM_weighted') && ~isempty(hydro.dM_weighted)
        Ainv = hydro.dM_weighted;
    else
        Ainv = hydro.M_weighted;
    end
    if isfield(hydro, 'M_twist') && ~isempty(hydro.M_twist)
        Minv_u = Ainv \ u_t;
        rhs_w  = u_w - hydro.C_twist_reverse * Minv_u;
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
