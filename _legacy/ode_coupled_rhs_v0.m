function dZ = ode_coupled_rhs_v0(t, Z, params)
%ODE_COUPLED_RHS_V0  Frozen reference copy (pre-Direction-1) — historical archive only
%
%  Description:
%    Bit-for-bit copy of ode_coupled_rhs.m at HEAD before the Direction 1
%    (LU decomposition-reuse) optimisation.  This version calls
%    build_coupled_lin_sys_v0 internally, so the entire legacy solve path
%    (including hydro.M_weighted \ u_t inside apply_Ainv_u_v0) is self-
%    contained.  It originally served the v5.5.5 invariance/benchmark
%    harnesses; those harnesses have since been removed after validation
%    (see CHANGELOG.md [5.5.5] "Removed").  The file is retained only as
%    a historical reference to the pre-optimisation behaviour.
%
%    Never invoke this function from the live simulator.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-21

%% Extract frequently used parameters
N_segs = params.N_segments;
sigma  = params.sigma;
Eg     = params.Eg;
g_hat  = params.g_hat(:);

%% Steps 1-7: Build coupled linear system (geometry + hydrodynamics + clamped BC)
has_bg_flow = ~strcmp(params.bg_flow.type, 'off');
if has_bg_flow
    [lin_sys, geom, hydro] = build_coupled_lin_sys_v0(Z, params);
else
    [lin_sys, geom] = build_coupled_lin_sys_v0(Z, params);
end

%% Step 8: Compute elastic restoring moments (RHS)
n_state  = 6 + 3 * sum(N_segs);
R_global = zeros(n_state, 1);
R_global(1:3) = -Eg * g_hat;

row_offset = 6;
for k = 1:params.M_pili
    N_k = N_segs(k);
    [kappa1, kappa2, kappa3] = compute_pili_elastic_moments( ...
        geom.d1_basis{k}, geom.d2_basis{k}, geom.d3_basis{k}, ...
        N_k, sigma, params.kb0{k});
    idx_k = row_offset + (1:3*N_k)';
    R_global(idx_k) = reshape([kappa1'; kappa2'; kappa3'], [], 1);

    %% Step 8b: Subtract active internal moments (per-pili, time-varying)
    if ~strcmp(params.am_cfg{k}.type, 'off')
        d1_k  = geom.d1_basis{k};
        d2_k  = geom.d2_basis{k};
        d3_k  = geom.d3_basis{k};
        [d1N_k, d2N_k, d3N_k] = compute_node_directors(d1_k, d2_k, d3_k, N_k);
        M_act = coupled_active_moment(t, k, N_k, d1_k, d2_k, d3_k, ...
                                      d1N_k, d2N_k, d3N_k, params.am_cfg{k});
        R_global(idx_k) = R_global(idx_k) - M_act(:);
    end

    row_offset = row_offset + 3 * N_k;
end

%% Step 8c: Subtract background flow contribution
if has_bg_flow
    u_bg = coupled_background_flow(geom.X_all, geom, params);
    back = apply_Ainv_u_v0(u_bg, hydro);
    R_global = R_global - back;
end

%% Apply clamped BC RHS zeros
col_offset = 6;
for k = 1:params.M_pili
    N_k = N_segs(k);
    r1 = col_offset + 1;
    R_global(r1)   = 0;
    R_global(r1+1) = 0;
    R_global(r1+2) = 0;
    col_offset = col_offset + 3 * N_k;
end

%% Step 9: Solve
z_dot = lin_sys \ R_global;

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

function back = apply_Ainv_u_v0(u_bg, hydro)
%APPLY_AINV_U_V0  Legacy (pre-Direction-1) B_global * Ainv * u_bg
    n_t = hydro.n_trans;
    u_t = u_bg(1:n_t);
    u_w = u_bg(n_t+1:end);

    if ~isempty(hydro.A_twist_total)
        Aw = hydro.A_twist_total * u_w;
        if ~isempty(hydro.M_rotlet)
            u_t = u_t - hydro.M_rotlet * Aw;
        end
        Au_t = hydro.M_weighted \ u_t;
        Au   = [Au_t; Aw];
    else
        Au = hydro.M_weighted \ u_t;
    end

    back = hydro.B_global * Au;
end
