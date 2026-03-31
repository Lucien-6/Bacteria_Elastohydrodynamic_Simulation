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
%  Updated: 2026-02-14 -- RSS rotlet segment coupling + vectorised assembly
%           2026-02-14 -- BUG FIX: added rot_k in Q base-velocity mapping
%           2026-03-31 -- Refactored: Steps 1-7 extracted to build_coupled_lin_sys

%% Extract frequently used parameters
N_segs = params.N_segments;
sigma  = params.sigma;
Eg     = params.Eg;
g_hat  = params.g_hat(:);

%% Steps 1-7: Build coupled linear system (geometry + hydrodynamics + clamped BC)
[lin_sys, geom] = build_coupled_lin_sys(Z, params);

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
    for j = 1:N_k
        R_global(row_offset + 3*(j-1) + 1) = kappa1(j);
        R_global(row_offset + 3*(j-1) + 2) = kappa2(j);
        R_global(row_offset + 3*(j-1) + 3) = kappa3(j);
    end
    row_offset = row_offset + 3 * N_k;
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
