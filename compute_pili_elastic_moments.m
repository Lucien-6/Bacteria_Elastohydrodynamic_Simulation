function [kappa1, kappa2, kappa3] = compute_pili_elastic_moments( ...
        d1, d2, d3, N_k, sigma, kb0_k)
%COMPUTE_PILI_ELASTIC_MOMENTS  Elastic curvature deviation for one pilus
%
%  Description:
%    Computes elastic restoring moments (curvature deviations from
%    intrinsic curvature) using the Kirchhoff rod model.  Follows
%    Walker (2020) dz_free_space.m lines 219-242.
%
%  Syntax:
%    [kappa1, kappa2, kappa3] = compute_pili_elastic_moments( ...
%        d1, d2, d3, N_k, sigma, kb0_k)
%
%  Inputs:
%    d1    - [N_k x 3] First director (rot_k basis)
%    d2    - [N_k x 3] Second director
%    d3    - [N_k x 3] Third director (tangent)
%    N_k   - Number of segments
%    sigma - Poisson's ratio
%    kb0_k - [N_k+1 x 3] Intrinsic curvature at nodes
%
%  Outputs:
%    kappa1 - [N_k x 1] d1-curvature deviation
%    kappa2 - [N_k x 1] d2-curvature deviation
%    kappa3 - [N_k x 1] Twist deviation / (1 + sigma)
%
%  Note:
%    The root segment (index 1) elastic moment is set to zero here as a
%    placeholder for the clamped-end constraint.  It is overwritten by the
%    clamped boundary condition equations in ode_coupled_rhs (lin_sys Step 7)
%    and must not be interpreted as a free-end condition.
%
%  References:
%    [1] Walker et al. (2020), Physical Review E, Eqs. (11)-(12).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Standardised formatting

%% Approximate directors at nodes (Walker's crude method)
d1NC = [d1; d1(end,:)];
d2NC = [d2; d2(end,:)];
d3NC = [d3; d3(end,:)];

%% Finite-difference derivatives
d1NCds = [d1NC(2,:)-d1NC(1,:); d1NC(2:end,:)-d1NC(1:end-1,:)] * N_k;
d2NCds = [d2NC(2,:)-d2NC(1,:); d2NC(2:end,:)-d2NC(1:end-1,:)] * N_k;
d3NCds = [d3NC(2,:)-d3NC(1,:); d3NC(2:end,:)-d3NC(1:end-1,:)] * N_k;

%% Curvature components minus intrinsic curvature
kappa1_full = dot(d3NC, d2NCds, 2) - kb0_k(:, 1);
kappa2_full = dot(d1NC, d3NCds, 2) - kb0_k(:, 2);
kappa3_full = dot(d2NC, d1NCds, 2) - kb0_k(:, 3);

%% Set root-segment moment to zero (placeholder, overwritten by clamped BC)
kappa1_full(1) = 0;
kappa2_full(1) = 0;
kappa3_full(1) = 0;

%% Return segment values
kappa1 = kappa1_full(1:N_k);
kappa2 = kappa2_full(1:N_k);
kappa3 = kappa3_full(1:N_k) / (1 + sigma);

end
