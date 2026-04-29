function [kappa1, kappa2, kappa3] = compute_pili_elastic_moments( ...
        d1, d2, d3, N_k, sigma, kb0_k, root_coupling_k)
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
%    [kappa1, kappa2, kappa3] = compute_pili_elastic_moments( ...
%        d1, d2, d3, N_k, sigma, kb0_k, root_coupling_k)
%
%  Inputs:
%    d1    - [N_k x 3] First director (rot_k basis)
%    d2    - [N_k x 3] Second director
%    d3    - [N_k x 3] Third director (tangent)
%    N_k   - Number of segments
%    sigma - Poisson's ratio
%    kb0_k - [N_k+1 x 3] Intrinsic curvature at nodes
%    root_coupling_k - (Optional) 'clamped' (default) | 'free_torsion'
%                      Root boundary condition type for this pilus.
%
%  Outputs:
%    kappa1 - [N_k x 1] d1-curvature deviation
%    kappa2 - [N_k x 1] d2-curvature deviation
%    kappa3 - [N_k x 1] Twist deviation / (1 + sigma)
%
%  Note:
%    For clamped coupling, all three root-segment moments are set to zero
%    as placeholders overwritten by the clamped boundary condition equations
%    in ode_coupled_rhs (lin_sys Step 7).
%    For free_torsion coupling, kappa1 and kappa2 are still zeroed (theta_1
%    and phi_1 are constrained), but kappa3(1) is preserved — the torsional
%    elastic moment drives the free psi_1 dynamics.
%
%  References:
%    [1] Walker et al. (2020), Physical Review E, Eqs. (11)-(12).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06

if nargin < 7 || isempty(root_coupling_k)
    root_coupling_k = 'clamped';
end

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

%% Set root-segment moments based on root coupling type
% kappa1(1) and kappa2(1) are always zeroed (theta_1, phi_1 constrained)
kappa1_full(1) = 0;
kappa2_full(1) = 0;
if strcmp(root_coupling_k, 'free_torsion')
    % Preserve kappa3(1): torsional elastic moment drives psi_1 dynamics
else
    kappa3_full(1) = 0;  % clamped: overwritten by BC
end

%% Return segment values
kappa1 = kappa1_full(1:N_k);
kappa2 = kappa2_full(1:N_k);
kappa3 = kappa3_full(1:N_k) / (1 + sigma);

end
