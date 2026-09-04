function [kappa1, kappa2, kappa3] = compute_pili_elastic_moments( ...
        d1, d2, d3, N_k, sigma, kb0_k, root_coupling_k) %#ok<INUSD>  (retained for interface stability, see Note)
%COMPUTE_PILI_ELASTIC_MOMENTS  Elastic curvature deviation for one pilus
%
%  Description:
%    Computes elastic restoring moments (curvature deviations from
%    intrinsic curvature) using the Kirchhoff rod model.  Follows
%    Walker (2020) dz_free_space.m lines 219-242.
%
%    The raw nodal curvature is produced by the shared operator
%    compute_pili_curvature: node-centred director frames from
%    compute_node_directors (SVD polar mean of adjacent segment frames,
%    the same nodal geometry used by build_coupled_lin_sys and
%    coupled_active_moment) projected onto arc-length derivatives that
%    are formed by differencing the two SEGMENT frames meeting at each
%    node.  Both factors are therefore collocated at the node and the
%    stencil is second-order accurate (v8.3.0; the earlier nodal-director
%    backward difference was centred half a segment off the node and
%    reduced the scheme to first order).  init_coupled_params uses the
%    same operator for kb0, so kappa - kb0 vanishes bit-exactly in the
%    reference configuration.
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
%    Row 1 (root node, s = 0) is zeroed for ALL three components and BOTH
%    root-coupling modes (v8.4.0).
%      * clamped: the three root rows are replaced by kinematic constraints
%        in build_coupled_lin_sys Step 7, so the values are placeholders.
%      * free_torsion: theta_1/phi_1 rows are again kinematic constraints.
%        The psi_1 row is the d3 moment balance of the WHOLE rod about the
%        root, in which the internal moment at s = 0 is the torque the
%        body exerts on the rod through the bearing/motor, i.e. the motor
%        torque supplied in ode_coupled_rhs (zero for a free bearing).  The
%        Kirchhoff twist moment GJ*kappa_3 is the internal moment carried
%        across a cross-section INSIDE the rod; the discrete row-1 value is
%        evaluated between segments 1 and 2 and is already the internal
%        moment used by the node-2 balance.  Adding it to row 1 (pre-8.4.0)
%        double-counted it: subtracting the node-2 row showed that the motor
%        torque was then balanced by the drag of segment 1 alone, so the
%        steady rotor speed scaled as N_k times the physical value.
%    The root_coupling_k argument is retained for interface stability.
%
%  References:
%    [1] Walker et al. (2020), Physical Review E, Eqs. (11)-(12).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-09-04 (v8.3.0: node-collocated second-order stencil via
%                       compute_pili_curvature;
%                       v8.4.0: root twist row zeroed for free_torsion)

%% Raw nodal curvature via the shared node-centred operator
%  compute_pili_curvature returns [(N_k+1) x 3] = [kappa_1, kappa_2, kappa_3]
%  using compute_node_directors for the nodal frames and a segment-
%  difference derivative collocated with each node (see its header).
kappa_raw = compute_pili_curvature(d1, d2, d3, N_k);

%% Curvature components minus intrinsic curvature
kappa1_full = kappa_raw(:, 1) - kb0_k(:, 1);
kappa2_full = kappa_raw(:, 2) - kb0_k(:, 2);
kappa3_full = kappa_raw(:, 3) - kb0_k(:, 3);

%% Zero the root row (node 1) for both root-coupling modes
%  clamped      : all three rows are replaced by kinematic constraints.
%  free_torsion : theta_1/phi_1 rows are kinematic constraints; the psi_1
%                 row is the whole-rod d3 moment balance whose boundary
%                 moment is the body->rod (motor) torque added in
%                 ode_coupled_rhs, NOT the interior Kirchhoff twist moment
%                 (see header Note).
kappa1_full(1) = 0;
kappa2_full(1) = 0;
kappa3_full(1) = 0;

%% Return segment values
kappa1 = kappa1_full(1:N_k);
kappa2 = kappa2_full(1:N_k);
kappa3 = kappa3_full(1:N_k) / (1 + sigma);

end
