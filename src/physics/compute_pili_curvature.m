function kappa_raw = compute_pili_curvature(d1, d2, d3, N_k)
%COMPUTE_PILI_CURVATURE  Raw nodal curvature / twist of one discrete pilus
%
%  Description:
%    Evaluates the three Kirchhoff curvature components
%        kappa_1 = d3 . d2'    kappa_2 = d1 . d3'    kappa_3 = d2 . d1'
%    at the N_k+1 nodes of a pilus discretised into N_k equal segments
%    (ds = 1/N_k in pilus-length units).  This is the SINGLE curvature
%    operator shared by the run-time elastic moments
%    (compute_pili_elastic_moments) and by the intrinsic-curvature
%    initialisation (init_coupled_params), so that kappa_live - kb0 is
%    bit-exact zero in the reference configuration.
%
%    Stencil (node-centred, second order):
%      * Director values d_i at node i (i = 1..N_k+1) come from
%        compute_node_directors (SVD polar mean of the two adjacent segment
%        frames; linear extrapolation at the two end nodes).
%      * Arc-length derivatives at interior node i (i = 2..N_k) are formed
%        by differencing the two SEGMENT frames that meet at that node:
%            d'(s_i) ~ (d_seg(i) - d_seg(i-1)) / ds
%        Segment i-1 and segment i are centred at s_i -/+ ds/2, so this
%        difference is a genuine central stencil about node i, i.e.
%        O(ds^2) accurate and collocated with the nodal directors it is
%        projected onto.  (Differencing the NODAL directors themselves
%        would centre the derivative at s_i - ds/2, half a segment away
%        from the node, and reduce the projection to first order.)
%      * The end nodes (1 and N_k+1) replicate the nearest interior
%        stencil (one-sided).  Node 1 is overwritten by the root boundary
%        condition downstream and node N_k+1 is never used (moment-free
%        tip), so their formal accuracy is irrelevant.
%
%  Syntax:
%    kappa_raw = compute_pili_curvature(d1, d2, d3, N_k)
%
%  Inputs:
%    d1, d2, d3 - [N_k x 3] Segment material directors (rot_k basis or
%                 lab frame; the result is frame-invariant).
%    N_k        - Number of segments (>= 1).
%
%  Outputs:
%    kappa_raw  - [(N_k+1) x 3] Nodal curvature/twist in units of
%                 1/L_pili: columns are [kappa_1, kappa_2, kappa_3].
%
%  Notes:
%    No intrinsic curvature is subtracted and no boundary rows are
%    zeroed here; callers apply those rules.
%
%  References:
%    [1] Walker et al. (2020), Physical Review E, Eqs. (11)-(12).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-09-04

%% Node-centred directors (shared SVD-projected interpolation)
[d1N, d2N, d3N] = compute_node_directors(d1, d2, d3, N_k);

%% Segment-difference derivatives collocated with interior nodes
if N_k < 2
    % A single segment carries no bending information: return zeros so
    % downstream code sees a straight, untwisted rod.
    kappa_raw = zeros(N_k + 1, 3);
    return;
end

dd1 = (d1(2:end, :) - d1(1:end-1, :)) * N_k;   % [(N_k-1) x 3] at nodes 2..N_k
dd2 = (d2(2:end, :) - d2(1:end-1, :)) * N_k;
dd3 = (d3(2:end, :) - d3(1:end-1, :)) * N_k;

% Pad to N_k+1 rows: node 1 borrows the node-2 stencil, node N_k+1 borrows
% the node-N_k stencil (one-sided ends, see header).
d1ds = [dd1(1, :); dd1; dd1(end, :)];
d2ds = [dd2(1, :); dd2; dd2(end, :)];
d3ds = [dd3(1, :); dd3; dd3(end, :)];

%% Curvature components (Kirchhoff definitions)
kappa_raw = [dot(d3N, d2ds, 2), ...
             dot(d1N, d3ds, 2), ...
             dot(d2N, d1ds, 2)];

end
