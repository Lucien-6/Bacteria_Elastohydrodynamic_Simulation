function [d1N, d2N, d3N] = compute_node_directors(d1, d2, d3, N_k)
%COMPUTE_NODE_DIRECTORS  Interpolate segment directors to node directors via SVD
%
%  Description:
%    Computes node-centred director frames by averaging adjacent segment
%    director frames and projecting onto the nearest proper rotation
%    matrix via SVD polar decomposition.
%
%    Interior nodes use the arithmetic mean of adjacent segment frames.
%    End nodes use linear extrapolation from the two nearest segments.
%    Each averaged matrix is projected to SO(3) via R = U*V' from the
%    SVD, with a determinant correction to ensure det(R) > 0 (proper
%    rotation, no reflections).
%
%  Syntax:
%    [d1N, d2N, d3N] = compute_node_directors(d1, d2, d3, N_k)
%
%  Inputs:
%    d1  - [N_k x 3] First segment director (normal)
%    d2  - [N_k x 3] Second segment director (binormal)
%    d3  - [N_k x 3] Third segment director (tangent)
%    N_k - Number of segments
%
%  Outputs:
%    d1N - [(N_k+1) x 3] First nodal director
%    d2N - [(N_k+1) x 3] Second nodal director
%    d3N - [(N_k+1) x 3] Third nodal director
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-14
%  Updated: 2026-11-24 — vectorised batch pagesvd polar projection.

%% Build per-segment rotation matrices [3 x 3 x N_k]
segBases = zeros(3, 3, N_k);
segBases(:, 1, :) = d1';
segBases(:, 2, :) = d2';
segBases(:, 3, :) = d3';

%% Interpolate / extrapolate to nodes
meanRots = zeros(3, 3, N_k + 1);

% First node: linear extrapolation (boundary-safe for N_k = 1)
meanRots(:, :, 1) = 1.5 * segBases(:, :, 1) - 0.5 * segBases(:, :, min(2, N_k));

% Interior nodes: arithmetic mean of adjacent segments
if N_k >= 2
    meanRots(:, :, 2:end-1) = 0.5 * (segBases(:, :, 1:end-1) + segBases(:, :, 2:end));
end

% Last node: linear extrapolation
meanRots(:, :, end) = 1.5 * segBases(:, :, end) - 0.5 * segBases(:, :, max(1, end - 1));

%% Project to nearest proper rotation via batched SVD polar decomposition
[U_svd, ~, V_svd] = pagesvd(meanRots);
R_near = pagemtimes(U_svd, 'none', V_svd, 'transpose');

% Explicit 3x3 determinant per page; flip last column of U where det < 0
det_R = R_near(1,1,:) .* (R_near(2,2,:).*R_near(3,3,:) - R_near(2,3,:).*R_near(3,2,:)) ...
      - R_near(1,2,:) .* (R_near(2,1,:).*R_near(3,3,:) - R_near(2,3,:).*R_near(3,1,:)) ...
      + R_near(1,3,:) .* (R_near(2,1,:).*R_near(3,2,:) - R_near(2,2,:).*R_near(3,1,:));
neg_idx = find(det_R(:) < 0);
if ~isempty(neg_idx)
    U_svd(:, 3, neg_idx) = -U_svd(:, 3, neg_idx);
    R_near(:, :, neg_idx) = pagemtimes(U_svd(:, :, neg_idx), 'none', ...
                                       V_svd(:, :, neg_idx), 'transpose');
end

%% Extract nodal directors
d1N = reshape(R_near(:, 1, :), [3, N_k + 1])';
d2N = reshape(R_near(:, 2, :), [3, N_k + 1])';
d3N = reshape(R_near(:, 3, :), [3, N_k + 1])';

end
