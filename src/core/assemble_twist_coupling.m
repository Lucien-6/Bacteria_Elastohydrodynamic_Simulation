function [C_forward, C_reverse, W_trans_diag, W_twist_diag] = ...
    assemble_twist_coupling(M_rotlet, D_twist, quad_weights, segment_lengths)
%ASSEMBLE_TWIST_COUPLING  Build weighted-reciprocal twist coupling blocks
%
%  C_forward maps axial torque density to translational velocity.  The
%  reverse block is its adjoint under translational quadrature weights and
%  actual segment-length weights:
%      C_reverse = W_m \ (C_forward' * W_f).
%  W_f is represented by its diagonal to avoid constructing a large dense
%  matrix.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13
%  Updated: 2026-11-24 — exploit block-diagonal structure of D_twist (avoid dense GEMM).

quad_weights = quad_weights(:);
segment_lengths = segment_lengths(:);
W_trans_diag = repelem(quad_weights, 3);
W_twist_diag = segment_lengths;

n_seg = size(D_twist, 2);
if size(D_twist, 1) == 3 * n_seg && size(M_rotlet, 2) == 3 * n_seg && n_seg > 0
    % D_twist is block-diagonal with one 3x1 director per segment column;
    % contract each 3-column rotlet triplet directly instead of dense GEMM.
    idx_base = (0:n_seg-1) * (3 * n_seg + 3);
    d3x = D_twist(idx_base + 1);
    d3y = D_twist(idx_base + 2);
    d3z = D_twist(idx_base + 3);
    C_forward = 0.5 * (M_rotlet(:, 1:3:end) .* d3x ...
                     + M_rotlet(:, 2:3:end) .* d3y ...
                     + M_rotlet(:, 3:3:end) .* d3z);
else
    C_forward = 0.5 * (M_rotlet * D_twist);
end

if size(C_forward, 1) ~= numel(W_trans_diag) || ...
        size(C_forward, 2) ~= numel(W_twist_diag)
    error('assemble_twist_coupling:DimensionMismatch', ...
        ['Coupling dimensions must match three translational components per ', ...
         'quadrature point and one twist component per segment.']);
end
if any(~isfinite(W_trans_diag)) || any(W_trans_diag <= 0) || ...
        any(~isfinite(W_twist_diag)) || any(W_twist_diag <= 0)
    error('assemble_twist_coupling:InvalidWeights', ...
        'Translational and twist weights must be finite and strictly positive.');
end

% Right multiplication by W_f scales columns of C_forward'.  Left solution
% with diagonal W_m scales rows, without allocating either dense matrix.
C_reverse = C_forward' .* W_trans_diag';
C_reverse = C_reverse ./ W_twist_diag;
end
