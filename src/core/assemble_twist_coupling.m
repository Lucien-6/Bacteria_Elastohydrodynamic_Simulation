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

quad_weights = quad_weights(:);
segment_lengths = segment_lengths(:);
W_trans_diag = repelem(quad_weights, 3);
W_twist_diag = segment_lengths;

C_forward = 0.5 * (M_rotlet * D_twist);

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
C_reverse = bsxfun(@times, C_forward', W_trans_diag');
C_reverse = bsxfun(@rdivide, C_reverse, W_twist_diag);
end
