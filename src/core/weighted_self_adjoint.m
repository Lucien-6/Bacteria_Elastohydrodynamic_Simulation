function A_corr = weighted_self_adjoint(A, weights)
%WEIGHTED_SELF_ADJOINT  Project a square operator onto a weighted adjoint
%
%  A_corr = 0.5 * (A + W \ (A' * W)), where W = diag(weights).
%  Consequently, W*A_corr is self-adjoint.  For uniform weights this
%  reduces to the ordinary symmetric projection 0.5*(A + A').
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

weights = weights(:);
n = size(A, 1);

if size(A, 2) ~= n || numel(weights) ~= n
    error('weighted_self_adjoint:DimensionMismatch', ...
        'A must be square and numel(weights) must equal size(A, 1).');
end
if any(~isfinite(weights)) || any(weights <= 0)
    error('weighted_self_adjoint:InvalidWeights', ...
        'All weights must be finite and strictly positive.');
end

% Scale columns by the source measure and rows by the inverse target
% measure without allocating a dense diagonal matrix.
weighted_adjoint = (A' .* weights') ./ weights;
A_corr = 0.5 * (A + weighted_adjoint);
end
