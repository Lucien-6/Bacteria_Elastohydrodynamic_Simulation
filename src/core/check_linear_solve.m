function backward_residual = check_linear_solve( ...
    A, x, b, context, warning_threshold, error_threshold)
%CHECK_LINEAR_SOLVE  Validate finiteness and normwise backward residual
%
%  The scale-independent diagnostic is
%      ||A*x-b||_inf / (||A||_inf*||x||_inf + ||b||_inf).
%  By default, residuals above 1e-10 warn and residuals above 1e-6 error.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

if nargin < 4 || isempty(context)
    context = 'linear solve';
end
if nargin < 5 || isempty(warning_threshold)
    warning_threshold = 1e-10;
end
if nargin < 6 || isempty(error_threshold)
    error_threshold = 1e-6;
end

if any(~isfinite(A(:))) || any(~isfinite(x(:))) || any(~isfinite(b(:)))
    error('numerical_diagnostics:NonFiniteValue', ...
        'Non-finite value detected in %s.', context);
end
if size(A, 2) ~= size(x, 1) || size(A, 1) ~= size(b, 1) || ...
        size(x, 2) ~= size(b, 2)
    error('numerical_diagnostics:DimensionMismatch', ...
        'Incompatible dimensions in %s.', context);
end
if ~(warning_threshold >= 0 && error_threshold > warning_threshold)
    error('numerical_diagnostics:InvalidThresholds', ...
        'Residual thresholds must satisfy 0 <= warning < error.');
end

residual = A * x - b;
numerator = norm(residual, inf);
denominator = norm(A, inf) * norm(x, inf) + norm(b, inf);
if denominator == 0
    if numerator == 0
        backward_residual = 0;
    else
        backward_residual = Inf;
    end
else
    backward_residual = numerator / denominator;
end
if isfinite(backward_residual)
    numerical_diagnostics_state('residual', backward_residual);
end

if ~isfinite(backward_residual)
    error('numerical_diagnostics:NonFiniteResidual', ...
        'Non-finite backward residual detected in %s.', context);
elseif backward_residual > error_threshold
    check_matrix_condition(A, sprintf('%s condition', context));
    error('numerical_diagnostics:LargeBackwardResidual', ...
        'Backward residual %.3e exceeds %.3e in %s.', ...
        backward_residual, error_threshold, context);
elseif backward_residual > warning_threshold
    check_matrix_condition(A, sprintf('%s condition', context));
    warning('numerical_diagnostics:ElevatedBackwardResidual', ...
        'Backward residual %.3e exceeds %.3e in %s.', ...
        backward_residual, warning_threshold, context);
end
end
