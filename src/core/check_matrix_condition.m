function condition_reciprocal = check_matrix_condition(A, context)
%CHECK_MATRIX_CONDITION  Validate finiteness and reciprocal condition number
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

if nargin < 2 || isempty(context)
    context = 'matrix';
end
if any(~isfinite(A(:)))
    error('numerical_diagnostics:NonFiniteMatrix', ...
        'Non-finite value detected in %s.', context);
end
if size(A, 1) ~= size(A, 2)
    error('numerical_diagnostics:NonSquareMatrix', ...
        'Condition check requires a square matrix in %s.', context);
end

condition_reciprocal = rcond(A);
if ~isfinite(condition_reciprocal)
    error('numerical_diagnostics:NonFiniteCondition', ...
        'Non-finite rcond detected in %s.', context);
elseif condition_reciprocal < 1e-14
    error('numerical_diagnostics:SingularMatrix', ...
        'rcond %.3e is below 1e-14 in %s.', ...
        condition_reciprocal, context);
elseif condition_reciprocal < 1e-10
    warning('numerical_diagnostics:IllConditionedMatrix', ...
        'rcond %.3e is below 1e-10 in %s.', ...
        condition_reciprocal, context);
end
numerical_diagnostics_state('rcond', condition_reciprocal);
end
