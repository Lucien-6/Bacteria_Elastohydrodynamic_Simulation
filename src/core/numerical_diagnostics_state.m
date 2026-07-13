function diagnostics = numerical_diagnostics_state(action, varargin)
%NUMERICAL_DIAGNOSTICS_STATE  Run-scoped numerical diagnostic accumulator
%
%  The ODE solver repeatedly invokes pure RHS functions that cannot return
%  auxiliary statistics.  This run-scoped accumulator records only scalar
%  health metrics and is reset explicitly by each simulation entry point.
%  Parallel workers have independent MATLAB processes and therefore
%  independent persistent states.
%
%  Syntax:
%    numerical_diagnostics_state('reset')
%    numerical_diagnostics_state('residual', value)
%    numerical_diagnostics_state('rcond', value)
%    numerical_diagnostics_state('hydro', reciprocity, dissipation)
%    diagnostics = numerical_diagnostics_state('snapshot')
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

persistent state
if isempty(state)
    state = default_state();
end

switch lower(char(action))
    case 'reset'
        state = default_state();

    case 'residual'
        value = require_scalar(varargin, 1, 'residual');
        state.max_backward_residual = max(state.max_backward_residual, value);
        state.linear_solve_count = state.linear_solve_count + 1;

    case 'rcond'
        value = require_scalar(varargin, 1, 'rcond');
        state.min_rcond = min(state.min_rcond, value);
        state.rcond_check_count = state.rcond_check_count + 1;

    case 'hydro'
        reciprocity = require_scalar(varargin, 1, 'reciprocity');
        dissipation = require_scalar(varargin, 2, 'dissipation');
        state.max_weighted_reciprocity_residual = max( ...
            state.max_weighted_reciprocity_residual, reciprocity);
        state.min_dissipation_eigenvalue = min( ...
            state.min_dissipation_eigenvalue, dissipation);
        state.hydrodynamic_check_count = ...
            state.hydrodynamic_check_count + 1;

    case 'snapshot'
        % No mutation.

    otherwise
        error('numerical_diagnostics:UnknownAction', ...
            'Unknown diagnostics action "%s".', char(action));
end

diagnostics = state;
end


function state = default_state()
state = struct( ...
    'backward_residual_warning_threshold', 1e-10, ...
    'backward_residual_error_threshold', 1e-6, ...
    'rcond_warning_threshold', 1e-10, ...
    'rcond_error_threshold', 1e-14, ...
    'reciprocity_error_threshold', 1e-10, ...
    'dissipation_error_relative_threshold', 1e-10, ...
    'max_backward_residual', 0, ...
    'min_rcond', Inf, ...
    'max_weighted_reciprocity_residual', 0, ...
    'min_dissipation_eigenvalue', Inf, ...
    'linear_solve_count', 0, ...
    'rcond_check_count', 0, ...
    'hydrodynamic_check_count', 0);
end


function value = require_scalar(args, index, name)
if numel(args) < index
    error('numerical_diagnostics:MissingMetric', ...
        'Missing %s metric.', name);
end
value = double(args{index});
if ~isscalar(value) || ~isreal(value) || ~isfinite(value)
    error('numerical_diagnostics:InvalidMetric', ...
        '%s metric must be a finite real scalar.', name);
end
end
