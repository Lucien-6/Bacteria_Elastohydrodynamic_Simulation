function diagnostics = check_coupled_system_health(Z, params, context)
%CHECK_COUPLED_SYSTEM_HEALTH  Condition, reciprocity, and dissipation checks
%
%  This intentionally expensive check is called at the initial state and
%  after an Euler-basis change, not at every ODE RHS evaluation.
%
%  rcond alone never aborts.  A gravity/unit probe solve is checked through
%  check_linear_solve so singularity is gated jointly with backward residual.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13
%  Updated: 2026-07-21

if nargin < 3 || isempty(context)
    context = 'coupled system';
end
if any(~isfinite(Z(:)))
    error('numerical_diagnostics:NonFiniteState', ...
        'Non-finite state detected in %s.', context);
end

[lin_sys, ~, hydro] = build_coupled_lin_sys(Z, params);
condition_reciprocal = check_matrix_condition(lin_sys, context);

% Joint solve-quality gate: usable Stokes systems may have tiny rcond.
R_probe = zeros(size(lin_sys, 1), 1);
if isfield(params, 'Eg') && isnumeric(params.Eg) && isscalar(params.Eg) && ...
        isfinite(params.Eg) && params.Eg ~= 0
    R_probe(3) = -params.Eg;
else
    R_probe(1) = 1;
end
x_probe = lin_sys \ R_probe;
backward_residual = check_linear_solve( ...
    lin_sys, x_probe, R_probe, sprintf('%s probe solve', context));

mobility = check_weighted_mobility_health( ...
    hydro.M_weighted, hydro.C_twist, hydro.C_twist_reverse, ...
    hydro.M_twist, hydro.W_trans_diag, hydro.W_twist_diag, context);

diagnostics = struct( ...
    'rcond', condition_reciprocal, ...
    'backward_residual', backward_residual, ...
    'weighted_reciprocity_residual', ...
        mobility.weighted_reciprocity_residual, ...
    'min_dissipation_eigenvalue', mobility.min_dissipation_eigenvalue, ...
    'dissipation_scale', mobility.dissipation_scale);
end
