function diagnostics = check_coupled_system_health(Z, params, context)
%CHECK_COUPLED_SYSTEM_HEALTH  Condition, reciprocity, and dissipation checks
%
%  This intentionally expensive check is called at the initial state and
%  after an Euler-basis change, not at every ODE RHS evaluation.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

if nargin < 3 || isempty(context)
    context = 'coupled system';
end
if any(~isfinite(Z(:)))
    error('numerical_diagnostics:NonFiniteState', ...
        'Non-finite state detected in %s.', context);
end

[lin_sys, ~, hydro] = build_coupled_lin_sys(Z, params);
condition_reciprocal = check_matrix_condition(lin_sys, context);
mobility = check_weighted_mobility_health( ...
    hydro.M_weighted, hydro.C_twist, hydro.C_twist_reverse, ...
    hydro.M_twist, hydro.W_trans_diag, hydro.W_twist_diag, context);

diagnostics = struct( ...
    'rcond', condition_reciprocal, ...
    'weighted_reciprocity_residual', ...
        mobility.weighted_reciprocity_residual, ...
    'min_dissipation_eigenvalue', mobility.min_dissipation_eigenvalue, ...
    'dissipation_scale', mobility.dissipation_scale);
end
