function metrics = check_weighted_mobility_health( ...
    Mtt, C_forward, C_reverse, M_twist, w_force, w_twist, context)
%CHECK_WEIGHTED_MOBILITY_HEALTH  Reciprocity and dissipation diagnostics
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13
%  Updated: 2026-07-20

if nargin < 7 || isempty(context)
    context = 'weighted mobility';
end
w_force = w_force(:);
w_twist = w_twist(:);
if any(~isfinite(Mtt(:))) || any(~isfinite(w_force)) || ...
        any(w_force <= 0)
    error('numerical_diagnostics:NonFiniteMatrix', ...
        'Invalid translational mobility or weights in %s.', context);
end
if ~isequal(size(Mtt), [numel(w_force), numel(w_force)])
    error('numerical_diagnostics:DimensionMismatch', ...
        'Translational mobility and force weights mismatch in %s.', context);
end

weighted_Mtt = Mtt .* w_force;
reciprocity_residual = 0;
grand_weighted = weighted_Mtt;
if ~isempty(M_twist)
    if any(~isfinite(M_twist(:))) || any(~isfinite(C_forward(:))) || ...
            any(~isfinite(C_reverse(:))) || any(~isfinite(w_twist)) || ...
            any(w_twist <= 0)
        error('numerical_diagnostics:NonFiniteMatrix', ...
            'Invalid twist mobility, coupling, or weights in %s.', context);
    end

    cross_left = C_forward .* w_force;
    cross_right = C_reverse' .* w_twist';
    cross_scale = max([norm(cross_left, 'fro'), ...
                       norm(cross_right, 'fro'), realmin('double')]);
    cross_residual = norm(cross_left - cross_right, 'fro') / cross_scale;

    twist_left = M_twist .* w_twist;
    twist_right = M_twist' .* w_twist';
    twist_scale = max([norm(twist_left, 'fro'), ...
                       norm(twist_right, 'fro'), realmin('double')]);
    twist_residual = norm(twist_left - twist_right, 'fro') / twist_scale;
    reciprocity_residual = max(cross_residual, twist_residual);

    grand_weighted = [weighted_Mtt, cross_left; ...
                      cross_right', twist_left];
end

power_symmetric = 0.5 * (grand_weighted + grand_weighted');
dissipation_scale = max(norm(power_symmetric, 'fro'), realmin('double'));
min_dissipation = smallest_symmetric_eigenvalue(power_symmetric, context);

if reciprocity_residual > 1e-10
    error('numerical_diagnostics:ReciprocityFailure', ...
        ['Weighted hydrodynamic reciprocity residual %.3e exceeds ', ...
         '1e-10 in %s.'], reciprocity_residual, context);
elseif reciprocity_residual > 1e-12
    warning('numerical_diagnostics:ReciprocityRoundoff', ...
        'Weighted reciprocity residual %.3e in %s.', ...
        reciprocity_residual, context);
end
if min_dissipation < -1e-10 * dissipation_scale
    error('numerical_diagnostics:NegativeDissipation', ...
        ['Minimum symmetric weighted-mobility eigenvalue %.3e is ', ...
         'materially negative relative to scale %.3e in %s.'], ...
        min_dissipation, dissipation_scale, context);
elseif min_dissipation < 0
    warning('numerical_diagnostics:DissipationRoundoff', ...
        'Round-off-scale negative dissipation eigenvalue %.3e in %s.', ...
        min_dissipation, context);
end
numerical_diagnostics_state('hydro', reciprocity_residual, min_dissipation);

metrics = struct( ...
    'weighted_reciprocity_residual', reciprocity_residual, ...
    'min_dissipation_eigenvalue', min_dissipation, ...
    'dissipation_scale', dissipation_scale);
end


function lambda_min = smallest_symmetric_eigenvalue(A, context)
%SMALLEST_SYMMETRIC_EIGENVALUE  Robust lambda_min for dissipation diagnostics
%
%  RSM/RSS weighted mobility operators are dense.  Iterative eigs with
%  'smallestreal' is unreliable on dense, clustered Stokes spectra and must
%  never abort a physically valid solve merely because the eigensolver fails.

n = size(A, 1);
use_iterative = issparse(A) && n > 600;
if use_iterative
    options.issym = true;
    options.isreal = true;
    options.tol = 1e-9;
    options.maxit = 1000;
    [~, value, flag] = eigs(A, 1, 'smallestreal', options);
    if flag == 0
        lambda_min = real(value(1, 1));
        return;
    end
    warning('numerical_diagnostics:DissipationEigensolveFallback', ...
        ['Iterative smallest-eigenvalue solve did not converge in %s; ', ...
         'falling back to dense eig.'], context);
end

lambda_min = real(min(eig(full(A), 'vector')));
end
