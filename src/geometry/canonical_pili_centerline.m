function [points_canonical, chord_length, parameter_nodes] = ...
        canonical_pili_centerline(points_raw, n_segments)
%CANONICAL_PILI_CENTERLINE  Build an equal-chord PCHIP centreline
%
%  Description:
%    Interpolates a three-dimensional input polyline with a parametric PCHIP
%    curve and selects n_segments+1 monotonically ordered parameter values
%    whose adjacent points have one common chord length.  For a trial chord
%    length, each next parameter is the first forward root of the chord
%    equation.  An outer scalar root determines the common chord for which
%    the final chord terminates exactly at the original endpoint.
%
%    This is an equal-chord construction, not equal-arc-length sampling.
%    The first and last output points are exactly the input endpoints.
%
%  Syntax:
%    points_canonical = canonical_pili_centerline(points_raw, n_segments)
%    [points_canonical, chord_length, parameter_nodes] = ...
%        canonical_pili_centerline(points_raw, n_segments)
%
%  Inputs:
%    points_raw  - [N x 3] finite, real centreline points, N >= 2
%    n_segments  - Positive integer target segment count
%
%  Outputs:
%    points_canonical - [(n_segments+1) x 3] equal-chord points
%    chord_length     - Common chord length in input coordinate units
%    parameter_nodes  - Strictly increasing PCHIP parameters in [0, 1]
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

validateattributes(points_raw, {'numeric'}, ...
    {'real', 'finite', '2d', 'ncols', 3}, mfilename, 'points_raw', 1);
validateattributes(n_segments, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'}, ...
    mfilename, 'n_segments', 2);

if size(points_raw, 1) < 2
    error('canonical_pili_centerline:TooFewPoints', ...
        'points_raw must contain at least two three-dimensional points.');
end

points_raw = double(points_raw);
origin = points_raw(1, :);
centred = points_raw - origin;
scale = max(sqrt(sum(centred.^2, 2)));
if ~isfinite(scale) || scale <= 0
    error('canonical_pili_centerline:DegenerateCurve', ...
        'The input centreline has zero spatial extent.');
end

points = centred / scale;
raw_chords = sqrt(sum(diff(points, 1, 1).^2, 2));
degenerate_tol = 1e-13;
if any(raw_chords <= degenerate_tol)
    error('canonical_pili_centerline:DegenerateInputSegment', ...
        ['The input centreline contains coincident or numerically ', ...
         'degenerate consecutive points.']);
end

path_length = sum(raw_chords);
endpoint_distance = norm(points(end, :) - points(1, :));
if endpoint_distance <= degenerate_tol
    error('canonical_pili_centerline:CoincidentEndpoints', ...
        'The input centreline endpoints must be distinct.');
end

parameter_raw = [0; cumsum(raw_chords)] / path_length;
if any(diff(parameter_raw) <= 0)
    error('canonical_pili_centerline:NonMonotoneParameter', ...
        'The PCHIP parameter must be strictly increasing.');
end

pp = cell(1, 3);
coefs = cell(1, 3);
for dim = 1:3
    pp{dim} = pchip(parameter_raw, points(:, dim));
    [breaks_dim, coefs{dim}] = unmkpp(pp{dim});
    if dim == 1
        breaks = breaks_dim;
    elseif any(breaks_dim ~= breaks)
        error('canonical_pili_centerline:InternalPchipMismatch', ...
            'PCHIP coordinate breakpoints are inconsistent.');
    end
end

if n_segments == 1
    points_canonical = points_raw([1, end], :);
    chord_length = norm(points_raw(end, :) - points_raw(1, :));
    parameter_nodes = [0; 1];
    return;
end

% Find a valid positive-residual lower bracket.
chord_low = path_length / n_segments * 1e-6;
[res_low, valid_low, nodes_low] = evaluate_trial( ...
    chord_low, n_segments, pp, coefs, breaks, points(end, :));
for attempt = 1:40
    if valid_low && res_low > 0
        break;
    end
    chord_low = chord_low / 2;
    [res_low, valid_low, nodes_low] = evaluate_trial( ...
        chord_low, n_segments, pp, coefs, breaks, points(end, :));
end
if ~valid_low || res_low <= 0
    error('canonical_pili_centerline:LowerBracketFailure', ...
        'Could not establish the lower bracket for the common chord root.');
end

best_residual = abs(res_low);
best_chord = chord_low;
best_nodes = nodes_low;

% Increase the trial chord until the residual is nonpositive or the
% requested number of forward chord roots can no longer be placed.
chord_high = max(path_length / n_segments, 2 * chord_low);
have_upper = false;
for attempt = 1:60
    [res_high, valid_high, nodes_high] = evaluate_trial( ...
        chord_high, n_segments, pp, coefs, breaks, points(end, :));
    if valid_high && abs(res_high) < best_residual
        best_residual = abs(res_high);
        best_chord = chord_high;
        best_nodes = nodes_high;
    end
    if ~valid_high || res_high <= 0
        have_upper = true;
        break;
    end
    chord_high = min(2 * chord_high, 2 * path_length);
end
if ~have_upper
    error('canonical_pili_centerline:UpperBracketFailure', ...
        'Could not establish the upper bracket for the common chord root.');
end

% Outer common-chord root.  Invalid trials belong to the upper side because
% their chord is too large to place all required intermediate points.
outer_tol = 1e-14;
for iteration = 1:150
    chord_mid = 0.5 * (chord_low + chord_high);
    [res_mid, valid_mid, nodes_mid] = evaluate_trial( ...
        chord_mid, n_segments, pp, coefs, breaks, points(end, :));

    if valid_mid
        if abs(res_mid) < best_residual
            best_residual = abs(res_mid);
            best_chord = chord_mid;
            best_nodes = nodes_mid;
        end
        if res_mid > 0
            chord_low = chord_mid;
        else
            chord_high = chord_mid;
        end
    else
        chord_high = chord_mid;
    end

    if valid_mid && abs(res_mid) <= outer_tol * max(1, chord_mid)
        break;
    end
    if chord_high - chord_low <= outer_tol * max(1, chord_mid)
        break;
    end
end

if best_residual > 5e-11 * max(1, best_chord)
    error('canonical_pili_centerline:OuterRootNoConvergence', ...
        ['The common chord root did not converge: normalized endpoint ', ...
         'residual %.3e.'], best_residual);
end

% Polish: among valid trials in the bracket, minimise segment-length spread.
endpoint_tol = 5e-11 * max(1, best_chord);
best_spread = Inf;
n_polish = 160;
for polish = 1:n_polish
    chord_trial = chord_low + (polish - 1) / (n_polish - 1) * ...
        (chord_high - chord_low);
    [res_trial, valid_trial, nodes_trial] = evaluate_trial( ...
        chord_trial, n_segments, pp, coefs, breaks, points(end, :));
    if ~valid_trial || abs(res_trial) > endpoint_tol
        continue;
    end

    trial_nodes = [nodes_trial; 1];
    trial_points = zeros(n_segments + 1, 3);
    for dim = 1:3
        trial_points(:, dim) = ppval(pp{dim}, trial_nodes);
    end
    trial_lengths = sqrt(sum(diff(trial_points, 1, 1).^2, 2));
    if any(~isfinite(trial_lengths)) || any(trial_lengths <= 0)
        continue;
    end
    rel_spread = (max(trial_lengths) - min(trial_lengths)) / ...
        mean(trial_lengths);
    if rel_spread < best_spread
        best_spread = rel_spread;
        best_chord = chord_trial;
        best_nodes = nodes_trial;
        best_residual = abs(res_trial);
    end
end

parameter_nodes = [best_nodes; 1];
if numel(parameter_nodes) ~= n_segments + 1 || ...
        any(diff(parameter_nodes) <= 0) || ...
        parameter_nodes(1) ~= 0 || parameter_nodes(end) ~= 1
    error('canonical_pili_centerline:InvalidParameterSequence', ...
        'Canonical PCHIP parameters are not strictly monotone from 0 to 1.');
end

points_normalized = zeros(n_segments + 1, 3);
for dim = 1:3
    points_normalized(:, dim) = ppval(pp{dim}, parameter_nodes);
end
points_canonical = points_normalized * scale + origin;
points_canonical(1, :) = points_raw(1, :);
points_canonical(end, :) = points_raw(end, :);

segment_lengths = sqrt(sum(diff(points_canonical, 1, 1).^2, 2));
chord_length = mean(segment_lengths);
relative_spread = (max(segment_lengths) - min(segment_lengths)) / chord_length;
if any(~isfinite(segment_lengths)) || any(segment_lengths <= 0) || ...
        relative_spread > 1e-10
    error('canonical_pili_centerline:EqualChordNoConvergence', ...
        ['Canonical segment lengths failed the equal-chord check ', ...
         '(relative spread %.3e).'], relative_spread);
end

end


function [residual, valid, parameter_nodes] = evaluate_trial( ...
        chord, n_segments, pp, coefs, breaks, endpoint)
%EVALUATE_TRIAL  Place intermediate points for one common-chord trial.

parameter_nodes = zeros(n_segments, 1);
parameter_nodes(1) = 0;
point_previous = zeros(1, 3);
for dim = 1:3
    point_previous(dim) = ppval(pp{dim}, 0);
end

valid = true;
for index = 2:n_segments
    [parameter_next, point_next, found] = next_chord_root( ...
        parameter_nodes(index - 1), point_previous, chord, ...
        pp, coefs, breaks);
    if ~found || parameter_next >= 1
        valid = false;
        residual = -Inf;
        parameter_nodes = [];
        return;
    end
    parameter_nodes(index) = parameter_next;
    point_previous = point_next;
end

residual = norm(endpoint - point_previous) - chord;

end


function [parameter_next, point_next, found] = next_chord_root( ...
        parameter_previous, point_previous, chord, pp, coefs, breaks)
%NEXT_CHORD_ROOT  Find the first forward PCHIP parameter one chord away.

n_pieces = numel(breaks) - 1;
parameter_tol = 1e-12;
piece_start = find(breaks <= parameter_previous + parameter_tol, 1, 'last');
piece_start = min(piece_start, n_pieces);

parameter_next = NaN;
point_next = nan(1, 3);
found = false;

for piece = piece_start:n_pieces
    local_left = max(0, parameter_previous - breaks(piece));
    local_right = breaks(piece + 1) - breaks(piece);

    polynomial = zeros(1, 7);
    for dim = 1:3
        component = coefs{dim}(piece, :);
        component(end) = component(end) - point_previous(dim);
        polynomial = polynomial + conv(component, component);
    end
    polynomial(end) = polynomial(end) - chord^2;

    coefficient_scale = max(abs(polynomial));
    first_nonzero = find(abs(polynomial) > ...
        1e-14 * max(1, coefficient_scale), 1, 'first');
    if isempty(first_nonzero)
        continue;
    end

    candidate_roots = roots(polynomial(first_nonzero:end));
    real_mask = abs(imag(candidate_roots)) <= ...
        1e-10 * max(1, abs(real(candidate_roots)));
    candidates = real(candidate_roots(real_mask));
    candidates = candidates(candidates > local_left + parameter_tol & ...
        candidates <= local_right + parameter_tol);
    if isempty(candidates)
        continue;
    end

    local_root = min(candidates);
    parameter_next = min(1, breaks(piece) + local_root);
    for dim = 1:3
        point_next(dim) = ppval(pp{dim}, parameter_next);
    end

    root_error = abs(norm(point_next - point_previous) - chord);
    if root_error <= 1e-10 * max(1, chord)
        found = true;
        return;
    end
end

end
