function [Body_str, body_meta] = custom_body_building(custom_cfg)
%CUSTOM_BODY_BUILDING  Construct custom bacterial body from user point cloud
%
%  Description:
%    Processes a user-provided custom body point cloud file (CSV or MAT),
%    validates the data (including duplicate-point removal), computes
%    geometric properties (mean spacing, volume via triangulation),
%    automatically triangulates the surface for visualisation and volume
%    calculation, and processes user-defined attachment sites for pili
%    connection.
%
%    In the Regularised Stokeslet Method (RSM), force points and
%    evaluation points are collocated, so only one set of surface
%    coordinates is needed.  The triangulated surface mesh is used for
%    visualisation (trisurf) and enclosed-volume computation (divergence
%    theorem); it is NOT used in the hydrodynamic kernel evaluation.
%
%  Syntax:
%    [Body_str, body_meta] = custom_body_building(custom_cfg)
%
%  Input:
%    custom_cfg - Struct with fields:
%      .body_points_file   - Path to CSV/MAT file with N x 3 point cloud [m]
%      .attachment_sites   - Struct array / cell array (optional)
%          .position       - [3 x 1] Attachment position in body-local frame [m]
%          .direction      - [3 x 1] Unit outward direction at attachment
%          .label          - String label (e.g. 'left_pole') (optional)
%      .L_0_override       - Scalar characteristic length [m] (optional)
%      .L_0_suggest_mode   - 'volume' (default) | 'extent' (optional)
%                             volume : equivalent sphere radius from V_body
%                             extent : half of the maximum bounding-box extent
%      .duplicate_tol      - Scalar relative tolerance for point deduplication
%                             (default: 1e-6 of max extent)
%
%  Outputs:
%    Body_str   - Structure with surface point coordinates:
%                   .rxB, .ryB, .rzB  Surface point coordinates [N_body x 1] [m]
%    body_meta  - Structure with geometric metadata:
%                   .N_body           Number of surface points (post-dedup)
%                   .volume           Enclosed volume [m^3] (from triangulation)
%                   .faces            Triangulation face indices [F x 3]
%                   .vertices         Triangulation vertices [V x 3] [m]
%                   .attachment_sites Processed attachment sites struct array
%                   .h_mean           Mean nearest-neighbour spacing [m]
%                   .h_min            Minimum nearest-neighbour spacing [m]
%                   .L_0_suggested    Suggested characteristic length [m]
%                   .L_0_suggest_mode Source of L_0_suggested ('volume'/'extent')
%                   .alpha_shape      alphaShape handle (if available; [] otherwise)
%                   .triangulation_method 'alphaShape' | 'convhulln'
%                   .n_duplicates_removed Integer count of deduplicated points
%
%  References:
%    [1] Cortez et al. (2005), Phys. Fluids 17, 031504.
%
%  See also: read_geometry_csv_mat, capsule_body_building, normalise_json_struct_array
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-28
%  Updated: 2026-05-04 — Duplicate-point removal, volume-based L_0
%                        default, Euler characteristic check, expanded
%                        alpha scan, cell/struct-array compatible sites.

%% Read point cloud data
pts_raw = read_geometry_csv_mat(custom_cfg.body_points_file);
N_raw   = size(pts_raw, 1);

fprintf('  [Custom Body] Read %d points from: %s\n', ...
    N_raw, custom_cfg.body_points_file);

%% Deduplicate (near-)identical points (M4)
if isfield(custom_cfg, 'duplicate_tol') && ~isempty(custom_cfg.duplicate_tol)
    dup_tol = custom_cfg.duplicate_tol;
else
    dup_tol = 1e-6 * max(max(pts_raw) - min(pts_raw));
    dup_tol = max(dup_tol, 1e-15);
end
pts = uniquetol(pts_raw, dup_tol, 'ByRows', true, 'DataScale', 1);
n_dup_removed = N_raw - size(pts, 1);
if n_dup_removed > 0
    warning('custom_body_building:DuplicatePoints', ...
        ['%d duplicate / near-duplicate point(s) removed (tol = %.3e m). ', ...
         'Verify the source point cloud.'], n_dup_removed, dup_tol);
end

N_body = size(pts, 1);

%% Validate point cloud quality
if N_body < 12
    error('custom_body_building:InsufficientPoints', ...
          ['Custom body requires at least 12 unique surface points, got %d ', ...
           '(after deduplication of %d raw points).'], N_body, N_raw);
end

%% Compute nearest-neighbour spacing (characteristic mesh spacing)
%  Vectorised O(N^2); acceptable for typical N < 5000.
D2 = sum((reshape(pts, [N_body, 1, 3]) - ...
          reshape(pts, [1, N_body, 3])).^2, 3);
D2(1:N_body+1:end) = Inf;   % Exclude self-distances
min_dists = sqrt(min(D2, [], 2));

h_min  = min(min_dists);
h_mean = mean(min_dists);

% Safety check: after deduplication this should never trigger, but belt & braces
if ~(h_min > 0 && isfinite(h_min))
    error('custom_body_building:DegenerateSpacing', ...
        ['Minimum nearest-neighbour spacing is %.3e m, which is non-positive ', ...
         'or non-finite.  Point cloud may contain residual duplicate points ', ...
         'below the deduplication tolerance (%.3e m).'], h_min, dup_tol);
end

fprintf('  [Custom Body] Spacing: h_min = %.4e m, h_mean = %.4e m\n', ...
    h_min, h_mean);

%% Automatically triangulate surface (for visualisation + volume)
[faces, vertices, volume, tri_method, shp] = auto_triangulate_surface(pts);

fprintf('  [Custom Body] Triangulation (%s): %d vertices, %d faces\n', ...
    tri_method, size(vertices, 1), size(faces, 1));
fprintf('  [Custom Body] Enclosed volume: %.4e m^3\n', volume);

% Euler characteristic / closed-surface diagnostic
check_closed_surface(faces, size(vertices, 1));

%% Process attachment sites (handle struct / cell / struct-array)
sites = [];
if isfield(custom_cfg, 'attachment_sites') && ~isempty(custom_cfg.attachment_sites)
    sites_raw = normalise_json_struct_array(custom_cfg.attachment_sites);
    n_sites   = numel(sites_raw);

    % Copy to working struct array and validate / normalise each site
    sites = sites_raw;
    for s = 1:n_sites
        if ~isfield(sites(s), 'position') || isempty(sites(s).position)
            error('custom_body_building:MissingPosition', ...
                'Attachment site %d: .position is required.', s);
        end
        if ~isfield(sites(s), 'direction') || isempty(sites(s).direction)
            error('custom_body_building:MissingDirection', ...
                'Attachment site %d: .direction is required.', s);
        end

        pos     = sites(s).position(:);
        dir_vec = sites(s).direction(:);

        dir_norm = norm(dir_vec);
        if dir_norm < 1e-30
            error('custom_body_building:ZeroDirection', ...
                  'Attachment site %d has zero direction vector.', s);
        end
        sites(s).direction = dir_vec / dir_norm;
        sites(s).position  = pos;

        if ~isfield(sites(s), 'label') || isempty(sites(s).label)
            sites(s).label = sprintf('site_%d', s);
        else
            sites(s).label = char(sites(s).label);
        end
    end

    fprintf('  [Custom Body] Attachment sites: %d defined\n', n_sites);
    for s = 1:n_sites
        fprintf(['    Site %d (''%s''): pos = [%.3e, %.3e, %.3e] m, ', ...
                 'dir = [%.3f, %.3f, %.3f]\n'], ...
            s, sites(s).label, ...
            sites(s).position(1), sites(s).position(2), sites(s).position(3), ...
            sites(s).direction(1), sites(s).direction(2), sites(s).direction(3));
    end
else
    fprintf('  [Custom Body] No attachment sites defined.\n');
end

%% Compute suggested characteristic length (m9: prefer volume-based)
suggest_mode = 'volume';
if isfield(custom_cfg, 'L_0_suggest_mode') && ~isempty(custom_cfg.L_0_suggest_mode)
    suggest_mode = lower(char(custom_cfg.L_0_suggest_mode));
end

if isfield(custom_cfg, 'L_0_override') && ~isempty(custom_cfg.L_0_override)
    L_0_suggested = custom_cfg.L_0_override;
    suggest_mode  = 'override';
else
    extent_range  = max(pts) - min(pts);
    L_0_extent    = max(extent_range) / 2;
    if isfinite(volume) && volume > 0
        L_0_volume = (3 * volume / (4 * pi))^(1/3);
    else
        L_0_volume = NaN;
    end

    switch suggest_mode
        case 'volume'
            if isnan(L_0_volume)
                warning('custom_body_building:VolumeLengthFallback', ...
                    'L_0_suggest_mode = ''volume'' requested but volume is invalid; falling back to ''extent''.');
                L_0_suggested = L_0_extent;
                suggest_mode  = 'extent';
            else
                L_0_suggested = L_0_volume;
            end
        case 'extent'
            L_0_suggested = L_0_extent;
        otherwise
            error('custom_body_building:BadSuggestMode', ...
                'Unknown L_0_suggest_mode ''%s''. Expected ''volume'' or ''extent''.', suggest_mode);
    end
end

%% Pack output structure
Body_str = struct('rxB', pts(:,1), 'ryB', pts(:,2), 'rzB', pts(:,3));

body_meta.N_body                = N_body;
body_meta.volume                = volume;
body_meta.faces                 = faces;
body_meta.vertices              = vertices;
body_meta.attachment_sites      = sites;
body_meta.h_mean                = h_mean;
body_meta.h_min                 = h_min;
body_meta.L_0_suggested         = L_0_suggested;
body_meta.L_0_suggest_mode      = suggest_mode;
body_meta.alpha_shape           = shp;
body_meta.triangulation_method  = tri_method;
body_meta.n_duplicates_removed  = n_dup_removed;

end


%% ==================== Local Helper Functions ====================

function [faces, vertices, volume, method, shp] = auto_triangulate_surface(pts)
%AUTO_TRIANGULATE_SURFACE  Reconstruct triangulated surface from scattered points
%
%  Uses alphaShape with a wide logarithmic alpha search.  Falls back to
%  convhulln if alphaShape fails or produces a non-manifold result.
%  Computes enclosed volume via the divergence theorem.

    shp    = [];
    method = 'convhulln';

    % --- Try alphaShape with expanded logarithmic alpha scan (m12) ---
    [faces, vertices, alpha_ok, shp_try] = try_alpha_shape(pts);
    if alpha_ok
        method = 'alphaShape';
        shp    = shp_try;
    else
        fprintf(['  [Custom Body] alphaShape failed or non-manifold; ', ...
                 'falling back to convhulln.\n']);
        try
            faces    = convhulln(pts);
            vertices = pts;
        catch ME
            warning('custom_body_building:TriangulationFailed', ...
                ['Surface triangulation failed: %s. Visualisation will ', ...
                 'use scatter plot; volume = NaN.'], ME.message);
            faces    = zeros(0, 3);
            vertices = pts;
            volume   = NaN;
            return;
        end
    end

    % --- Compute enclosed volume ---
    if ~isempty(faces) && size(faces, 2) == 3
        volume = compute_triangulation_volume(vertices, faces);
    else
        volume = NaN;
    end

end


function [faces, vertices, ok, shp] = try_alpha_shape(pts)
%TRY_ALPHA_SHAPE  Attempt alphaShape surface reconstruction with expanded scan

    ok       = false;
    faces    = zeros(0, 3);
    vertices = pts;
    shp      = [];

    try
        shp = alphaShape(pts);

        % Try the default alpha first
        [faces, vertices] = boundaryFacets(shp);
        if ~isempty(faces) && is_manifold_tri(faces, size(vertices, 1))
            ok = true;
            return;
        end

        % Widened logarithmic scan: two orders of magnitude around alpha_crit
        alpha_crit = criticalAlpha(shp, 'all-points');
        if isfinite(alpha_crit) && alpha_crit > 0
            alpha_vals = logspace(log10(0.1 * alpha_crit), ...
                                   log10(10   * alpha_crit), 15);
            for ai = 1:numel(alpha_vals)
                shp.Alpha = alpha_vals(ai);
                [faces, vertices] = boundaryFacets(shp);
                if ~isempty(faces) && size(faces, 2) == 3 && ...
                        is_manifold_tri(faces, size(vertices, 1))
                    ok = true;
                    return;
                end
            end
        end

    catch
        ok = false;
    end

end


function ok = is_manifold_tri(faces, n_verts)
%IS_MANIFOLD_TRI  Basic manifold check: each edge shared by at most 2 faces

    if isempty(faces) || size(faces, 2) ~= 3
        ok = false;
        return;
    end

    edges = [faces(:, [1,2]); faces(:, [1,3]); faces(:, [2,3])];
    edges = sort(edges, 2);

    [unique_edges, ~, ic] = unique(edges, 'rows');
    counts = accumarray(ic, 1);

    ok = max(counts) <= 2 && ...
         all(unique_edges(:) >= 1) && ...
         all(unique_edges(:) <= n_verts);

end


function vol = compute_triangulation_volume(vertices, faces)
%COMPUTE_TRIANGULATION_VOLUME  Enclosed volume via divergence theorem
%
%  V = (1/6) * | sum_over_faces det([v1; v2; v3]) |
%
%  The absolute value accommodates triangulations whose global orientation
%  is inward rather than outward.  Inconsistently-oriented meshes (mix
%  of inward/outward) still give a non-physical result; such meshes
%  trigger the manifold / Euler-characteristic diagnostics elsewhere.

    v1 = vertices(faces(:, 1), :);
    v2 = vertices(faces(:, 2), :);
    v3 = vertices(faces(:, 3), :);

    det_vals = v1(:,1).*(v2(:,2).*v3(:,3) - v3(:,2).*v2(:,3)) ...
             - v2(:,1).*(v1(:,2).*v3(:,3) - v3(:,2).*v1(:,3)) ...
             + v3(:,1).*(v1(:,2).*v2(:,3) - v2(:,2).*v1(:,3));

    vol = abs(sum(det_vals)) / 6;

end


function check_closed_surface(faces, n_verts)
%CHECK_CLOSED_SURFACE  Diagnostic: Euler characteristic and face/vertex ratio.
%
%  For a genus-0 closed surface (sphere topology): V - E + F == 2, and
%  F == 2*(V - 2) from Euler's formula combined with 2E == 3F for a
%  triangulation.  Deviations indicate holes, self-intersections, or
%  non-manifold geometry that will corrupt volume calculation.

    if isempty(faces) || size(faces, 2) ~= 3
        return;
    end

    F = size(faces, 1);
    edges = [faces(:, [1,2]); faces(:, [1,3]); faces(:, [2,3])];
    edges = sort(edges, 2);
    E = size(unique(edges, 'rows'), 1);
    V = n_verts;

    chi = V - E + F;    % Euler characteristic
    if chi ~= 2
        warning('custom_body_building:NonClosedSurface', ...
            ['Triangulation has Euler characteristic chi = %d (expected 2 ', ...
             'for a genus-0 closed surface). V=%d, E=%d, F=%d. The body ', ...
             'mesh may have holes or be non-manifold; enclosed volume ', ...
             'may be inaccurate.'], chi, V, E, F);
    end

    % Face-count sanity: for closed sphere topology F ~= 2*(V - 2)
    F_expected = 2 * (V - 2);
    if F_expected > 0 && abs(F - F_expected) > 0.25 * F_expected
        warning('custom_body_building:UnexpectedFaceCount', ...
            ['Face count F = %d deviates substantially from the closed-', ...
             'surface expectation 2*(V-2) = %d. The triangulation may ', ...
             'have extra / missing faces.'], F, F_expected);
    end

end
