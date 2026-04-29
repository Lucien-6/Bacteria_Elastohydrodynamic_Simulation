function [Body_str, body_meta] = custom_body_building(custom_cfg)
%CUSTOM_BODY_BUILDING  Construct custom bacterial body from user point cloud
%
%  Description:
%    Processes a user-provided custom body point cloud file (CSV or MAT),
%    validates the data, computes geometric properties (mean spacing, volume),
%    automatically triangulates the surface for visualisation, and processes
%    user-defined attachment sites for pili connection.
%
%    In the Regularised Stokeslet Method (RSM), force points and evaluation
%    points are collocated, so only one set of surface coordinates is needed.
%    The triangulated surface mesh is generated solely for visualisation
%    (trisurf in run_simulation) and is NOT used in the hydrodynamic
%    computation.
%
%  Syntax:
%    [Body_str, body_meta] = custom_body_building(custom_cfg)
%
%  Input:
%    custom_cfg - Struct with fields:
%      .body_points_file   - Path to CSV/MAT file with N x 3 point cloud [m]
%      .attachment_sites   - [S x 1] struct array (optional, default: empty)
%          .position       - [3 x 1] Attachment position in body-local frame [m]
%          .direction      - [3 x 1] Unit outward direction at attachment [m]
%          .label          - String label (e.g. 'left_pole') (optional)
%      .L_0_override       - Scalar characteristic length [m] (optional)
%
%  Outputs:
%    Body_str   - Structure with surface point coordinates:
%                   .rxB, .ryB, .rzB  Surface point coordinates [N_body x 1] [m]
%    body_meta  - Structure with geometric metadata:
%                   .N_body           Number of surface points
%                   .volume           Enclosed volume [m^3] (from triangulation)
%                   .faces            Triangulation face indices [F x 3]
%                   .vertices         Triangulation vertices [V x 3] [m]
%                   .attachment_sites Processed attachment sites struct array
%                   .h_mean           Mean nearest-neighbour spacing [m]
%                   .h_min            Minimum nearest-neighbour spacing [m]
%                   .L_0_suggested    Suggested characteristic length [m]
%
%  References:
%    [1] Cortez et al. (2005), Phys. Fluids 17, 031504.
%
%  See also: read_geometry_csv_mat, capsule_body_building
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-28

%% Read point cloud data
pts = read_geometry_csv_mat(custom_cfg.body_points_file);
N_body = size(pts, 1);

fprintf('  [Custom Body] Read %d points from: %s\n', N_body, custom_cfg.body_points_file);

%% Validate point cloud quality
if N_body < 12
    error('custom_body_building:InsufficientPoints', ...
          'Custom body requires at least 12 surface points, got %d.', N_body);
end

%% Compute nearest-neighbour spacing (characteristic mesh spacing)
D2 = sum((reshape(pts, [N_body, 1, 3]) - ...
          reshape(pts, [1, N_body, 3])).^2, 3);
D2(1:N_body+1:end) = Inf;   % Exclude self-distances
min_dists = sqrt(min(D2, [], 2));

h_min  = min(min_dists);
h_mean = mean(min_dists);

fprintf('  [Custom Body] Spacing: h_min = %.4e m, h_mean = %.4e m\n', h_min, h_mean);

%% Automatically triangulate surface for visualisation
[faces, vertices, volume] = auto_triangulate_surface(pts);

fprintf('  [Custom Body] Triangulation: %d vertices, %d faces\n', ...
    size(vertices, 1), size(faces, 1));
fprintf('  [Custom Body] Enclosed volume: %.4e m^3\n', volume);

%% Process attachment sites
sites = [];
if isfield(custom_cfg, 'attachment_sites') && ~isempty(custom_cfg.attachment_sites)
    n_sites = numel(custom_cfg.attachment_sites);
    sites = custom_cfg.attachment_sites;

    % Validate each site
    for s = 1:n_sites
        pos = sites(s).position(:);
        dir_vec = sites(s).direction(:);

        % Normalise direction
        dir_norm = norm(dir_vec);
        if dir_norm < 1e-30
            error('custom_body_building:ZeroDirection', ...
                  'Attachment site %d has zero direction vector.', s);
        end
        sites(s).direction = dir_vec / dir_norm;

        % Ensure position is a column vector
        sites(s).position = pos;

        % Default label if not provided
        if ~isfield(sites(s), 'label') || isempty(sites(s).label)
            sites(s).label = sprintf('site_%d', s);
        end
    end

    fprintf('  [Custom Body] Attachment sites: %d defined\n', n_sites);
    for s = 1:n_sites
        fprintf('    Site %d (''%s''): pos = [%.3e, %.3e, %.3e] m, dir = [%.3f, %.3f, %.3f]\n', ...
            s, sites(s).label, ...
            sites(s).position(1), sites(s).position(2), sites(s).position(3), ...
            sites(s).direction(1), sites(s).direction(2), sites(s).direction(3));
    end
else
    fprintf('  [Custom Body] No attachment sites defined.\n');
end

%% Compute suggested characteristic length
if isfield(custom_cfg, 'L_0_override') && ~isempty(custom_cfg.L_0_override)
    L_0_suggested = custom_cfg.L_0_override;
else
    % Use half the maximum extent as the default characteristic length
    extent_range = max(pts) - min(pts);
    L_0_suggested = max(extent_range) / 2;
end

%% Pack output structure
Body_str = struct('rxB', pts(:,1), 'ryB', pts(:,2), 'rzB', pts(:,3));

body_meta.N_body           = N_body;
body_meta.volume            = volume;
body_meta.faces             = faces;
body_meta.vertices          = vertices;
body_meta.attachment_sites  = sites;
body_meta.h_mean            = h_mean;
body_meta.h_min             = h_min;
body_meta.L_0_suggested     = L_0_suggested;

end


%% ==================== Local Helper Functions ====================

function [faces, vertices, volume] = auto_triangulate_surface(pts)
%AUTO_TRIANGULATE_SURFACE  Reconstruct triangulated surface from scattered points
%
%  Uses alphaShape with volume-guided alpha selection.  Falls back to
%  convhulln if alphaShape fails or produces a non-manifold result.
%  Computes enclosed volume from the triangulation.

    N = size(pts, 1);

    % --- Try alphaShape with automatic alpha selection ---
    [faces, vertices, volume, alpha_ok] = try_alpha_shape(pts);

    % --- Fallback: convex hull via convhulln ---
    if ~alpha_ok
        fprintf('  [Custom Body] alphaShape failed or non-manifold; falling back to convhulln.\n');
        try
            faces = convhulln(pts);
            vertices = pts;
            alpha_ok = true;
        catch ME
            warning('custom_body_building:TriangulationFailed', ...
                'Surface triangulation failed: %s. Visualisation will use scatter plot.', ME.message);
            faces = zeros(0, 3);
            vertices = pts;
            volume = NaN;
            return;
        end
    end

    % --- Compute enclosed volume from triangulation ---
    if alpha_ok && ~isempty(faces) && size(faces, 2) == 3
        volume = compute_triangulation_volume(vertices, faces);
    else
        volume = NaN;
    end

end


function [faces, vertices, volume, ok] = try_alpha_shape(pts)
%TRY_ALPHA_SHAPE  Attempt alphaShape surface reconstruction with auto alpha

    ok = false;
    faces = zeros(0, 3);
    vertices = pts;
    volume = NaN;

    try
        shp = alphaShape(pts);

        % Try default alpha first
        [faces, vertices] = boundaryFacets(shp);

        if ~isempty(faces) && is_manifold_tri(faces, size(vertices, 1))
            ok = true;
            return;
        end

        % If default fails, scan alpha values
        alpha_crit = criticalAlpha(shp, 'all-points');
        alpha_vals = linspace(0.5 * alpha_crit, 2.0 * alpha_crit, 10);

        for ai = 1:numel(alpha_vals)
            shp.Alpha = alpha_vals(ai);
            [faces, vertices] = boundaryFacets(shp);

            if ~isempty(faces) && size(faces, 2) == 3 && ...
                    is_manifold_tri(faces, size(vertices, 1))
                ok = true;
                return;
            end
        end

    catch
        % alphaShape not available or failed
        ok = false;
    end

end


function ok = is_manifold_tri(faces, n_verts)
%IS_MANIFOLD_TRI  Basic manifold check: each edge shared by at most 2 faces

    if isempty(faces) || size(faces, 2) ~= 3
        ok = false;
        return;
    end

    % Build edge list (sorted vertex pairs)
    edges = [faces(:, [1,2]); faces(:, [1,3]); faces(:, [2,3])];
    edges = sort(edges, 2);

    % Count edge occurrences
    [unique_edges, ~, ic] = unique(edges, 'rows');
    counts = accumarray(ic, 1);

    % Manifold: every edge shared by at most 2 faces
    ok = max(counts) <= 2 && all(unique_edges(:) >= 1) && all(unique_edges(:) <= n_verts);

end


function vol = compute_triangulation_volume(vertices, faces)
%COMPUTE_TRIANGULATION_VOLUME  Compute enclosed volume via divergence theorem
%
%  V = (1/6) * sum_over_faces( det([v1; v2; v3]) )
%  where v1, v2, v3 are the 3 vertices of each triangular face.

    v1 = vertices(faces(:, 1), :);
    v2 = vertices(faces(:, 2), :);
    v3 = vertices(faces(:, 3), :);

    % Signed volume contribution of each tetrahedron (origin to face)
    det_vals = v1(:,1).*(v2(:,2).*v3(:,3) - v3(:,2).*v2(:,3)) ...
             - v2(:,1).*(v1(:,2).*v3(:,3) - v3(:,2).*v1(:,3)) ...
             + v3(:,1).*(v1(:,2).*v2(:,3) - v2(:,2).*v1(:,3));

    vol = abs(sum(det_vals)) / 6;

end
