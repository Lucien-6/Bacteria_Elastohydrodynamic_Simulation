function [errs_out, warns_out] = validate_custom_geometry(custom_cfg, errs, warns)
%VALIDATE_CUSTOM_GEOMETRY  Validate custom geometry configuration
%
%  Description:
%    Performs comprehensive pre-flight validation of the custom geometry
%    configuration.  Checks:
%      - Body point cloud file existence and integrity
%      - Attachment sites (labels, positions, directions; near-surface and
%        outward-pointing checks)
%      - Pili centerline file integrity
%      - Site indices and deflection angles
%      - Rejection of the retired per-pilus bending_stiffness interface
%      - L_0 override, L_0_suggest_mode
%      - initial_orientation spec (quaternion/axis_angle/euler_zyx)
%
%    Errors and warnings are appended to the provided cell arrays and
%    returned so they can be aggregated with other validation checks by
%    the caller (validate_simulation_params).
%    Body points and attachment positions are expected in the same
%    body-local, COM-centred frame.  Pili files are read in their source
%    Cartesian input frame; validation does not transform coordinates.
%    All positional coordinates and lengths follow the project input
%    convention of SI metres; this validator performs no unit conversion.
%
%  Syntax:
%    [errs, warns] = validate_custom_geometry(custom_cfg, errs, warns)
%
%  Inputs:
%    custom_cfg - Struct with custom geometry configuration
%    errs       - Cell array to append error messages to
%    warns      - Cell array to append warning messages to
%
%  Outputs:
%    errs_out  - Updated cell array with error messages
%    warns_out - Updated cell array with warning messages
%
%  Note:
%    MATLAB passes cell arrays by value, so the updated arrays must be
%    returned and captured by the caller.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-28
%  Updated: 2026-05-04 — Uniform-segment check, site-on-surface diagnostics,
%                         cell/struct compat, initial_orientation.
%  Updated: 2026-07-13

%% ---- Body point cloud file ----
if ~isfield(custom_cfg, 'body_points_file') || isempty(custom_cfg.body_points_file)
    errs{end+1} = 'custom_geometry.body_points_file is required but not specified.';
    errs_out  = errs;
    warns_out = warns;
    return;
end

body_file = custom_cfg.body_points_file;
body_pts  = [];
if ~isfile(body_file)
    errs{end+1} = sprintf('Body point cloud file not found: %s', body_file);
else
    try
        body_pts = read_geometry_csv_mat(body_file);
        if size(body_pts, 1) < 12
            errs{end+1} = sprintf( ...
                'Body point cloud has only %d points; minimum 12 required for surface reconstruction.', ...
                size(body_pts, 1));
        end
        if size(body_pts, 2) ~= 3
            errs{end+1} = sprintf( ...
                'Body point cloud must have 3 columns, got %d.', size(body_pts, 2));
        end
    catch ME
        errs{end+1} = sprintf('Failed to read body point cloud file: %s', ME.message);
    end
end

%% ---- Attachment sites (handle cell/struct polymorphism) ----
n_sites = 0;
sites = [];
if isfield(custom_cfg, 'attachment_sites') && ~isempty(custom_cfg.attachment_sites)
    try
        sites = normalise_json_struct_array(custom_cfg.attachment_sites);
        n_sites = numel(sites);
    catch ME
        errs{end+1} = sprintf('attachment_sites could not be normalised: %s', ME.message);
        sites = [];
        n_sites = 0;
    end

    if n_sites > 0
        % Label uniqueness check
        labels = cell(1, n_sites);
        for s = 1:n_sites
            if isfield(sites(s), 'label') && ~isempty(sites(s).label)
                labels{s} = char(sites(s).label);
            else
                labels{s} = sprintf('site_%d', s);
            end
        end
        [u_labels, ~, ic_l] = unique(labels);
        if numel(u_labels) < n_sites
            counts_l  = accumarray(ic_l, 1);
            dup_labels = u_labels(counts_l > 1);
            warns{end+1} = sprintf( ...
                'Duplicate attachment site labels: %s. This may cause confusion.', ...
                strjoin(dup_labels, ', '));
        end

        % Precompute body surface metrics for site-on-surface checks
        have_body_pts = (~isempty(body_pts)) && (size(body_pts, 2) == 3) && ...
            (size(body_pts, 1) >= 12);

        for s = 1:n_sites
            pos_ok = true;  dir_ok = true;

            if ~isfield(sites(s), 'position') || isempty(sites(s).position)
                errs{end+1} = sprintf('Attachment site %d: .position is required.', s);
                pos_ok = false;
            elseif numel(sites(s).position) ~= 3
                errs{end+1} = sprintf( ...
                    'Attachment site %d: .position must be 3-element, got %d.', ...
                    s, numel(sites(s).position));
                pos_ok = false;
            end

            if ~isfield(sites(s), 'direction') || isempty(sites(s).direction)
                errs{end+1} = sprintf('Attachment site %d: .direction is required.', s);
                dir_ok = false;
            elseif numel(sites(s).direction) ~= 3
                errs{end+1} = sprintf( ...
                    'Attachment site %d: .direction must be 3-element, got %d.', ...
                    s, numel(sites(s).direction));
                dir_ok = false;
            elseif norm(sites(s).direction(:)) < 1e-30
                errs{end+1} = sprintf('Attachment site %d: .direction has zero magnitude.', s);
                dir_ok = false;
            end

            % Site-on-surface and outward-direction diagnostics (m7, m8)
            if pos_ok && dir_ok && have_body_pts
                [warn_msg_list] = check_site_geometry( ...
                    body_pts, sites(s).position(:), sites(s).direction(:), s, labels{s});
                for mm = 1:numel(warn_msg_list)
                    warns{end+1} = warn_msg_list{mm}; %#ok<AGROW>
                end
            end
        end
    end
end

%% ---- Pili configurations ----
pili_cfg = [];
n_pili = 0;
if isfield(custom_cfg, 'pili') && ~isempty(custom_cfg.pili)
    try
        pili_cfg = normalise_json_struct_array(custom_cfg.pili);
        n_pili   = numel(pili_cfg);
    catch ME
        errs{end+1} = sprintf('custom_geometry.pili could not be normalised: %s', ME.message);
        pili_cfg = [];
        n_pili   = 0;
    end

    for k = 1:n_pili
        % Points file
        if ~isfield(pili_cfg(k), 'points_file') || isempty(pili_cfg(k).points_file)
            errs{end+1} = sprintf('Custom pili %d: .points_file is required.', k);
        elseif ~isfile(pili_cfg(k).points_file)
            errs{end+1} = sprintf('Custom pili %d: centerline file not found: %s', ...
                k, pili_cfg(k).points_file);
        else
            try
                ppts = read_geometry_csv_mat(pili_cfg(k).points_file);
                if size(ppts, 1) < 2
                    errs{end+1} = sprintf( ...
                        'Custom pili %d: centerline has only %d points; minimum 2 required.', ...
                        k, size(ppts, 1));
                end
            catch ME
                errs{end+1} = sprintf('Custom pili %d: failed to read centerline file: %s', ...
                    k, ME.message);
            end
        end

        % Site index
        if ~isfield(pili_cfg(k), 'site_index') || isempty(pili_cfg(k).site_index)
            errs{end+1} = sprintf('Custom pili %d: .site_index is required.', k);
        elseif ~isnumeric(pili_cfg(k).site_index) || ...
                pili_cfg(k).site_index < 1 || ...
                pili_cfg(k).site_index > n_sites
            errs{end+1} = sprintf( ...
                'Custom pili %d: .site_index = %g is out of range [1, %d].', ...
                k, pili_cfg(k).site_index, n_sites);
        end

        % Deflection angles (optional, scalar check)
        if isfield(pili_cfg(k), 'deflection_z') && ~isempty(pili_cfg(k).deflection_z)
            if ~isnumeric(pili_cfg(k).deflection_z) || ~isscalar(pili_cfg(k).deflection_z)
                errs{end+1} = sprintf('Custom pili %d: deflection_z must be a scalar.', k);
            end
        end
        if isfield(pili_cfg(k), 'deflection_x') && ~isempty(pili_cfg(k).deflection_x)
            if ~isnumeric(pili_cfg(k).deflection_x) || ~isscalar(pili_cfg(k).deflection_x)
                errs{end+1} = sprintf('Custom pili %d: deflection_x must be a scalar.', k);
            end
        end

        % Per-pilus stiffness was a dead interface.  Material stiffness is
        % defined once from the global E and r_pili inputs.
        if isfield(pili_cfg(k), 'bending_stiffness') && ~isempty(pili_cfg(k).bending_stiffness)
            errs{end+1} = sprintf( ...
                ['Custom pili %d: bending_stiffness is not supported. ', ...
                 'Set the global material E and r_pili values instead.'], k);
        end
    end

    % Duplicate site_index
    if n_pili > 1
        valid_idx_mask = false(1, n_pili);
        site_indices   = nan(1, n_pili);
        for k = 1:n_pili
            if isfield(pili_cfg(k), 'site_index') && ...
                    isnumeric(pili_cfg(k).site_index) && ...
                    isscalar(pili_cfg(k).site_index)
                valid_idx_mask(k) = true;
                site_indices(k)   = pili_cfg(k).site_index;
            end
        end
        if any(valid_idx_mask)
            [u_sites, ~, ic_p] = unique(site_indices(valid_idx_mask));
            counts_p = accumarray(ic_p, 1);
            if any(counts_p > 1)
                dup_sites_p = u_sites(counts_p > 1);
                warns{end+1} = sprintf( ...
                    'Multiple custom pili share the same attachment site(s): %s. This may be physically unrealistic.', ...
                    mat2str(dup_sites_p));
            end
        end
    end
else
    warns{end+1} = 'Custom geometry mode has no pili defined. The simulation will run with a bare body.';
end

%% ---- L_0 override / L_0_suggest_mode ----
if isfield(custom_cfg, 'L_0_override') && ~isempty(custom_cfg.L_0_override)
    if ~isnumeric(custom_cfg.L_0_override) || ~isscalar(custom_cfg.L_0_override) || ...
            custom_cfg.L_0_override <= 0
        errs{end+1} = sprintf( ...
            'custom_geometry.L_0_override must be a positive scalar, got %s.', ...
            mat2str(custom_cfg.L_0_override));
    end
end

if isfield(custom_cfg, 'L_0_suggest_mode') && ~isempty(custom_cfg.L_0_suggest_mode)
    mode_str = char(custom_cfg.L_0_suggest_mode);
    if ~any(strcmpi(mode_str, {'volume', 'extent'}))
        errs{end+1} = sprintf( ...
            'custom_geometry.L_0_suggest_mode must be ''volume'' or ''extent'', got ''%s''.', ...
            mode_str);
    end
end

%% ---- initial_orientation (C2) ----
if isfield(custom_cfg, 'initial_orientation') && ~isempty(custom_cfg.initial_orientation)
    iospec = custom_cfg.initial_orientation;
    if ~isstruct(iospec)
        errs{end+1} = sprintf( ...
            'custom_geometry.initial_orientation must be a struct, got %s.', ...
            class(iospec));
    else
        try
            parse_initial_orientation(iospec);   % reject bad specs
        catch ME
            errs{end+1} = sprintf( ...
                'custom_geometry.initial_orientation is invalid: %s', ME.message);
        end
    end
end

%% Return updated cell arrays
errs_out  = errs;
warns_out = warns;

end


%% ==================== Local helper: site-on-surface diagnostic ====================

function warns_new = check_site_geometry(body_pts, pos, dir_vec, s_idx, label)
%CHECK_SITE_GEOMETRY  Warn if attachment site is far from surface or points inward.
%
%  - Finds the nearest body point and reports the distance relative to
%    the mean nearest-neighbour spacing h_mean.
%  - Estimates the local outward normal via the smallest-eigenvector of
%    the covariance of the nearest-K neighbours (PCA), oriented away
%    from the body centroid.  Warns when dir_vec makes an angle > 90 deg
%    with this normal (i.e. dir_vec points inward).

    warns_new = {};

    % Distance to nearest body point
    d2 = sum((body_pts - pos').^2, 2);
    [d2_min, idx_nn] = min(d2);
    d_min = sqrt(d2_min);

    % Mean body spacing (sample a subset for speed)
    Nmax = min(size(body_pts, 1), 600);
    sub_idx = round(linspace(1, size(body_pts, 1), Nmax));
    sub_pts = body_pts(sub_idx, :);
    D2 = sum((reshape(sub_pts, [Nmax, 1, 3]) - ...
              reshape(sub_pts, [1, Nmax, 3])).^2, 3);
    D2(1:Nmax+1:end) = Inf;
    h_mean_est = mean(sqrt(min(D2, [], 2)));

    if d_min > 5 * h_mean_est
        warns_new{end+1} = sprintf( ...
            ['Attachment site %d (''%s''): position is %.2e m from the ', ...
             'nearest body point (%.1fx h_mean=%.2e m). Likely off-surface.'], ...
            s_idx, label, d_min, d_min / h_mean_est, h_mean_est);
    end

    % Local outward-normal estimate via PCA on K nearest neighbours
    K = min(20, size(body_pts, 1));
    [~, nn_order] = sort(d2);
    nn_idx = nn_order(1:K);
    nbrs = body_pts(nn_idx, :);
    nbrs_c = nbrs - mean(nbrs, 1);
    [V_eig, D_eig] = eig(nbrs_c' * nbrs_c);
    [~, i_min]  = min(diag(D_eig));
    local_normal = V_eig(:, i_min);

    % Orient outward (away from body centroid)
    body_centroid = mean(body_pts, 1)';
    nearest_pt    = body_pts(idx_nn, :)';
    if dot(local_normal, nearest_pt - body_centroid) < 0
        local_normal = -local_normal;
    end

    dir_u = dir_vec / max(norm(dir_vec), 1e-30);
    cos_ang = dot(dir_u, local_normal);
    if cos_ang < 0
        warns_new{end+1} = sprintf( ...
            ['Attachment site %d (''%s''): direction has negative projection ', ...
             '(cos(theta) = %.3f) on the estimated local outward normal. ', ...
             'The specified direction appears to point INTO the body.'], ...
            s_idx, label, cos_ang);
    end

end
