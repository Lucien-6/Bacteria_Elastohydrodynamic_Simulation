function [errs_out, warns_out] = validate_custom_geometry(custom_cfg, errs, warns)
%VALIDATE_CUSTOM_GEOMETRY  Validate custom geometry configuration
%
%  Description:
%    Performs comprehensive validation of the custom geometry configuration
%    before construction.  Checks file existence and integrity, attachment
%    site consistency, and characteristic length validity.  Errors and
%    warnings are appended to the provided cell arrays and returned so
%    they can be aggregated with other validation checks by the caller
%    (validate_simulation_params).
%
%  Syntax:
%    [errs, warns] = validate_custom_geometry(custom_cfg, errs, warns)
%
%  Inputs:
%    custom_cfg - Struct with custom geometry configuration:
%      .body_points_file   - Path to body point cloud file
%      .attachment_sites   - [S x 1] struct array (optional)
%          .position       - [3 x 1] Attachment position [m]
%          .direction      - [3 x 1] Unit outward direction
%          .label          - String label
%      .pili               - [M x 1] struct array (optional)
%          .points_file    - Path to pili centerline file
%          .site_index     - Index into attachment_sites
%      .L_0_override       - Scalar [m] (optional)
%    errs  - Cell array to append error messages to
%    warns - Cell array to append warning messages to
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

%% ---- Body point cloud file ----
if ~isfield(custom_cfg, 'body_points_file') || isempty(custom_cfg.body_points_file)
    errs{end+1} = 'custom_geometry.body_points_file is required but not specified.';
    errs_out  = errs;
    warns_out = warns;
    return;   % Cannot continue without a body file
end

body_file = custom_cfg.body_points_file;
if ~isfile(body_file)
    errs{end+1} = sprintf('Body point cloud file not found: %s', body_file);
else
    % Try to read and validate the file
    try
        pts = read_geometry_csv_mat(body_file);
        if size(pts, 1) < 12
            errs{end+1} = sprintf( ...
                'Body point cloud has only %d points; minimum 12 required for surface reconstruction.', ...
                size(pts, 1));
        end
        if size(pts, 2) ~= 3
            errs{end+1} = sprintf( ...
                'Body point cloud must have 3 columns, got %d.', size(pts, 2));
        end
    catch ME
        errs{end+1} = sprintf('Failed to read body point cloud file: %s', ME.message);
    end
end

%% ---- Attachment sites ----
n_sites = 0;
if isfield(custom_cfg, 'attachment_sites') && ~isempty(custom_cfg.attachment_sites)
    sites = custom_cfg.attachment_sites;
    n_sites = numel(sites);

    % Check label uniqueness
    labels = cell(1, n_sites);
    for s = 1:n_sites
        if isfield(sites(s), 'label') && ~isempty(sites(s).label)
            labels{s} = char(sites(s).label);
        else
            labels{s} = sprintf('site_%d', s);
        end
    end
    [unique_labels, ~, ic] = unique(labels);
    if numel(unique_labels) < n_sites
        counts = accumarray(ic, 1);
        dup_mask = counts > 1;
        dup_labels = unique_labels(dup_mask);
        warns{end+1} = sprintf( ...
            'Duplicate attachment site labels: %s. This may cause confusion.', ...
            strjoin(dup_labels, ', '));
    end

    % Validate each site
    for s = 1:n_sites
        if ~isfield(sites(s), 'position') || isempty(sites(s).position)
            errs{end+1} = sprintf('Attachment site %d: .position is required.', s);
        elseif numel(sites(s).position) ~= 3
            errs{end+1} = sprintf('Attachment site %d: .position must be 3-element, got %d.', ...
                s, numel(sites(s).position));
        end

        if ~isfield(sites(s), 'direction') || isempty(sites(s).direction)
            errs{end+1} = sprintf('Attachment site %d: .direction is required.', s);
        elseif numel(sites(s).direction) ~= 3
            errs{end+1} = sprintf('Attachment site %d: .direction must be 3-element, got %d.', ...
                s, numel(sites(s).direction));
        elseif norm(sites(s).direction(:)) < 1e-30
            errs{end+1} = sprintf('Attachment site %d: .direction has zero magnitude.', s);
        end
    end
end

%% ---- Pili configurations ----
if isfield(custom_cfg, 'pili') && ~isempty(custom_cfg.pili)
    pili_cfg = custom_cfg.pili;
    n_pili = numel(pili_cfg);

    for k = 1:n_pili
        % Check points file
        if ~isfield(pili_cfg(k), 'points_file') || isempty(pili_cfg(k).points_file)
            errs{end+1} = sprintf('Custom pili %d: .points_file is required.', k);
        elseif ~isfile(pili_cfg(k).points_file)
            errs{end+1} = sprintf('Custom pili %d: centerline file not found: %s', ...
                k, pili_cfg(k).points_file);
        else
            % Try to read file
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

        % Check site index
        if ~isfield(pili_cfg(k), 'site_index') || isempty(pili_cfg(k).site_index)
            errs{end+1} = sprintf('Custom pili %d: .site_index is required.', k);
        elseif ~isnumeric(pili_cfg(k).site_index) || ...
                pili_cfg(k).site_index < 1 || ...
                pili_cfg(k).site_index > n_sites
            errs{end+1} = sprintf( ...
                'Custom pili %d: .site_index = %d is out of range [1, %d].', ...
                k, pili_cfg(k).site_index, n_sites);
        end

        % Check deflection angles if provided
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

        % Check custom bending stiffness if provided
        if isfield(pili_cfg(k), 'bending_stiffness') && ~isempty(pili_cfg(k).bending_stiffness)
            if ~isnumeric(pili_cfg(k).bending_stiffness) || ...
                    ~isscalar(pili_cfg(k).bending_stiffness) || ...
                    pili_cfg(k).bending_stiffness <= 0
                errs{end+1} = sprintf( ...
                    'Custom pili %d: bending_stiffness must be a positive scalar.', k);
            end
        end
    end

    % Check for duplicate site indices
    if n_pili > 1
        site_indices = [pili_cfg.site_index];
        [u_sites, ~, ic_p] = unique(site_indices);
        counts_p = accumarray(ic_p, 1);
        if any(counts_p > 1)
            dup_sites_p = u_sites(counts_p > 1);
            warns{end+1} = sprintf( ...
                'Multiple custom pili share the same attachment site(s): %s. This may be physically unrealistic.', ...
                mat2str(dup_sites_p));
        end
    end
else
    warns{end+1} = 'Custom geometry mode has no pili defined. The simulation will run with a bare body.';
end

%% ---- L_0 override ----
if isfield(custom_cfg, 'L_0_override') && ~isempty(custom_cfg.L_0_override)
    if ~isnumeric(custom_cfg.L_0_override) || ~isscalar(custom_cfg.L_0_override) || ...
            custom_cfg.L_0_override <= 0
        errs{end+1} = sprintf( ...
            'custom_geometry.L_0_override must be a positive scalar, got %s.', ...
            mat2str(custom_cfg.L_0_override));
    end
end

%% Return updated cell arrays (MATLAB pass-by-value requires explicit return)
errs_out  = errs;
warns_out = warns;

end
