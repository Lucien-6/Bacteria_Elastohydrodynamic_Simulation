function [Pili_Matrix, pili_ext] = parse_pili_config(pili_cfg)
%PARSE_PILI_CONFIG  Convert named pili struct array to Pili_Matrix [8 x M]
%
%  Description:
%    Translates a user-friendly struct array (with named fields) into the
%    [8 x M] numeric Pili_Matrix expected by downstream functions such as
%    pili_building and init_coupled_params.
%
%    Optionally returns a second output with per-pili extended fields
%    (intrinsic_curvature, active_moment) that are not encoded in the
%    numeric Pili_Matrix.
%
%  Syntax:
%    Pili_Matrix = parse_pili_config(pili_cfg)
%    [Pili_Matrix, pili_ext] = parse_pili_config(pili_cfg)
%
%  Input:
%    pili_cfg - [1 x M] struct array.  Each element may contain:
%      .type           - Morphology type: string ('linear','parabolic',
%                        'circular','sinusoidal','spiral') OR integer 1-5
%      .length_um      - Pilus length [um] (must be > 0)
%      .deflection_z   - Deflection angle from Z-axis [rad] (default: 0)
%      .deflection_x   - Deflection angle from X-axis [rad] (default: 0)
%      .shape_params   - [1x3] shape parameters (default: [NaN NaN NaN])
%      .site           - Attachment site index, 1-14
%      .intrinsic_curvature - (Optional) 'frozen' | struct with .type etc.
%      .active_moment       - (Optional) struct('type','off') | etc.
%
%  Output:
%    Pili_Matrix - [8 x M] numeric matrix:
%      Row 1: type, Row 2: length [um], Row 3: deflection Z [rad],
%      Row 4: deflection X [rad], Rows 5-7: shape params, Row 8: site
%    pili_ext    - [1 x M] struct array with .intrinsic_curvature and
%                  .active_moment per pili (default: 'frozen' / 'off')
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17

if isempty(pili_cfg)
    Pili_Matrix = zeros(8, 0);
    pili_ext    = struct('intrinsic_curvature', {}, 'active_moment', {});
    return;
end

type_map = struct( ...
    'linear',     1, ...
    'parabolic',  2, ...
    'circular',   3, ...
    'sinusoidal', 4, ...
    'spiral',     5);

M = numel(pili_cfg);
Pili_Matrix = zeros(8, M);

for k = 1:M
    pc = pili_cfg(k);

    % --- Type ---
    if isfield(pc, 'type') && ~isempty(pc.type)
        if ischar(pc.type) || isstring(pc.type)
            type_str = lower(char(pc.type));
            if isfield(type_map, type_str)
                Pili_Matrix(1, k) = type_map.(type_str);
            else
                error('parse_pili_config:InvalidType', ...
                    'Pili %d: unknown type ''%s''. Valid: linear, parabolic, circular, sinusoidal, spiral.', ...
                    k, type_str);
            end
        elseif isnumeric(pc.type) && isscalar(pc.type)
            if pc.type < 1 || pc.type > 5 || floor(pc.type) ~= pc.type
                error('parse_pili_config:InvalidType', ...
                    'Pili %d: numeric type must be integer 1-5, got %g.', k, pc.type);
            end
            Pili_Matrix(1, k) = pc.type;
        else
            error('parse_pili_config:InvalidType', ...
                'Pili %d: type must be a string or integer 1-5.', k);
        end
    else
        error('parse_pili_config:MissingType', 'Pili %d: field ''type'' is required.', k);
    end

    % --- Length ---
    if isfield(pc, 'length_um') && ~isempty(pc.length_um)
        if pc.length_um <= 0
            error('parse_pili_config:InvalidLength', ...
                'Pili %d: length_um must be > 0, got %g.', k, pc.length_um);
        end
        Pili_Matrix(2, k) = pc.length_um;
    else
        error('parse_pili_config:MissingLength', 'Pili %d: field ''length_um'' is required.', k);
    end

    % --- Deflection angles (optional, default 0) ---
    if isfield(pc, 'deflection_z') && ~isempty(pc.deflection_z)
        Pili_Matrix(3, k) = pc.deflection_z;
    else
        Pili_Matrix(3, k) = 0;
    end

    if isfield(pc, 'deflection_x') && ~isempty(pc.deflection_x)
        Pili_Matrix(4, k) = pc.deflection_x;
    else
        Pili_Matrix(4, k) = 0;
    end

    % --- Shape parameters (optional, default [NaN NaN NaN]) ---
    if isfield(pc, 'shape_params') && ~isempty(pc.shape_params)
        sp = pc.shape_params(:)';
        if numel(sp) > 3
            error('parse_pili_config:InvalidShapeParams', ...
                'Pili %d: shape_params must have at most 3 elements.', k);
        end
        sp_full = [NaN NaN NaN];
        sp_full(1:numel(sp)) = sp;
        Pili_Matrix(5:7, k) = sp_full';
    else
        Pili_Matrix(5:7, k) = [NaN; NaN; NaN];
    end

    % --- Attachment site ---
    if isfield(pc, 'site') && ~isempty(pc.site)
        if pc.site < 1 || pc.site > 14 || floor(pc.site) ~= pc.site
            error('parse_pili_config:InvalidSite', ...
                'Pili %d: site must be integer 1-14, got %g.', k, pc.site);
        end
        Pili_Matrix(8, k) = pc.site;
    else
        error('parse_pili_config:MissingSite', 'Pili %d: field ''site'' is required.', k);
    end

    % --- Extended fields (pass-through for init_coupled_params) ---
    if isfield(pc, 'intrinsic_curvature') && ~isempty(pc.intrinsic_curvature)
        pili_ext(k).intrinsic_curvature = pc.intrinsic_curvature;
    else
        pili_ext(k).intrinsic_curvature = 'frozen';
    end

    if isfield(pc, 'active_moment') && ~isempty(pc.active_moment)
        pili_ext(k).active_moment = pc.active_moment;
    else
        pili_ext(k).active_moment = struct('type', 'off');
    end
end

end
