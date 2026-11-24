function [lo_um, hi_um, cm_um] = flowfield_frame_window(geom, params, flow_opts)
%FLOWFIELD_FRAME_WINDOW  One-frame axis window around the centre of mass
%
%  Description:
%    The configured cube of side grid_range_um is centred on the body
%    centre of mass.  Each side is then pushed out so the body surface and
%    every flagellum node lie inside, with a margin of at least 0.5 um.
%
%  Syntax:
%    [lo_um, hi_um, cm_um] = flowfield_frame_window(geom, params, flow_opts)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-13

if nargin < 3 || isempty(flow_opts), flow_opts = struct(); end
spacing_um = opt_num(flow_opts, 'grid_spacing_um', 0.05);
grid_range_um = opt_num(flow_opts, 'grid_range_um', 10.0);
half_range = grid_range_um / 2;
L_to_um = params.scales.L_0 * 1e6;

cm = geom.X_cm(:);
cm_um = cm * L_to_um;

n_body = size(geom.X_body, 1);
n_pili = 0;
for k = 1:numel(geom.X_pili)
    n_pili = n_pili + size(geom.X_pili{k}, 1);
end
use_mesh = isfield(params, 'body_type') && strcmp(params.body_type, 'custom') ...
    && isfield(params, 'body_verts_local') && ~isempty(params.body_verts_local);
n_mesh = 0;
if use_mesh
    n_mesh = size(params.body_verts_local, 1);
end
pts = zeros(n_body + n_pili + n_mesh, 3);
pts(1:n_body, :) = geom.X_body;
row = n_body + 1;
for k = 1:numel(geom.X_pili)
    nk = size(geom.X_pili{k}, 1);
    pts(row:row + nk - 1, :) = geom.X_pili{k};
    row = row + nk;
end
if use_mesh
    pts(row:end, :) = (geom.R_body * params.body_verts_local.' + cm).';
end

pts_um = pts * L_to_um;
bb_min = min(pts_um, [], 1);
bb_max = max(pts_um, [], 1);
pad = max(0.5, 4 * spacing_um);
lo_um = min(cm_um.' - half_range, bb_min - pad);
hi_um = max(cm_um.' + half_range, bb_max + pad);
end


function v = opt_num(s, name, default_val)
v = default_val;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
end
end
