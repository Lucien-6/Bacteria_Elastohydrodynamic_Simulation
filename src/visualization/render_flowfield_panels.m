function out = render_flowfield_panels(frame, opts)
%RENDER_FLOWFIELD_PANELS  Draw the six-panel flow-field figure
%
%  Description:
%    Draws the XZ/YZ/XY flow field on a fixed 2500 x 1400 canvas.  The top
%    row holds the 3D views (a-c) and the bottom row the 2D views (d-f).  Static
%    calls export PNG, or the existing PNG/PDF/HTML bundle.  Animation
%    calls create the objects once and then update them.
%
%  Syntax:
%    render_flowfield_panels(frame, opts)
%    handles = render_flowfield_panels(frame, struct('action','create'))
%    handles = render_flowfield_panels(frame, struct('action','update','handles',handles))
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28
%  Updated: 2026-10-30 — body interior left blank; axes cover the cell;
%           panel letters outside the axes; vertical colour bar beside (d);
%           light direction [-1 1 -1]; time is the only figure title
%  Updated: 2026-11-11 — added the XY plane; 2 x 3 layout, colour bar beside (f)
%  Updated: 2026-11-13 — 2D panels draw in-plane streamlines instead of arrows
%  Updated: 2026-11-13 — no white halo; each curve runs to the panel edge
%  Updated: 2026-11-24 — one-time static styling, cached local body mesh &
%           per-frame mesh/pili reuse, 2x2 texturemap surface corners,
%           precomputed streamline reuse, optional skip_drawnow

if nargin < 2 || isempty(opts), opts = struct(); end
action = 'static';
if isfield(opts, 'action') && ~isempty(opts.action)
    action = char(opts.action);
end

switch action
    case 'static'
        visible = opt_value(opts, 'visible', 'on');
        handles = create_figure(frame, visible);
        export_static(handles.fig, opts);
        if opt_value(opts, 'close_after', true)
            close(handles.fig);
            out = [];
        else
            out = handles;
        end
    case 'create'
        out = create_figure(frame, opt_value(opts, 'visible', 'off'));
    case 'update'
        skip_drawnow = logical(opt_value(opts, 'skip_drawnow', false));
        out = update_figure(opts.handles, frame, skip_drawnow);
    case 'export_png'
        exportgraphics(opts.figure, opts.png_path, ...
            'Resolution', opt_value(opts, 'dpi', 600));
        out = [];
    otherwise
        error('render_flowfield_panels:BadAction', ...
            'Unknown action ''%s''.', action);
end

end


function handles = create_figure(frame, visible)
% Top row: 3D views (a) XZ, (b) YZ, (c) XY.
% Bottom row: 2D views (d) XZ, (e) YZ, (f) XY.
pos_a = [0.040, 0.56, 0.28, 0.30];
pos_b = [0.355, 0.56, 0.26, 0.30];
pos_c = [0.660, 0.56, 0.26, 0.30];
pos_d = [0.040, 0.13, 0.28, 0.36];
pos_e = [0.355, 0.13, 0.26, 0.36];
pos_f = [0.660, 0.13, 0.26, 0.36];

fig = figure('Units', 'pixels', 'Position', [40, 40, 2500, 1400], ...
    'Color', 'w', 'Visible', visible);
set(fig, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 12, ...
    'DefaultAxesLineWidth', 0.6, 'DefaultTextFontName', 'Arial', ...
    'DefaultTextInterpreter', 'none');

ax_a = axes('Parent', fig, 'Position', pos_a);
ax_b = axes('Parent', fig, 'Position', pos_b);
ax_c = axes('Parent', fig, 'Position', pos_c);
ax_d = axes('Parent', fig, 'Position', pos_d);
ax_e = axes('Parent', fig, 'Position', pos_e);
ax_f = axes('Parent', fig, 'Position', pos_f);

handles.fig = fig;
handles.ax_a = ax_a;
handles.ax_b = ax_b;
handles.ax_c = ax_c;
handles.ax_d = ax_d;
handles.ax_e = ax_e;
handles.ax_f = ax_f;
handles.pos = struct('a', pos_a, 'b', pos_b, 'c', pos_c, ...
    'd', pos_d, 'e', pos_e, 'f', pos_f);
handles.created = false;
handles.cmax = NaN;
handles.last_spans = [];
handles.body_local_mesh = build_body_local_mesh(frame.planes.body_vis);

handles = update_figure(handles, frame, false);

cb = colorbar(ax_f);
cb.Location = 'manual';
cb.Position = colorbar_position(ax_f, pos_f);
cb.FontName = 'Arial';
cb.FontSize = 12;
cb.Label.String = 'Velocity magnitude (μm/s)';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 14;
cb.Label.Rotation = 90;
cb.Label.Interpreter = 'none';
cb.TickLabelInterpreter = 'none';
handles.cb = cb;
restore_positions(handles);

handles.time = annotation(fig, 'textbox', [0.08, 0.909, 0.84, 0.045], ...
    'String', time_string(frame.time_s), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontName', 'Arial', ...
    'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none');
restore_positions(handles);
drawnow;
end


function handles = update_figure(handles, frame, skip_drawnow)
plane_names = {'xz', 'yz', 'xy'};
is_first = ~handles.created;

if is_first
    for n = 1:numel(plane_names)
        spec = plane_spec(plane_names{n});
        style_3d(handles.(spec.ax3), frame, plane_names{n});
        style_2d(handles.(spec.ax2), frame, plane_names{n});
    end
end

% Transform 3D body surface and pili coordinates ONCE per frame
if ~isfield(handles, 'body_local_mesh') || isempty(handles.body_local_mesh)
    handles.body_local_mesh = build_body_local_mesh(frame.planes.body_vis);
end
body_mesh = transform_body_mesh(handles.body_local_mesh, ...
    frame.geom.R_body, frame.planes.X_cm_um);
M = frame.params.M_pili;
pili_um = cell(M, 1);
for k = 1:M
    pili_um{k} = frame.geom.X_pili{k} * frame.L_to_um;
end

lim = frame.planes.axis_lim_um;
spans = diff(lim, 1, 2);
spans_changed = is_first || isempty(handles.last_spans) ...
    || any(abs(spans - handles.last_spans) > 1e-12);
cmax_changed = is_first || ~isequal(handles.cmax, frame.cmax);

for n = 1:numel(plane_names)
    handles = paint_3d(handles, frame, plane_names{n}, body_mesh, pili_um, ...
        spans_changed, cmax_changed);
end
for n = 1:numel(plane_names)
    handles = paint_2d(handles, frame, plane_names{n}, pili_um, ...
        spans_changed, cmax_changed);
end

if isfield(handles, 'time') && isgraphics(handles.time)
    handles.time.String = time_string(frame.time_s);
end
if spans_changed
    restore_positions(handles);
end
handles.last_spans = spans;
handles.cmax = frame.cmax;
handles.created = true;
if ~skip_drawnow
    drawnow;
end
end


function export_static(fig, opts)
dpi = opt_value(opts, 'dpi', 600);
bundle = opt_value(opts, 'bundle_export', false);
if bundle
    export_figure_png_pdf_html(fig, opts.png_path, opts.pdf_path, opts.html_path);
elseif isfield(opts, 'png_path') && ~isempty(opts.png_path)
    exportgraphics(fig, opts.png_path, 'Resolution', dpi);
end
end


function spec = plane_spec(plane_name)
%PLANE_SPEC  Axes handles, panel letters, and labels of one plane
switch plane_name
    case 'xz'
        spec.ax3 = 'ax_a';  spec.tag3 = 'a';
        spec.ax2 = 'ax_d';  spec.tag2 = 'd';
        spec.xlab = 'X (μm)';  spec.ylab = 'Z (μm)';
    case 'yz'
        spec.ax3 = 'ax_b';  spec.tag3 = 'b';
        spec.ax2 = 'ax_e';  spec.tag2 = 'e';
        spec.xlab = 'Y (μm)';  spec.ylab = 'Z (μm)';
    case 'xy'
        spec.ax3 = 'ax_c';  spec.tag3 = 'c';
        spec.ax2 = 'ax_f';  spec.tag2 = 'f';
        spec.xlab = 'X (μm)';  spec.ylab = 'Y (μm)';
    otherwise
        error('render_flowfield_panels:BadPlane', ...
            'Unknown plane ''%s''.', plane_name);
end
spec.title2 = [upper(plane_name), ' plane'];
spec.title3 = ['3D — ', upper(plane_name), ' plane'];
end


function handles = paint_3d(handles, frame, plane_name, body_mesh, pili_um, ...
    spans_changed, cmax_changed)
spec = plane_spec(plane_name);
ax = handles.(spec.ax3);
tag = spec.tag3;
plane = frame.planes.(plane_name);
[C, A] = display_image(plane.Umag, plane.mask);
[Xs, Ys, Zs] = plane_surface_corners(plane, frame.planes.X_cm_um, plane_name);
body_color = [0.20, 0.65, 0.33];

if ~handles.created
    hold(ax, 'on');
    handles.(['surf_' tag]) = surface(ax, Xs, Ys, Zs, C, ...
        'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
        'FaceAlpha', 'texturemap', 'AlphaData', A);
    handles.(['body_' tag]) = render_body_mesh(ax, body_mesh, body_color, 0.40);
    handles.(['pili_' tag]) = draw_pili_lines(ax, pili_um);
    handles.(['tag_' tag]) = panel_label(ax, tag);
    % Propagation direction is [-1, 1, -1], so the source sits at [1, -1, 1].
    light(ax, 'Position', [1, -1, 1], 'Style', 'infinite');
    lighting(ax, 'gouraud');
    material(ax, 'dull');
    hold(ax, 'off');
else
    set(handles.(['surf_' tag]), 'XData', Xs, 'YData', Ys, ...
        'ZData', Zs, 'CData', C, 'AlphaData', A);
    update_body_mesh(handles.(['body_' tag]), body_mesh);
    update_pili_lines(handles.(['pili_' tag]), pili_um);
end
apply_3d_limits(ax, frame, spans_changed);
if cmax_changed
    clim(ax, [0, frame.cmax]);
end
end


function handles = paint_2d(handles, frame, plane_name, pili_um, ...
    spans_changed, cmax_changed)
spec = plane_spec(plane_name);
ax = handles.(spec.ax2);
tag = spec.tag2;
plane = frame.planes.(plane_name);
switch plane_name
    case 'xz'
        vary_col = 2;  h_col = 1;  v_col = 3;
    case 'yz'
        vary_col = 1;  h_col = 2;  v_col = 3;
    case 'xy'
        vary_col = 3;  h_col = 1;  v_col = 2;
end
xlab = spec.xlab;
ylab = spec.ylab;
[C, A] = display_image(plane.Umag, plane.mask);
if isfield(plane, 'stream') && isstruct(plane.stream) && ~isempty(plane.stream)
    geom = plane.stream;
else
    switch plane_name
        case 'xz'
            Uh = plane.Ux;  Uv = plane.Uz;
        case 'yz'
            Uh = plane.Uy;  Uv = plane.Uz;
        case 'xy'
            Uh = plane.Ux;  Uv = plane.Uy;
    end
    geom = plane_stream_geometry(frame, plane, plane_name, Uh, Uv);
end
[cx, cy] = pili_crossings(pili_um, frame.planes.X_cm_um, vary_col, h_col, v_col);
ink = [0.12, 0.12, 0.12];

if ~handles.created
    hold(ax, 'on');
    handles.(['img_' tag]) = imagesc(ax, plane.h_vec, plane.v_vec, C);
    set(handles.(['img_' tag]), 'AlphaData', A);
    set(ax, 'YDir', 'normal');
    % Curve, arrowhead, then the body so the outline stays visible where a
    % curve stops on the mask.
    handles.(['line_k_' tag]) = plot(ax, geom.x, geom.y, '-', ...
        'Color', ink, 'LineWidth', 0.75, 'Clipping', 'on');
    handles.(['head_' tag]) = patch(ax, 'Faces', geom.faces, ...
        'Vertices', geom.vertices, 'FaceColor', ink, 'EdgeColor', 'none', ...
        'Clipping', 'on');
    handles.(['outline_' tag]) = fill(ax, [0, 1, 0], [0, 0, 1], 'w', ...
        'EdgeColor', [0.25, 0.25, 0.25], 'LineWidth', 1.0);
    handles.(['cross_' tag]) = plot(ax, nan, nan, 'o', 'MarkerSize', 5, ...
        'LineWidth', 1.0, 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', [0.25, 0.25, 0.25]);
    handles.(['tag_' tag]) = panel_label(ax, tag);
    xlabel(ax, xlab, 'FontName', 'Arial', 'FontSize', 14, 'Interpreter', 'none');
    ylabel(ax, ylab, 'FontName', 'Arial', 'FontSize', 14, 'Interpreter', 'none');
    hold(ax, 'off');
else
    set(handles.(['img_' tag]), 'XData', [plane.h_vec(1), plane.h_vec(end)], ...
        'YData', [plane.v_vec(1), plane.v_vec(end)], ...
        'CData', C, 'AlphaData', A);
    set(handles.(['line_k_' tag]), 'XData', geom.x, 'YData', geom.y);
    set(handles.(['head_' tag]), 'Faces', geom.faces, ...
        'Vertices', geom.vertices);
end
set_xy(handles.(['outline_' tag]), plane.cap_h, plane.cap_v);
set_xy(handles.(['cross_' tag]), cx, cy);
apply_2d_limits(ax, frame, plane_name, spans_changed);
if cmax_changed
    clim(ax, [0, frame.cmax]);
end
end


function [Xs, Ys, Zs] = plane_surface_corners(plane, X_cm_um, plane_name)
%PLANE_SURFACE_CORNERS  2x2 corner grid for texture-mapped 3D plane surface
if isfield(plane, 'HG') && ~isempty(plane.HG) && ...
        size(plane.HG, 1) >= 2 && size(plane.HG, 2) >= 2
    Xs = plane.HG([1, end], [1, end]);
    Ys = plane.VG([1, end], [1, end]);
    Zs = plane.ZG([1, end], [1, end]);
    return
end
h0 = plane.h_vec(1);  h1 = plane.h_vec(end);
v0 = plane.v_vec(1);  v1 = plane.v_vec(end);
cm = X_cm_um(:);
switch plane_name
    case 'xz'
        Xs = [h0, h1; h0, h1];
        Ys = [cm(2), cm(2); cm(2), cm(2)];
        Zs = [v0, v0; v1, v1];
    case 'yz'
        Xs = [cm(1), cm(1); cm(1), cm(1)];
        Ys = [h0, h1; h0, h1];
        Zs = [v0, v0; v1, v1];
    case 'xy'
        Xs = [h0, h1; h0, h1];
        Ys = [v0, v0; v1, v1];
        Zs = [cm(3), cm(3); cm(3), cm(3)];
end
end


function style_3d(ax, frame, plane_name)
grid(ax, 'on');
box(ax, 'on');
set(ax, 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 12, ...
    'LabelFontSizeMultiplier', 1);
view(ax, 45, 25);
camproj(ax, 'orthographic');
colormap(ax, jet(256));
clim(ax, [0, frame.cmax]);
xlabel(ax, 'X (μm)', 'FontName', 'Arial', 'FontSize', 14, 'Interpreter', 'none');
ylabel(ax, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 14, 'Interpreter', 'none');
zlabel(ax, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 14, 'Interpreter', 'none');
spec = plane_spec(plane_name);
set_panel_title(ax, spec.title3);
end


function style_2d(ax, frame, plane_name)
box(ax, 'on');
set(ax, 'TickDir', 'out', 'Layer', 'top', 'FontName', 'Arial', 'FontSize', 12, ...
    'LabelFontSizeMultiplier', 1);
colormap(ax, jet(256));
clim(ax, [0, frame.cmax]);
spec = plane_spec(plane_name);
set_panel_title(ax, spec.title2);
end


function apply_3d_limits(ax, frame, spans_changed)
lim = frame.planes.axis_lim_um;
if spans_changed
    set(ax, 'XLim', lim(1, :), 'YLim', lim(2, :), 'ZLim', lim(3, :), ...
        'DataAspectRatio', [1, 1, 1], ...
        'PlotBoxAspectRatio', [diff(lim(1, :)), diff(lim(2, :)), diff(lim(3, :))]);
else
    set(ax, 'XLim', lim(1, :), 'YLim', lim(2, :), 'ZLim', lim(3, :));
end
end


function apply_2d_limits(ax, frame, plane_name, spans_changed)
lim = frame.planes.axis_lim_um;
switch plane_name
    case 'xz'
        xr = lim(1, :);  yr = lim(3, :);
    case 'yz'
        xr = lim(2, :);  yr = lim(3, :);
    case 'xy'
        xr = lim(1, :);  yr = lim(2, :);
end
if spans_changed
    set(ax, 'XLim', xr, 'YLim', yr, ...
        'DataAspectRatio', [1, 1, 1], ...
        'PlotBoxAspectRatio', [diff(xr), diff(yr), 1]);
else
    set(ax, 'XLim', xr, 'YLim', yr);
end
end


function lines = draw_pili_lines(ax, pili_um)
colors = pilus_colors();
M = numel(pili_um);
lines = gobjects(M, 1);
hold(ax, 'on');
for k = 1:M
    Xp = pili_um{k};
    pk = colors(mod(k - 1, size(colors, 1)) + 1, :);
    lines(k) = plot3(ax, Xp(:, 1), Xp(:, 2), Xp(:, 3), '-', ...
        'Color', pk, 'LineWidth', 1.0);
end
end


function update_pili_lines(lines, pili_um)
for k = 1:numel(lines)
    Xp = pili_um{k};
    set(lines(k), 'XData', Xp(:, 1), 'YData', Xp(:, 2), 'ZData', Xp(:, 3));
end
end


function [C, A] = display_image(Umag, mask)
C = Umag;
C(mask) = NaN;
A = single(~mask);
end


function h = panel_label(ax, tag)
%PANEL_LABEL  Bold panel letter just outside the top-left corner of the axes
h = text(ax, 0, 1.02, tag, 'Units', 'normalized', ...
    'Clipping', 'off', 'FontName', 'Arial', 'FontSize', 18, ...
    'FontWeight', 'bold', 'Interpreter', 'none', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
    'HitTest', 'off');
end


function set_panel_title(ax, str)
%SET_PANEL_TITLE  Subplot title in Arial, clear of the top of the axes box
h = title(ax, str, 'FontName', 'Arial', 'FontSize', 15, ...
    'FontWeight', 'bold', 'Interpreter', 'none');
set(h, 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Position', [0.5, 1.04, 0]);
end


function p = colorbar_position(ax, pos_d)
%COLORBAR_POSITION  Vertical colour bar just to the right of panel (f)
%  Equal scaling leaves the square plot box inside a wider axes rectangle.
%  The bar is placed against that plot box, not against the axes rectangle.
fig = ancestor(ax, 'figure');
fig_units = fig.Units;
fig.Units = 'pixels';
fp = fig.Position;
fig.Units = fig_units;

pbar = ax.PlotBoxAspectRatio;
target = pbar(1) / max(pbar(2), eps);
aw = pos_d(3) * fp(3);
ah = pos_d(4) * fp(4);
if aw / max(ah, eps) > target
    bw = ah * target;
    bh = ah;
else
    bw = aw;
    bh = aw / max(target, eps);
end
bx = pos_d(1) * fp(3) + (aw - bw) / 2;
by = pos_d(2) * fp(4) + (ah - bh) / 2;
% A little wider than before: about five times the axes frame line.
ppi = get(groot, 'ScreenPixelsPerInch');
if ~(isfinite(ppi) && ppi > 0)
    ppi = 96;
end
gap_px = 5 * ax.LineWidth * ppi / 72;
if gap_px < 0.008 * fp(3)
    gap_px = 0.008 * fp(3);
end
bar_w = 0.016 * fp(3);
p = [(bx + bw + gap_px) / fp(3), by / fp(4), bar_w / fp(3), bh / fp(4)];
end


function geom = plane_stream_geometry(frame, plane, plane_name, Uh, Uv)
%PLANE_STREAM_GEOMETRY  Resolve spacing and length, then place the curves
span_h = abs(plane.h_vec(end) - plane.h_vec(1));
span_v = abs(plane.v_vec(end) - plane.v_vec(1));
spacing_um = 0.045 * min(span_h, span_v);
length_um = 4 * hypot(span_h, span_v);
if isfield(frame.planes, 'stream_spacing_um') && ~isempty(frame.planes.stream_spacing_um)
    spacing_um = frame.planes.stream_spacing_um;
end
if isfield(frame.planes, 'stream_length_um') && ~isempty(frame.planes.stream_length_um)
    length_um = frame.planes.stream_length_um;
end
cm = frame.planes.X_cm_um(:);
switch plane_name
    case 'xz'
        origin = [cm(1); cm(3)];
    case 'yz'
        origin = [cm(2); cm(3)];
    case 'xy'
        origin = [cm(1); cm(2)];
    otherwise
        error('render_flowfield_panels:BadPlane', ...
            'Unknown plane ''%s''.', plane_name);
end
geom = flowfield_plane_streamlines(plane.h_vec, plane.v_vec, Uh, Uv, ...
    plane.mask, origin, spacing_um, length_um, plane_name, true);
end


function [cx, cy] = pili_crossings(pili_um, X_cm_um, vary_col, h_col, v_col)
cm = X_cm_um(:);
cx = [];
cy = [];
for k = 1:numel(pili_um)
    Xp = pili_um{k};
    d = Xp(:, vary_col) - cm(vary_col);
    hit = find(d(1:end - 1) .* d(2:end) < 0);
    if ~isempty(hit)
        t = d(hit) ./ (d(hit) - d(hit + 1));
        cx = [cx; Xp(hit, h_col) + t .* (Xp(hit + 1, h_col) - Xp(hit, h_col))]; %#ok<AGROW>
        cy = [cy; Xp(hit, v_col) + t .* (Xp(hit + 1, v_col) - Xp(hit, v_col))]; %#ok<AGROW>
    end
end
if isempty(cx)
    cx = nan;
    cy = nan;
end
end


function set_xy(h, x, y)
if isempty(x)
    set(h, 'XData', nan, 'YData', nan);
else
    set(h, 'XData', x(:)', 'YData', y(:)');
end
end


function restore_positions(handles)
if ~isfield(handles, 'pos'), return; end
handles.ax_a.Position = handles.pos.a;
handles.ax_b.Position = handles.pos.b;
handles.ax_c.Position = handles.pos.c;
handles.ax_d.Position = handles.pos.d;
handles.ax_e.Position = handles.pos.e;
handles.ax_f.Position = handles.pos.f;
if isfield(handles, 'cb') && isgraphics(handles.cb)
    handles.cb.Location = 'manual';
    handles.cb.Position = colorbar_position(handles.ax_f, handles.pos.f);
end
end


function s = time_string(time_s)
s = sprintf('Time = %.4f s', time_s);
end


function colors = pilus_colors()
colors = [0.20 0.40 0.80;
          0.85 0.33 0.10;
          0.47 0.67 0.19;
          0.93 0.69 0.13;
          0.49 0.18 0.56;
          0.30 0.75 0.93;
          0.64 0.08 0.18;
          0.00 0.45 0.74;
          0.80 0.40 0.60;
          0.50 0.50 0.50];
end


function value = opt_value(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end


function local_mesh = build_body_local_mesh(body_vis)
%BUILD_BODY_LOCAL_MESH  Precompute body-local 3D surface mesh once per figure
local_mesh.type = body_vis.type;
switch body_vis.type
    case 'capsule'
        a_um = body_vis.a_um;
        b_um = body_vis.b_um;
        n_circ = 36;
        n_cap = 12;
        n_cyl = 10;
        theta_cap = linspace(0, 2 * pi, n_circ);
        phi_L = linspace(pi, pi / 2, n_cap);
        [TH_L, PH_L] = meshgrid(theta_cap, phi_L);
        XL = b_um * cos(PH_L) - (a_um - b_um);
        YL = b_um * sin(PH_L) .* cos(TH_L);
        ZL = b_um * sin(PH_L) .* sin(TH_L);
        x_cyl = linspace(-(a_um - b_um), (a_um - b_um), n_cyl);
        [TH_C, XC_grid] = meshgrid(theta_cap, x_cyl);
        YC = b_um * cos(TH_C);
        ZC = b_um * sin(TH_C);
        phi_R = linspace(pi / 2, 0, n_cap);
        [TH_R, PH_R] = meshgrid(theta_cap, phi_R);
        XR = b_um * cos(PH_R) + (a_um - b_um);
        YR = b_um * sin(PH_R) .* cos(TH_R);
        ZR = b_um * sin(PH_R) .* sin(TH_R);
        capsule_X0 = [XL; XC_grid; XR];
        capsule_Y0 = [YL; YC; YR];
        capsule_Z0 = [ZL; ZC; ZR];
        local_mesh.mesh_size = size(capsule_X0);
        local_mesh.flat = [capsule_X0(:).'; capsule_Y0(:).'; capsule_Z0(:).'];
    case 'custom'
        local_mesh.verts_um_t = (body_vis.verts_m * 1e6).';
        local_mesh.faces = body_vis.faces;
    otherwise
        error('render_flowfield_panels:UnknownBodyType', ...
            'Unknown body_vis.type = ''%s''.', body_vis.type);
end
end


function body_mesh = transform_body_mesh(local_mesh, R_body, X_cm_um)
body_mesh.type = local_mesh.type;
switch local_mesh.type
    case 'capsule'
        pts_3d = R_body * local_mesh.flat + X_cm_um(:);
        ms = local_mesh.mesh_size;
        body_mesh.Xsurf = reshape(pts_3d(1, :), ms);
        body_mesh.Ysurf = reshape(pts_3d(2, :), ms);
        body_mesh.Zsurf = reshape(pts_3d(3, :), ms);
    case 'custom'
        pts_lab = R_body * local_mesh.verts_um_t + X_cm_um(:);
        body_mesh.X = pts_lab(1, :).';
        body_mesh.Y = pts_lab(2, :).';
        body_mesh.Z = pts_lab(3, :).';
        body_mesh.faces = local_mesh.faces;
    otherwise
        error('render_flowfield_panels:UnknownBodyType', ...
            'Unknown body_vis.type = ''%s''.', local_mesh.type);
end
end


function h = render_body_mesh(ax, body_mesh, body_color, body_alpha)
switch body_mesh.type
    case 'capsule'
        h = surf(ax, body_mesh.Xsurf, body_mesh.Ysurf, body_mesh.Zsurf, ...
            'FaceColor', body_color, 'FaceAlpha', body_alpha, ...
            'EdgeColor', 'none', 'FaceLighting', 'gouraud');
    case 'custom'
        if isempty(body_mesh.faces) || size(body_mesh.faces, 2) ~= 3
            h = scatter3(ax, body_mesh.X, body_mesh.Y, body_mesh.Z, ...
                8, body_color, 'filled', 'MarkerFaceAlpha', body_alpha);
        else
            h = trisurf(body_mesh.faces, body_mesh.X, body_mesh.Y, body_mesh.Z, ...
                'FaceColor', body_color, 'FaceAlpha', body_alpha, ...
                'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'Parent', ax);
        end
    otherwise
        error('render_flowfield_panels:UnknownBodyType', ...
            'Unknown body_mesh.type = ''%s''.', body_mesh.type);
end
end


function update_body_mesh(h, body_mesh)
switch body_mesh.type
    case 'capsule'
        set(h, 'XData', body_mesh.Xsurf, 'YData', body_mesh.Ysurf, ...
            'ZData', body_mesh.Zsurf);
    case 'custom'
        if isprop(h, 'Vertices')
            set(h, 'Vertices', [body_mesh.X, body_mesh.Y, body_mesh.Z]);
        else
            set(h, 'XData', body_mesh.X, 'YData', body_mesh.Y, 'ZData', body_mesh.Z);
        end
    otherwise
        error('render_flowfield_panels:UnknownBodyType', ...
            'Unknown body_mesh.type = ''%s''.', body_mesh.type);
end
end
