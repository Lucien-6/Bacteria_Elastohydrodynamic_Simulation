function out = render_flowfield_panels(frame, opts)
%RENDER_FLOWFIELD_PANELS  Draw the four-panel flow-field figure
%
%  Description:
%    Draws the XZ/YZ flow field on a fixed 1800 x 1400 canvas.  Static
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
%           light direction [-1 1 0]; time is the only figure title

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
        out = update_figure(opts.handles, frame);
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
pos_a = [0.055, 0.56, 0.40, 0.30];
pos_b = [0.50,  0.56, 0.36, 0.30];
pos_c = [0.055, 0.13, 0.40, 0.36];
pos_d = [0.50,  0.13, 0.36, 0.36];

fig = figure('Units', 'pixels', 'Position', [40, 40, 1800, 1400], ...
    'Color', 'w', 'Visible', visible);
set(fig, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 12, ...
    'DefaultAxesLineWidth', 0.6, 'DefaultTextFontName', 'Arial', ...
    'DefaultTextInterpreter', 'none');

ax_a = axes('Parent', fig, 'Position', pos_a);
ax_b = axes('Parent', fig, 'Position', pos_b);
ax_c = axes('Parent', fig, 'Position', pos_c);
ax_d = axes('Parent', fig, 'Position', pos_d);

handles.fig = fig;
handles.ax_a = ax_a;
handles.ax_b = ax_b;
handles.ax_c = ax_c;
handles.ax_d = ax_d;
handles.pos = struct('a', pos_a, 'b', pos_b, 'c', pos_c, 'd', pos_d);
handles.created = false;
handles = update_figure(handles, frame);

cb = colorbar(ax_d);
cb.Location = 'manual';
cb.Position = colorbar_position(ax_d, pos_d);
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


function handles = update_figure(handles, frame)
style_3d(handles.ax_a, frame, 'xz');
style_3d(handles.ax_b, frame, 'yz');
style_2d(handles.ax_c, frame, 'xz');
style_2d(handles.ax_d, frame, 'yz');

handles = paint_3d(handles, frame, 'xz');
handles = paint_3d(handles, frame, 'yz');
handles = paint_2d(handles, frame, 'xz');
handles = paint_2d(handles, frame, 'yz');

if isfield(handles, 'time') && isgraphics(handles.time)
    handles.time.String = time_string(frame.time_s);
end
restore_positions(handles);
handles.created = true;
drawnow;
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


function handles = paint_3d(handles, frame, plane_name)
if strcmp(plane_name, 'xz')
    ax = handles.ax_a;
    tag = 'a';
    plane = frame.planes.xz;
else
    ax = handles.ax_b;
    tag = 'b';
    plane = frame.planes.yz;
end
[C, A] = display_image(plane.Umag, plane.mask);
mesh = build_body_mesh(frame.planes.body_vis, frame.params, ...
    frame.geom.R_body, frame.planes.X_cm_um);
body_color = [0.20, 0.65, 0.33];

if ~handles.created
    hold(ax, 'on');
    handles.(['surf_' tag]) = surface(ax, plane.HG, plane.VG, plane.ZG, C, ...
        'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
        'FaceAlpha', 'texturemap', 'AlphaData', A);
    handles.(['body_' tag]) = render_body_mesh(ax, mesh, body_color, 0.40);
    handles.(['pili_' tag]) = draw_pili_lines(ax, frame);
    handles.(['tag_' tag]) = panel_label(ax, tag);
    % Propagation direction is [-1, 1, 0], so the source sits at [1, -1, 0].
    light(ax, 'Position', [1, -1, 0], 'Style', 'infinite');
    lighting(ax, 'gouraud');
    material(ax, 'dull');
    hold(ax, 'off');
else
    set(handles.(['surf_' tag]), 'XData', plane.HG, 'YData', plane.VG, ...
        'ZData', plane.ZG, 'CData', C, 'AlphaData', A);
    update_body_mesh(handles.(['body_' tag]), mesh);
    update_pili_lines(handles.(['pili_' tag]), frame);
end
apply_3d_limits(ax, frame);
clim(ax, [0, frame.cmax]);
if strcmp(plane_name, 'xz')
    set_panel_title(ax, '3D — XZ plane');
else
    set_panel_title(ax, '3D — YZ plane');
end
end


function handles = paint_2d(handles, frame, plane_name)
if strcmp(plane_name, 'xz')
    ax = handles.ax_c;
    tag = 'c';
    plane = frame.planes.xz;
    Uh = plane.Ux;
    Uv = plane.Uz;
    Gh = plane.HG;
    Gv = plane.ZG;
    vary_col = 2;
    h_col = 1;
    v_col = 3;
    xlab = 'X (μm)';
else
    ax = handles.ax_d;
    tag = 'd';
    plane = frame.planes.yz;
    Uh = plane.Uy;
    Uv = plane.Uz;
    Gh = plane.VG;
    Gv = plane.ZG;
    vary_col = 1;
    h_col = 2;
    v_col = 3;
    xlab = 'Y (μm)';
end
[C, A] = display_image(plane.Umag, plane.mask);
arrow = 1.5 * 0.55 * frame.planes.grid_spacing_um * frame.planes.quiver_skip;
[qx, qy, qu, qv] = quiver_samples(Gh, Gv, Uh, Uv, plane.mask, ...
    frame.planes.quiver_skip, arrow);
[cx, cy] = pili_crossings(frame, vary_col, h_col, v_col);

if ~handles.created
    hold(ax, 'on');
    handles.(['img_' tag]) = imagesc(ax, plane.h_vec, plane.v_vec, C);
    set(handles.(['img_' tag]), 'AlphaData', A);
    set(ax, 'YDir', 'normal');
    handles.(['outline_' tag]) = fill(ax, [0, 1, 0], [0, 0, 1], 'w', ...
        'EdgeColor', [0.25, 0.25, 0.25], 'LineWidth', 1.0);
    handles.(['cross_' tag]) = plot(ax, nan, nan, 'o', 'MarkerSize', 5, ...
        'LineWidth', 1.0, 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', [0.25, 0.25, 0.25]);
    handles.(['quiver_' tag]) = quiver(ax, qx, qy, qu, qv, 0, ...
        'Color', [0.12, 0.12, 0.12], 'LineWidth', 0.75, ...
        'MaxHeadSize', 0.45, 'AutoScale', 'off');
    handles.(['tag_' tag]) = panel_label(ax, tag);
    xlabel(ax, xlab, 'FontName', 'Arial', 'FontSize', 14, 'Interpreter', 'none');
    ylabel(ax, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 14, 'Interpreter', 'none');
    hold(ax, 'off');
else
    set(handles.(['img_' tag]), 'XData', [plane.h_vec(1), plane.h_vec(end)], ...
        'YData', [plane.v_vec(1), plane.v_vec(end)], ...
        'CData', C, 'AlphaData', A);
    delete(handles.(['quiver_' tag]));
    hold(ax, 'on');
    handles.(['quiver_' tag]) = quiver(ax, qx, qy, qu, qv, 0, ...
        'Color', [0.12, 0.12, 0.12], 'LineWidth', 0.75, ...
        'MaxHeadSize', 0.45, 'AutoScale', 'off');
    hold(ax, 'off');
end
set_xy(handles.(['outline_' tag]), plane.cap_h, plane.cap_v);
set_xy(handles.(['cross_' tag]), cx, cy);
apply_2d_limits(ax, frame, plane_name);
clim(ax, [0, frame.cmax]);
colormap(ax, jet(256));
if strcmp(plane_name, 'xz')
    set_panel_title(ax, 'XZ plane');
else
    set_panel_title(ax, 'YZ plane');
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
if strcmp(plane_name, 'xz')
    set_panel_title(ax, '3D — XZ plane');
else
    set_panel_title(ax, '3D — YZ plane');
end
end


function style_2d(ax, frame, plane_name)
box(ax, 'on');
set(ax, 'TickDir', 'out', 'Layer', 'top', 'FontName', 'Arial', 'FontSize', 12, ...
    'LabelFontSizeMultiplier', 1);
colormap(ax, jet(256));
clim(ax, [0, frame.cmax]);
if strcmp(plane_name, 'xz')
    set_panel_title(ax, 'XZ plane');
else
    set_panel_title(ax, 'YZ plane');
end
end


function apply_3d_limits(ax, frame)
lim = frame.planes.axis_lim_um;
xlim(ax, lim(1, :));
ylim(ax, lim(2, :));
zlim(ax, lim(3, :));
set(ax, 'DataAspectRatio', [1, 1, 1], ...
    'PlotBoxAspectRatio', [diff(lim(1, :)), diff(lim(2, :)), diff(lim(3, :))]);
end


function apply_2d_limits(ax, frame, plane_name)
lim = frame.planes.axis_lim_um;
if strcmp(plane_name, 'xz')
    xr = lim(1, :);
    yr = lim(3, :);
else
    xr = lim(2, :);
    yr = lim(3, :);
end
xlim(ax, xr);
ylim(ax, yr);
set(ax, 'DataAspectRatio', [1, 1, 1], ...
    'PlotBoxAspectRatio', [diff(xr), diff(yr), 1]);
end


function lines = draw_pili_lines(ax, frame)
colors = pilus_colors();
M = frame.params.M_pili;
lines = gobjects(M, 1);
hold(ax, 'on');
for k = 1:M
    Xp = frame.geom.X_pili{k} * frame.L_to_um;
    pk = colors(mod(k - 1, size(colors, 1)) + 1, :);
    lines(k) = plot3(ax, Xp(:, 1), Xp(:, 2), Xp(:, 3), '-', ...
        'Color', pk, 'LineWidth', 1.0);
end
end


function update_pili_lines(lines, frame)
for k = 1:numel(lines)
    Xp = frame.geom.X_pili{k} * frame.L_to_um;
    set(lines(k), 'XData', Xp(:, 1), 'YData', Xp(:, 2), 'ZData', Xp(:, 3));
end
end


function [C, A] = display_image(Umag, mask)
C = Umag;
C(mask) = NaN;
A = double(~mask);
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
%COLORBAR_POSITION  Vertical colour bar just to the right of panel (d)
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


function [qx, qy, qu, qv] = quiver_samples(Gh, Gv, Uh, Uv, mask, skip, arrow)
qh = 1:skip:size(Uh, 2);
qv_idx = 1:skip:size(Uh, 1);
Uh_s = Uh(qv_idx, qh);
Uv_s = Uv(qv_idx, qh);
Gh_s = Gh(qv_idx, qh);
Gv_s = Gv(qv_idx, qh);
mk = mask(qv_idx, qh);
spd = sqrt(Uh_s.^2 + Uv_s.^2);
scale = max(spd(:), [], 'omitnan');
ok = isfinite(spd) & spd > 1e-15 * scale & ~mk;
qx = Gh_s(ok);
qy = Gv_s(ok);
if isempty(qx)
    qx = nan; qy = nan; qu = nan; qv = nan;
    return
end
qu = Uh_s(ok) ./ spd(ok) * arrow;
qv = Uv_s(ok) ./ spd(ok) * arrow;
end


function [cx, cy] = pili_crossings(frame, vary_col, h_col, v_col)
cm = frame.planes.X_cm_um(:);
cx = [];
cy = [];
for k = 1:frame.params.M_pili
    Xp = frame.geom.X_pili{k} * frame.L_to_um;
    d = Xp(:, vary_col) - cm(vary_col);
    for j = 1:size(Xp, 1) - 1
        if d(j) * d(j + 1) < 0
            t = d(j) / (d(j) - d(j + 1));
            cx(end + 1) = Xp(j, h_col) + t * (Xp(j + 1, h_col) - Xp(j, h_col)); %#ok<AGROW>
            cy(end + 1) = Xp(j, v_col) + t * (Xp(j + 1, v_col) - Xp(j, v_col)); %#ok<AGROW>
        end
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
if isfield(handles, 'cb') && isgraphics(handles.cb)
    handles.cb.Location = 'manual';
    handles.cb.Position = colorbar_position(handles.ax_d, handles.pos.d);
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


function body_mesh = build_body_mesh(body_vis, params, R_body, X_cm_um)
body_mesh.type = body_vis.type;
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
        mesh_size = size(capsule_X0);
        capsule_flat = [capsule_X0(:).'; capsule_Y0(:).'; capsule_Z0(:).'];
        pts_3d = R_body * capsule_flat + X_cm_um(:);
        body_mesh.Xsurf = reshape(pts_3d(1, :), mesh_size);
        body_mesh.Ysurf = reshape(pts_3d(2, :), mesh_size);
        body_mesh.Zsurf = reshape(pts_3d(3, :), mesh_size);
    case 'custom'
        verts_um = body_vis.verts_m * 1e6;
        pts_lab = R_body * verts_um.' + X_cm_um(:);
        body_mesh.X = pts_lab(1, :).';
        body_mesh.Y = pts_lab(2, :).';
        body_mesh.Z = pts_lab(3, :).';
        body_mesh.faces = body_vis.faces;
    otherwise
        error('render_flowfield_panels:UnknownBodyType', ...
            'Unknown body_vis.type = ''%s''.', body_vis.type);
end
assert(isstruct(params));
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
