function check_flowfield_streamlines()
%CHECK_FLOWFIELD_STREAMLINES  Acceptance checks for planar streamlines
%
%  Description:
%    Checks option parsing, curve placement (direction, length, separation,
%    one closed orbit, body exclusion), and the 2D panel drawing (object
%    order, no quiver marks, update in place).  Writes a preview PNG to the
%    temp directory and prints its path.  Does not need a saved simulation.
%
%  Syntax:
%    check_flowfield_streamlines()
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-13

this_dir = fileparts(mfilename('fullpath'));
proj = fileparts(fileparts(this_dir));
addpath(this_dir);
addpath(fullfile(proj, 'src', 'config'));

fprintf('1. Streamline option checks\n');
check_overrides();
fprintf('2. Uniform flow: length, direction, separation\n');
check_uniform();
fprintf('3. Closed orbit is drawn once\n');
check_vortex();
fprintf('4. Curves stay outside the body\n');
check_mask();
fprintf('5. Rejected inputs\n');
check_rejects();
fprintf('6. Config round-trip\n');
check_config(proj);
fprintf('7. Panel drawing, update, and preview\n');
png = check_panels();
fprintf('Preview PNG: %s\n', png);
fprintf('All streamline checks passed.\n');
end


function check_overrides()
[spacing_um, length_um] = flowfield_stream_overrides(struct());
assert(isempty(spacing_um) && isempty(length_um));
[spacing_um, length_um] = flowfield_stream_overrides(struct('stream_spacing_um', []));
assert(isempty(spacing_um) && isempty(length_um));
[spacing_um, length_um] = flowfield_stream_overrides(struct( ...
    'stream_spacing_um', 1.5, 'stream_length_um', 3));
assert(spacing_um == 1.5 && length_um == 3);
assert_error(@() flowfield_stream_overrides(struct('stream_spacing_um', 0)), ...
    'flowfield_stream_overrides:BadValue');
assert_error(@() flowfield_stream_overrides(struct('stream_spacing_um', -2)), ...
    'flowfield_stream_overrides:BadValue');
assert_error(@() flowfield_stream_overrides(struct( ...
    'stream_spacing_um', 2, 'stream_length_um', 0.5)), ...
    'flowfield_stream_overrides:BadLength');
end


function check_uniform()
h = -5:0.1:5;
v = h;
[H, V] = meshgrid(h, v);
Uh = ones(size(H));
Uv = zeros(size(V));
mask = false(size(H));
spacing_um = 0.9;
length_um = 4;
curves = flowfield_plane_streamlines(h, v, Uh, Uv, mask, [0; 0], ...
    spacing_um, length_um, 'xz');
n_expected = floor((h(end) - h(1)) / spacing_um);
assert(curves.n_lines >= n_expected);
assert(isequal(size(curves.faces), [curves.n_lines, 3]));
n_full = 0;
for i = 1:curves.n_lines
    poly = curves.polylines{i};
    arc = sum(hypot(diff(poly(:, 1)), diff(poly(:, 2))));
    assert(max(poly(:, 2)) - min(poly(:, 2)) < 1e-8);
    assert(all(diff(poly(:, 1)) > -1e-8));
    assert(arc + 1e-9 >= curves.stop_um);
    assert(arc <= 2 * length_um + 1e-6);
    tip = curves.vertices(curves.faces(i, 1), :);
    base = 0.5 * (curves.vertices(curves.faces(i, 2), :) ...
        + curves.vertices(curves.faces(i, 3), :));
    direction = tip - base;
    assert(direction(1) > 0 && abs(direction(2)) <= 1e-8 * direction(1));
    head = hypot(direction(1), direction(2));
    if abs(arc - 2 * length_um) <= 1e-6
        n_full = n_full + 1;
        assert(abs(head - 0.35 * spacing_um) < 1e-6);
    end
end
assert(n_full == n_expected, ...
    'Expected %d full lines, got %d.', n_expected, n_full);
sep = separation_across(curves.polylines, 0.5 * spacing_um);
assert(sep + 1e-6 >= spacing_um, ...
    'Parallel lines are %.4f um apart; spacing is %.4f um.', sep, spacing_um);
end


function check_vortex()
h = -4:0.05:4;
v = h;
[H, V] = meshgrid(h, v);
curves = flowfield_plane_streamlines(h, v, -V, H, false(size(H)), ...
    [0; 0], 1.0, 12, 'xy');
n_loop = 0;
for i = 1:curves.n_lines
    poly = curves.polylines{i};
    radius = mean(hypot(poly(:, 1), poly(:, 2)));
    if radius >= 1.2
        continue
    end
    circ = 2 * pi * radius;
    arc = sum(hypot(diff(poly(:, 1)), diff(poly(:, 2))));
    gap = hypot(poly(1, 1) - poly(end, 1), poly(1, 2) - poly(end, 2));
    assert(arc > 0.65 * circ && arc < 1.5 * circ, ...
        'Loop radius %.3f um, arc %.3f um, circumference %.3f um.', ...
        radius, arc, circ);
    assert(gap < 1.2, 'Loop gap %.3f um at radius %.3f um.', gap, radius);
    n_loop = n_loop + 1;
end
assert(n_loop >= 1, 'No closed orbit was kept.');
end


function check_mask()
h = -5:0.1:5;
v = h;
[H, V] = meshgrid(h, v);
mask = hypot(H, V) <= 1;
curves = flowfield_plane_streamlines(h, v, ones(size(H)), zeros(size(H)), ...
    mask, [0; 0], 0.8, 3, 'xz');
assert(curves.n_lines >= 4, 'Only %d lines around the body.', curves.n_lines);
near_body = false;
clear_line = false;
dh = h(2) - h(1);
for i = 1:curves.n_lines
    poly = curves.polylines{i};
    ih = min(max(round((poly(:, 1) - h(1)) / dh) + 1, 1), numel(h));
    iv = min(max(round((poly(:, 2) - v(1)) / dh) + 1, 1), numel(v));
    assert(~any(mask(sub2ind(size(mask), iv, ih))), ...
        'Line %d enters the body.', i);
    if any(hypot(poly(:, 1), poly(:, 2)) < 1.5)
        near_body = true;
    end
    if min(abs(poly(:, 2))) > 1.2
        arc = sum(hypot(diff(poly(:, 1)), diff(poly(:, 2))));
        if arc > 5
            assert(abs(arc - 6) <= curves.step_um + 1e-6, ...
                'Clear line arc %.4f um.', arc);
            clear_line = true;
        end
    end
    tip = curves.vertices(curves.faces(i, 1), :);
    base = 0.5 * (curves.vertices(curves.faces(i, 2), :) ...
        + curves.vertices(curves.faces(i, 3), :));
    assert(tip(1) > base(1));
end
assert(near_body, 'No curve reaches the body.');
assert(clear_line, 'No curve clears the body.');
sep = separation_across(curves.polylines, 0.5 * 0.8);
assert(sep + 1e-6 >= 0.8, ...
    'Parallel lines are %.4f um apart.', sep);
end


function check_rejects()
h = -5:0.1:5;
v = h;
[H, ~] = meshgrid(h, v);
Uh = ones(size(H));
Uv = zeros(size(H));
mask = false(size(H));
assert_error(@() flowfield_plane_streamlines(h, v, Uh, Uv, mask, [0; 0], ...
    0, 4, 'xz'), 'flowfield_plane_streamlines:BadSpacing');
assert_error(@() flowfield_plane_streamlines(h, v, Uh, Uv, mask, [0; 0], ...
    20, 4, 'xz'), 'flowfield_plane_streamlines:BadSpacing');
assert_error(@() flowfield_plane_streamlines(h, v, Uh, Uv, mask, [0; 0], ...
    0.01, 4, 'xz'), 'flowfield_plane_streamlines:BadSpacing');
assert_error(@() flowfield_plane_streamlines(h, v, Uh, Uv, mask, [0; 0], ...
    1, 0.2, 'xz'), 'flowfield_plane_streamlines:BadLength');
bad_h = [0, 1, 2, 2.2, 4];
bad_v = 0:4;
assert_error(@() flowfield_plane_streamlines(bad_h, bad_v, ...
    ones(5, 5), zeros(5, 5), false(5, 5), [0; 0], 1, 2, 'xz'), ...
    'flowfield_plane_streamlines:BadGrid');
Uh(3, 3) = Inf;
assert_error(@() flowfield_plane_streamlines(h, v, Uh, Uv, mask, [0; 0], ...
    1, 2, 'xz'), 'flowfield_plane_streamlines:BadVelocity');
end


function check_config(proj)
example = fullfile(proj, 'configs', 'config_example.json');
cfg = load_sim_config(example);
assert(isscalar(cfg.flowfield_stream_spacing) && isnan(cfg.flowfield_stream_spacing));
assert(isscalar(cfg.flowfield_stream_length) && isnan(cfg.flowfield_stream_length));

cfg.flowfield_stream_spacing = 1.25;
cfg.flowfield_stream_length = 4;
tmp = [tempname, '.json'];
cleaner = onCleanup(@() delete_if(tmp));
save_sim_config(cfg, tmp);
loaded = load_sim_config(tmp);
assert(abs(loaded.flowfield_stream_spacing - 1.25) < 1e-12);
assert(abs(loaded.flowfield_stream_length - 4) < 1e-12);

cfg.flowfield_stream_spacing = NaN;
cfg.flowfield_stream_length = NaN;
save_sim_config(cfg, tmp);
loaded = load_sim_config(tmp);
assert(isnan(loaded.flowfield_stream_spacing));
assert(isnan(loaded.flowfield_stream_length));
clear cleaner

raw = jsondecode(fileread(example));
raw.flowfield.stream_length_um = 0;
bad = [tempname, '.json'];
bad_cleaner = onCleanup(@() delete_if(bad));
fid = fopen(bad, 'w');
assert(fid > 0);
fwrite(fid, jsonencode(raw));
fclose(fid);
assert_error(@() load_sim_config(bad), 'load_sim_config:BadFlowfield');
clear bad_cleaner
end


function png = check_panels()
frame = synthetic_frame();
handles = render_flowfield_panels(frame, struct('action', 'create', 'visible', 'off'));
closer = onCleanup(@() close_if(handles.fig));
tags = {'d', 'e', 'f'};
names = {'xz', 'yz', 'xy'};
for n = 1:3
    expected = expected_geometry(frame, names{n});
    dark = handles.(['line_k_' tags{n}]);
    head = handles.(['head_' tags{n}]);
    assert(isequaln(dark.XData(:), expected.x(:)));
    assert(isequaln(dark.YData(:), expected.y(:)));
    assert(max(abs(dark.Color - [0.12, 0.12, 0.12])) < 1e-12);
    assert(dark.LineWidth == 0.75);
    assert(~isfield(handles, ['line_w_' tags{n}]));
    assert(isequaln(head.Vertices, expected.vertices));
    assert(isempty(findobj(handles.(['ax_' tags{n}]), 'Type', 'quiver')));
end
assert_stack(handles);
assert(~isempty(findobj(handles.ax_a, 'Type', 'surface')));

n_child = numel(handles.ax_d.Children);
x_before = handles.line_k_d.XData;
frame.planes.xz.Ux = -frame.planes.xz.Ux;
handles = render_flowfield_panels(frame, struct( ...
    'action', 'update', 'handles', handles));
assert(numel(handles.ax_d.Children) == n_child);
expected = expected_geometry(frame, 'xz');
assert(isequaln(handles.line_k_d.XData(:), expected.x(:)));
assert(~isequaln(handles.line_k_d.XData, x_before));
assert_arrows_sign(handles.head_d, -1);
assert_stack(handles);

png = fullfile(tempdir, 'flowfield_streamline_acceptance.png');
render_flowfield_panels([], struct('action', 'export_png', ...
    'figure', handles.fig, 'png_path', png, 'dpi', 40));
clear closer
end


function assert_stack(handles)
%ASSERT_STACK  Body outline and markers sit above the streamlines
ch = handles.ax_d.Children;
stack_index = @(h) find(ch == h, 1);
order = [stack_index(handles.tag_d), stack_index(handles.cross_d), stack_index(handles.outline_d), ...
    stack_index(handles.head_d), stack_index(handles.line_k_d), stack_index(handles.img_d)];
assert(all(isfinite(order)) && issorted(order, 'strictascend'), ...
    '2D object order is label, crossings, body, arrow, curve, image.');
end


function assert_arrows_sign(head, sign_x)
verts = head.Vertices;
faces = head.Faces;
for i = 1:size(faces, 1)
    tri = verts(faces(i, :), :);
    direction = tri(1, :) - 0.5 * (tri(2, :) + tri(3, :));
    assert(sign_x * direction(1) > 0);
end
end


function curves = expected_geometry(frame, plane_name)
%EXPECTED_GEOMETRY  Same component map as plane_stream_geometry
plane = frame.planes.(plane_name);
switch plane_name
    case 'xz'
        Uh = plane.Ux;
        Uv = plane.Uz;
        origin = [0; 0];
    case 'yz'
        Uh = plane.Uy;
        Uv = plane.Uz;
        origin = [0; 0];
    case 'xy'
        Uh = plane.Ux;
        Uv = plane.Uy;
        origin = [0; 0];
end
curves = flowfield_plane_streamlines(plane.h_vec, plane.v_vec, Uh, Uv, ...
    plane.mask, origin, frame.planes.stream_spacing_um, ...
    frame.planes.stream_length_um, plane_name);
end


function frame = synthetic_frame()
h = -5:0.25:5;
v = h;
[H, V] = meshgrid(h, v);
mask = hypot(H, V) <= 1;
th = linspace(0, 2 * pi, 64);
cap_h = cos(th);
cap_v = sin(th);
z0 = zeros(size(H));
frame.planes.xz = pack_plane(h, v, H, z0, V, ones(size(H)), z0, z0, ...
    mask, cap_h, cap_v);
frame.planes.yz = pack_plane(h, v, z0, H, V, z0, z0, ones(size(H)), ...
    mask, cap_h, cap_v);
frame.planes.xy = pack_plane(h, v, H, V, z0, -V, H, z0, ...
    mask, cap_h, cap_v);
frame.planes.body_vis = struct('type', 'capsule', 'a_um', 1.0, 'b_um', 0.35);
frame.planes.axis_lim_um = [-5, 5; -5, 5; -5, 5];
frame.planes.grid_spacing_um = 0.25;
frame.planes.quiver_skip = 5;
frame.planes.stream_spacing_um = 1.2;
frame.planes.stream_length_um = 2.5;
frame.planes.X_cm_um = [0; 0; 0];
frame.geom = struct('R_body', eye(3), 'X_pili', {{}});
frame.params = struct('M_pili', 0);
frame.L_to_um = 1;
frame.cmax = 1;
frame.time_s = 1.25;
frame.run_id = 'acceptance';
frame.scale_note = '';
end


function plane = pack_plane(h, v, HG, VG, ZG, Ux, Uy, Uz, mask, cap_h, cap_v)
plane.h_vec = h;
plane.v_vec = v;
plane.HG = HG;
plane.VG = VG;
plane.ZG = ZG;
plane.Ux = Ux;
plane.Uy = Uy;
plane.Uz = Uz;
plane.Umag = ones(size(mask));
plane.Umag(mask) = NaN;
plane.mask = mask;
plane.cap_h = cap_h;
plane.cap_v = cap_v;
end


function sep = separation_across(parts, min_shift)
%SEPARATION_ACROSS  Minimum distance between curves on different rows
sep = inf;
for a = 1:numel(parts) - 1
    ya = mean(parts{a}(:, 2));
    A = parts{a};
    for b = a + 1:numel(parts)
        if abs(mean(parts{b}(:, 2)) - ya) < min_shift
            continue
        end
        B = parts{b};
        d2 = (A(:, 1) - B(:, 1).').^2 + (A(:, 2) - B(:, 2).').^2;
        sep = min(sep, sqrt(min(d2(:))));
    end
end
end


function assert_error(fn, identifier)
threw = false;
ME = [];
try
    fn();
catch ME
    threw = true;
end
if ~threw
    error('check_flowfield_streamlines:MissedError', ...
        'Expected error %s.', identifier);
end
if ~strcmp(ME.identifier, identifier)
    error('check_flowfield_streamlines:WrongError', ...
        'Expected %s, got %s (%s).', identifier, ME.identifier, ME.message);
end
end


function close_if(fig)
if ~isempty(fig) && isgraphics(fig)
    close(fig);
end
end


function delete_if(path)
if isfile(path)
    delete(path);
end
end
