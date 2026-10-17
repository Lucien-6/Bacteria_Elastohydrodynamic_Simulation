function info = render_motion_animation(results, params, output_dir, opts)
%RENDER_MOTION_ANIMATION  Render the 3-panel 3-D motion animation (3D + XZ + YZ)
%
%  Description:
%    Renders the bacterium motion animation from a finished simulation and
%    writes it to a video file in output_dir.  The layout is a 1 x 3 figure:
%      left   : 3-D view (azimuth 45 deg, elevation 25 deg) with gravity arrow
%      centre : XZ projection (view along Y)
%      right  : YZ projection (view along X)
%    Each panel shows the semi-transparent body surface, every pilus, the
%    centre-of-mass marker and the full-resolution centre-of-mass trajectory
%    up to the current frame.
%
%    Rendering strategy: every graphics object (body, pili, markers,
%    trajectory, border lines, lights) is created ONCE; each frame only
%    updates the XData/YData/ZData (or Vertices) of the existing objects and
%    the 3-D title string before drawnow / getframe / writeVideo.  This
%    replaces the former per-frame cla + re-creation of all objects (and of
%    lights, labels and limits), which dominated the render time.
%
%    The function is shared by run_simulation.m, run_simulation_HPC_Linux.m
%    and post_processing/render_animation_from_results.m.
%
%    Video profile (platform aware): 'MPEG-4' (.mp4) on Windows / macOS;
%    'Motion JPEG AVI' (.avi) on Linux, where MATLAB does not provide the
%    MPEG-4 profile.  Quality is 100 for both.
%
%  Syntax:
%    info = render_motion_animation(results, params, output_dir, opts)
%
%  Inputs:
%    results    - Results structure of run_simulation (fields used: t, q,
%                 X_cm_um, X_pili_um); all arrays must cover the same NT time
%                 points.
%    params     - Parameter structure (fields used: M_pili, geometry_mode,
%                 body_type, dim.major_axis, dim.minor_axis [m], and for
%                 custom bodies body_verts [m] and body_faces).
%    output_dir - Output directory (must exist).
%    opts       - Options structure:
%                   .fps            Video frame rate [frames/s]       (required)
%                   .skip           Frame stride over the NT samples  (required)
%                   .silent         Suppress progress bar (default: false)
%                   .logf           Logging handle, console + file
%                                   (default: fprintf to the console)
%                   .logf_file      Logging handle, file only (default: no-op)
%                   .file_basename  Video base name without extension
%                                   (default: 'bacteria_3d_motion')
%
%  Outputs:
%    info - Structure with fields:
%             .file       Full path of the video file
%             .filename   File name (with extension)
%             .n_frames   Number of rendered frames
%             .elapsed_s  Total wall time of this function [s]
%             .stage_s    Structure of per-stage wall times [s]:
%                           prepass, setup, update, drawnow, getframe,
%                           resize, write
%
%  Notes:
%    The capsule mesh resolution and the per-pilus colour palette below
%    intentionally mirror the STATIC VISUALISATION section of run_simulation.m
%    and run_simulation_HPC_Linux.m (snapshots of the summary figure); keep them
%    in sync when either changes.  The body and trajectory colours are defined
%    only here.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-17
%  Updated: 2026-10-17

%% ---- Options ----
if nargin < 4 || ~isstruct(opts)
    error('render_motion_animation:InvalidOptions', ...
          'opts must be a struct with fields fps and skip.');
end
if ~isfield(opts, 'fps') || ~isfield(opts, 'skip')
    error('render_motion_animation:InvalidOptions', ...
          'opts.fps and opts.skip are required.');
end
animation_fps  = opts.fps;
animation_skip = opts.skip;
silent         = logical(get_opt(opts, 'silent', false));
file_basename  = get_opt(opts, 'file_basename', 'bacteria_3d_motion');
logf           = get_opt(opts, 'logf', @(varargin) fprintf(varargin{:}));
logf_file      = get_opt(opts, 'logf_file', @(varargin) []);

t_anim_start = tic;
logf('\n[ANIM] Generating 3D motion animation...\n');

%% ---- Data ----
NT      = numel(results.t);
M_pili  = params.M_pili;
if size(results.q, 2) ~= NT || size(results.X_cm_um, 2) ~= NT
    error('render_motion_animation:InconsistentResults', ...
          'results.t, results.q and results.X_cm_um must cover the same time points.');
end

is_capsule = strcmp(params.geometry_mode, 'parametric') && strcmp(params.body_type, 'capsule');
if is_capsule
    a_um = params.dim.major_axis * 1e6;   % Semi-major axis [um]
    b_um = params.dim.minor_axis * 1e6;   % Semi-minor axis (cap radius) [um]
    [capsule_flat, mesh_size] = build_capsule_mesh(a_um, b_um);
    custom_verts_um = [];
else
    custom_verts_um = params.body_verts * 1e6;   % [V x 3] in um
    capsule_flat = [];
    mesh_size    = [];
    a_um         = [];
end

% --- Colour palette (per-pilus palette mirrors the static visualisation section) ---
body_color    = [0.20 0.65 0.33];        % Forest green
body_alpha    = 0.35;                    % Body transparency
traj_color    = [0.85 0.25 0.25 0.45];   % Deeper red, 45 % opacity
pili_palette = [0.20 0.40 0.80;
                0.85 0.33 0.10;
                0.47 0.67 0.19;
                0.64 0.08 0.18;
                0.49 0.18 0.56;
                0.93 0.69 0.13;
                0.30 0.75 0.93;
                0.60 0.60 0.60;
                0.00 0.45 0.74;
                0.85 0.47 0.66];
pili_colors = pili_palette(mod((1:M_pili)-1, size(pili_palette,1)) + 1, :);

stage_s = struct('prepass', 0, 'setup', 0, 'update', 0, 'drawnow', 0, ...
                 'getframe', 0, 'resize', 0, 'write', 0);

%% ---- Frame indices ----
frame_indices = 1 : animation_skip : NT;
if frame_indices(end) ~= NT
    frame_indices(end+1) = NT;      % Always include last frame
end
n_frames = length(frame_indices);
logf('[ANIM] Frames: %d (skip = %d from %d total time points)\n', ...
    n_frames, animation_skip, NT);

% Pre-extract CM trajectory for all animation frames [3 x n_frames]
cm_traj_um = results.X_cm_um(:, frame_indices);

% Full-resolution CM trajectory for smooth rendering (all NT points)
cm_traj_full_um = results.X_cm_um;

%% ---- Global axis limits (body + pili + CM over ALL animation frames) ----
t_stage = tic;
x_min = min(cm_traj_um(1,:));  x_max = max(cm_traj_um(1,:));
y_min = min(cm_traj_um(2,:));  y_max = max(cm_traj_um(2,:));
z_min = min(cm_traj_um(3,:));  z_max = max(cm_traj_um(3,:));
for idx = 1:n_frames
    ti   = frame_indices(idx);
    R_fi = quat_to_rotmat(results.q(:, ti));
    if is_capsule
        pts = R_fi * capsule_flat + results.X_cm_um(:, ti);          % [3 x N_mesh]
    else
        pts = R_fi * custom_verts_um' + results.X_cm_um(:, ti);      % [3 x N_verts]
    end
    x_min = min(x_min, min(pts(1,:)));  x_max = max(x_max, max(pts(1,:)));
    y_min = min(y_min, min(pts(2,:)));  y_max = max(y_max, max(pts(2,:)));
    z_min = min(z_min, min(pts(3,:)));  z_max = max(z_max, max(pts(3,:)));
end
for k = 1:M_pili
    Xp_all = results.X_pili_um{k}(:, :, frame_indices);               % [Nn x 3 x n_frames]
    mn = min(Xp_all, [], [1 3]);
    mx = max(Xp_all, [], [1 3]);
    x_min = min(x_min, mn(1));  x_max = max(x_max, mx(1));
    y_min = min(y_min, mn(2));  y_max = max(y_max, mx(2));
    z_min = min(z_min, mn(3));  z_max = max(z_max, mx(3));
end
clear Xp_all

% Build per-axis fixed limits with independent margins
x_span = x_max - x_min;
y_span = y_max - y_min;
z_span = z_max - z_min;
margin_x = max(0.10 * x_span, 0.15);
margin_y = max(0.10 * y_span, 0.15);
margin_z = max(0.10 * z_span, 0.15);
ax_lim_x = [x_min - margin_x, x_max + margin_x];
ax_lim_y = [y_min - margin_y, y_max + margin_y];
ax_lim_z = [z_min - margin_z, z_max + margin_z];

% Enforce minimum axis range
if is_capsule
    min_ax_range = 4 * 2 * a_um;
else
    body_extents = max(custom_verts_um) - min(custom_verts_um);
    min_ax_range = 4 * max(body_extents);
end
if diff(ax_lim_x) < min_ax_range
    cx = mean(ax_lim_x); ax_lim_x = cx + [-1, 1] * min_ax_range / 2;
end
if diff(ax_lim_y) < min_ax_range
    cy = mean(ax_lim_y); ax_lim_y = cy + [-1, 1] * min_ax_range / 2;
end
if diff(ax_lim_z) < min_ax_range
    cz = mean(ax_lim_z); ax_lim_z = cz + [-1, 1] * min_ax_range / 2;
end

logf('[ANIM] Global axis limits (per-axis): X=[%.2f, %.2f], Y=[%.2f, %.2f], Z=[%.2f, %.2f] um\n', ...
    ax_lim_x(1), ax_lim_x(2), ax_lim_y(1), ax_lim_y(2), ax_lim_z(1), ax_lim_z(2));
stage_s.prepass = toc(t_stage);

%% ---- Figure dimensions (minimise left/right whitespace) ----
% With DataAspectRatio=[1 1 1], MATLAB centres the equal-aspect data box
% inside the axes Position rectangle.  The 3-D bounding box is projected onto
% the screen at view(45,25) to derive the natural content aspect ratio, and
% the figure is sized accordingly so the data box fills it without padding.
% The axes inner area occupies approximately 70 % (width) x 78 % (height) of
% the figure in a single-axes 3-D plot with title + axis labels.
t_stage = tic;
az_v = 45;  el_v = 25;
proj_w_data = diff(ax_lim_x) * abs(cosd(az_v)) + diff(ax_lim_y) * abs(sind(az_v));
proj_h_data = (diff(ax_lim_x) * abs(sind(az_v)) + diff(ax_lim_y) * abs(cosd(az_v))) ...
              * abs(sind(el_v)) + diff(ax_lim_z) * abs(cosd(el_v));
axes_fill_w = 0.70;   % approx fraction of figure width used by axes
axes_fill_h = 0.78;   % approx fraction of figure height used by axes
anim_fig_h_ref = 1440;
fig_w_single   = round(anim_fig_h_ref * (proj_w_data / proj_h_data) ...
                       * (axes_fill_h / axes_fill_w));
% 1x3 panel layout: triple the single-panel width; enforce minimum and H.264 even
anim_fig_w = 2 * ceil(max(round(fig_w_single * 2), 1800) / 2);
anim_fig_h = 1200;

% --- Fit within codec-safe bounding box while preserving aspect ratio ----
anim_max_w = 3840;   % safe upper bound for H.264 width  (even)
anim_max_h = 2160;   % safe upper bound for H.264 height (even)
ideal_ar   = anim_fig_w / anim_fig_h;   % target aspect ratio
if anim_fig_w > anim_max_w || anim_fig_h > anim_max_h
    scale = min(anim_max_w / anim_fig_w, anim_max_h / anim_fig_h);
    anim_fig_w_new = 2 * ceil(max(round(anim_fig_w * scale), 1800) / 2);
    anim_fig_h_new = 2 * round(anim_fig_w_new / ideal_ar / 2);
    anim_fig_h_new = max(anim_fig_h_new, 480);
    logf('[ANIM] Ideal size %d x %d exceeds codec limit [%d x %d] -rescaled to %d x %d.\n', ...
        anim_fig_w, anim_fig_h, anim_max_w, anim_max_h, anim_fig_w_new, anim_fig_h_new);
    anim_fig_w = anim_fig_w_new;
    anim_fig_h = anim_fig_h_new;
end
logf('[ANIM] Figure size (logical): %d x %d px  (1x3 layout: 3D + XZ + YZ)\n', ...
    anim_fig_w, anim_fig_h);

%% ---- Figure and axes (created once, invisible) ----
fig_anim = figure('Units', 'pixels', 'Position', [50, 50, anim_fig_w, anim_fig_h], ...
                  'Color', 'w', 'Visible', 'off');
cleanup_fig = onCleanup(@() safe_close_figure(fig_anim));
set(fig_anim, 'DefaultAxesFontName',  'Arial', ...
              'DefaultAxesFontSize',  11, ...
              'DefaultAxesLineWidth', 0.8);

% Three axes in a single row: 3D (left), XZ (centre), YZ (right)
ax_3d = axes(fig_anim, 'Position', [0.03, 0.08, 0.30, 0.84]);
ax_xz = axes(fig_anim, 'Position', [0.36, 0.08, 0.30, 0.84]);
ax_yz = axes(fig_anim, 'Position', [0.69, 0.08, 0.30, 0.84]);
ax_all = [ax_3d, ax_xz, ax_yz];

% Pre-configure 3D axes so Position is settled before gravity indicator
configure_axes_view(ax_3d, ax_lim_x, ax_lim_y, ax_lim_z, 45, 25, 'on', false);
xlabel(ax_3d, 'X (μm)', 'FontName', 'Arial', 'FontSize', 12);
ylabel(ax_3d, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 12);
zlabel(ax_3d, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 12);
h_title_3d = title(ax_3d, sprintf('Time = %.4f s  (frame %d/%d)', 0.0, 1, n_frames), ...
    'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');

% Pre-configure XZ projection axes (view along Y axis)
configure_axes_view(ax_xz, ax_lim_x, ax_lim_y, ax_lim_z, 0, 0, 'off', true);
xlabel(ax_xz, 'X (μm)', 'FontName', 'Arial', 'FontSize', 11);
zlabel(ax_xz, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
title(ax_xz, 'XZ Projection', ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

% Pre-configure YZ projection axes (view along X axis)
configure_axes_view(ax_yz, ax_lim_x, ax_lim_y, ax_lim_z, 90, 0, 'off', true);
ylabel(ax_yz, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 11);
zlabel(ax_yz, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
title(ax_yz, 'YZ Projection', ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

drawnow;
add_gravity_indicator_outside(fig_anim, ax_3d, 13, 2.5);

%% ---- Probe actual physical pixel size and enforce codec limits ----
% On HiDPI / display scaling (e.g. 125%), getframe returns physical pixels
% larger than the logical Position (e.g. 2250x1366 vs 1800x1200).  The codec
% limit must be checked against the PHYSICAL size, not the logical one.
probe_frame = getframe(fig_anim);
[phys_h, phys_w, ~] = size(probe_frame.cdata);
logf('[ANIM] Probe frame: %d x %d physical px (logical %d x %d)\n', ...
    phys_w, phys_h, anim_fig_w, anim_fig_h);

if phys_w > anim_max_w || phys_h > anim_max_h
    shrink = min(anim_max_w / phys_w, anim_max_h / phys_h);
    new_log_w = 2 * floor(anim_fig_w * shrink / 2);
    new_log_h = 2 * floor(anim_fig_h * shrink / 2);
    new_log_w = max(new_log_w, 400);
    new_log_h = max(new_log_h, 300);
    logf('[ANIM] Physical size exceeds codec limit [%d x %d] -shrinking figure to %d x %d logical.\n', ...
        anim_max_w, anim_max_h, new_log_w, new_log_h);
    set(fig_anim, 'Position', [50, 50, new_log_w, new_log_h]);
    drawnow;
    probe_frame = getframe(fig_anim);
    [phys_h, phys_w, ~] = size(probe_frame.cdata);
    logf('[ANIM] Re-probe: %d x %d physical px\n', phys_w, phys_h);
end

% Hi-DPI scaling can keep the captured physical frame larger than the logical
% figure even after resizing. Clamp the encoded frame size to a codec-safe,
% even-valued box and resize captured frames on the fly if needed.
[vid_h, vid_w] = fit_anim_frame_size(phys_h, phys_w, anim_max_h, anim_max_w);
if vid_w ~= phys_w || vid_h ~= phys_h
    logf(['[ANIM] Captured frame still exceeds video codec constraints; ', ...
          'encoding will use %d x %d px and resize frames as needed.\n'], ...
          vid_w, vid_h);
end

%% ---- Open VideoWriter (platform-aware profile) ----
% 'MPEG-4' is documented by MathWorks as Windows/Mac-only and throws "The
% specified profile is not valid." on any Linux system.  'Motion JPEG AVI' is
% MATLAB's built-in, dependency-free profile that works on every platform, so
% it is used on Linux (e.g. headless HPC nodes).
if ispc || ismac
    anim_filename = [file_basename '.mp4'];
    anim_profile  = 'MPEG-4';
else
    anim_filename = [file_basename '.avi'];
    anim_profile  = 'Motion JPEG AVI';
end
anim_file = fullfile(output_dir, anim_filename);
vid = VideoWriter(anim_file, anim_profile);
cleanup_vid = onCleanup(@() safe_close_video(vid));
vid.FrameRate = animation_fps;
vid.Quality   = 100;
open(vid);

%% ---- Create all graphics objects ONCE (child order identical per axis) ----
% Per axis: body, pili, centre-of-mass marker, trajectory, [XZ/YZ border
% lines], then lights.  Only their data are updated inside the frame loop.
h_body = gobjects(1, 3);
h_pili = gobjects(M_pili, 3);
h_mark = gobjects(1, 3);
h_traj = gobjects(1, 3);

R_1  = quat_to_rotmat(results.q(:, 1));
cm_1 = results.X_cm_um(:, 1);

for ai = 1:3
    ax = ax_all(ai);
    hold(ax, 'on');

    % Body surface (dispatch by mesh type)
    if is_capsule
        pts_1 = R_1 * capsule_flat + cm_1;
        h_body(ai) = surf(ax, reshape(pts_1(1,:), mesh_size), ...
                              reshape(pts_1(2,:), mesh_size), ...
                              reshape(pts_1(3,:), mesh_size), ...
                          'FaceColor', body_color, 'FaceAlpha', body_alpha, ...
                          'EdgeColor', 'none', 'FaceLighting', 'gouraud');
    else
        verts_1 = (R_1 * custom_verts_um' + cm_1)';   % [V x 3]
        % 'Parent' is used instead of trisurf(ax, ...), see the STATIC
        % VISUALISATION section of run_simulation.m.
        h_body(ai) = trisurf(params.body_faces, ...
                             verts_1(:,1), verts_1(:,2), verts_1(:,3), ...
                             'FaceColor', body_color, 'FaceAlpha', body_alpha, ...
                             'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'Parent', ax);
    end

    for k = 1:M_pili
        Xp = results.X_pili_um{k}(:, :, 1);
        h_pili(k, ai) = plot3(ax, Xp(:,1), Xp(:,2), Xp(:,3), ...
                              '-', 'Color', pili_colors(k, :), 'LineWidth', 1.25);
    end

    h_mark(ai) = plot3(ax, cm_1(1), cm_1(2), cm_1(3), ...
                       'r+', 'MarkerSize', 5, 'LineWidth', 1);

    % Trajectory: nothing visible on the first frame (NaN point; plot3 with
    % empty data returns no handle), grows with the frame index afterwards
    h_traj(ai) = plot3(ax, NaN, NaN, NaN, ...
                       '-', 'Color', traj_color, 'LineWidth', 1.0);

    hold(ax, 'off');
end

% --- 3D panel: lights, lighting, material, final view ---
camlight(ax_3d, 'headlight');
camlight(ax_3d, 'left');
lighting(ax_3d, 'gouraud');
material(ax_3d, 'dull');
configure_axes_view(ax_3d, ax_lim_x, ax_lim_y, ax_lim_z, 45, 25, 'on', false);

% --- XZ panel: light, material, view, border lines ---
camlight(ax_xz, 'headlight');
lighting(ax_xz, 'gouraud');
material(ax_xz, 'dull');
configure_axes_view(ax_xz, ax_lim_x, ax_lim_y, ax_lim_z, 0, 0, 'off', true);
hold(ax_xz, 'on');
plot3(ax_xz, [ax_lim_x(2) ax_lim_x(2)], ax_lim_y([1 1]), ax_lim_z, ...
    'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot3(ax_xz, ax_lim_x, ax_lim_y([1 1]), [ax_lim_z(2) ax_lim_z(2)], ...
    'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold(ax_xz, 'off');

% --- YZ panel: light, material, view, border lines ---
camlight(ax_yz, 'headlight');
lighting(ax_yz, 'gouraud');
material(ax_yz, 'dull');
configure_axes_view(ax_yz, ax_lim_x, ax_lim_y, ax_lim_z, 90, 0, 'off', true);
hold(ax_yz, 'on');
plot3(ax_yz, ax_lim_x([1 1]), [ax_lim_y(2) ax_lim_y(2)], ax_lim_z, ...
    'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot3(ax_yz, ax_lim_x([1 1]), ax_lim_y, [ax_lim_z(2) ax_lim_z(2)], ...
    'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold(ax_yz, 'off');

stage_s.setup = toc(t_stage);

%% ---- Frame loop: update data only ----
for fi = 1:n_frames
    ti = frame_indices(fi);

    % --- Update geometry of the existing objects ---
    t_stage = tic;
    R_fi  = quat_to_rotmat(results.q(:, ti));
    cm_um = results.X_cm_um(:, ti);

    if is_capsule
        pts = R_fi * capsule_flat + cm_um;                       % [3 x N_mesh]
        Xs = reshape(pts(1,:), mesh_size);
        Ys = reshape(pts(2,:), mesh_size);
        Zs = reshape(pts(3,:), mesh_size);
    else
        verts = (R_fi * custom_verts_um' + cm_um)';              % [V x 3]
    end

    for ai = 1:3
        if is_capsule
            set(h_body(ai), 'XData', Xs, 'YData', Ys, 'ZData', Zs);
        else
            set(h_body(ai), 'Vertices', verts);
        end
    end

    for k = 1:M_pili
        Xp = results.X_pili_um{k}(:, :, ti);
        for ai = 1:3
            set(h_pili(k, ai), 'XData', Xp(:,1), 'YData', Xp(:,2), 'ZData', Xp(:,3));
        end
    end

    for ai = 1:3
        set(h_mark(ai), 'XData', cm_um(1), 'YData', cm_um(2), 'ZData', cm_um(3));
    end

    if fi > 1
        for ai = 1:3
            set(h_traj(ai), 'XData', cm_traj_full_um(1, 1:ti), ...
                            'YData', cm_traj_full_um(2, 1:ti), ...
                            'ZData', cm_traj_full_um(3, 1:ti));
        end
    end

    set(h_title_3d, 'String', sprintf('Time = %.4f s  (frame %d/%d)', ...
        results.t(ti), fi, n_frames));
    stage_s.update = stage_s.update + toc(t_stage);

    % --- Flush rendering ---
    t_stage = tic;
    drawnow;
    stage_s.drawnow = stage_s.drawnow + toc(t_stage);

    % --- Capture frame (no rect - natural physical pixels) ---
    t_stage = tic;
    frame = getframe(fig_anim);
    stage_s.getframe = stage_s.getframe + toc(t_stage);

    t_stage = tic;
    [fh, fw, ~] = size(frame.cdata);
    if fh ~= vid_h || fw ~= vid_w
        frame.cdata = normalize_anim_frame_cdata(frame.cdata, vid_h, vid_w);
    end
    stage_s.resize = stage_s.resize + toc(t_stage);

    t_stage = tic;
    writeVideo(vid, frame);
    stage_s.write = stage_s.write + toc(t_stage);

    % In-place progress bar (console only) - suppressed in silent mode
    if ~silent
        render_progress_bar('[ANIM]   ', fi, n_frames, t_anim_start);
    end

    % Log animation progress to file every 10% (file-only, no console echo)
    if mod(fi, max(1, floor(n_frames / 10))) == 0 || fi == n_frames
        logf_file('[ANIM] Frame %d / %d rendered\n', fi, n_frames);
    end
end

close(vid);
close(fig_anim);

%% ---- Summary ----
t_anim_elapsed = toc(t_anim_start);
vid_info = dir(anim_file);
logf('[ANIM] Animation saved: %s (%.2f MB, %d frames, %.1f s render)\n', ...
    anim_file, vid_info.bytes/1e6, n_frames, t_anim_elapsed);

stage_names  = {'prepass', 'setup', 'update', 'drawnow', 'getframe', 'resize', 'write'};
stage_labels = {'axis-limit pre-scan', 'figure + objects setup', 'data update (set)', ...
                'drawnow', 'getframe', 'frame resize', 'writeVideo'};
logf('[ANIM] Stage timing (total s | share | ms per frame):\n');
for is = 1:numel(stage_names)
    ts = stage_s.(stage_names{is});
    logf('[ANIM]   %-24s : %8.2f s | %5.1f %% | %8.1f ms\n', ...
        stage_labels{is}, ts, 100 * ts / max(t_anim_elapsed, eps), 1e3 * ts / n_frames);
end

info = struct('file', anim_file, 'filename', anim_filename, ...
              'n_frames', n_frames, 'elapsed_s', t_anim_elapsed, ...
              'stage_s', stage_s);

end


%% ==================== Local Functions ====================

function v = get_opt(S, fname, default_val)
%GET_OPT  Return S.(fname) if present and non-empty, else default_val
    if isstruct(S) && isfield(S, fname) && ~isempty(S.(fname))
        v = S.(fname);
    else
        v = default_val;
    end
end

function configure_axes_view(ax, lim_x, lim_y, lim_z, az, el, box_state, ortho)
%CONFIGURE_AXES_VIEW  Apply fixed limits, equal aspect, grid, box, ticks and view
    xlim(ax, lim_x);  ylim(ax, lim_y);  zlim(ax, lim_z);
    if ortho
        set(ax, 'DataAspectRatio', [1 1 1], 'TickDir', 'out', ...
            'Projection', 'orthographic');
    else
        set(ax, 'DataAspectRatio', [1 1 1], 'TickDir', 'out');
    end
    grid(ax, 'on');  box(ax, box_state);
    view(ax, az, el);
end

function [capsule_flat, mesh_size] = build_capsule_mesh(a_um, b_um)
%BUILD_CAPSULE_MESH  Parametric capsule surface mesh [um], flattened to [3 x N_mesh]
%  Mirrors the capsule mesh of the STATIC VISUALISATION section of
%  run_simulation.m / run_simulation_HPC_Linux.m.
    n_circ = 36;                   % Circumferential resolution
    n_cap  = 12;                   % Meridional resolution per hemisphere cap
    n_cyl  = 10;                   % Axial resolution for cylinder section

    theta_cap = linspace(0, 2*pi, n_circ);

    % --- Left hemisphere (pole -> equator) ---
    phi_L = linspace(pi, pi/2, n_cap);
    [TH_L, PH_L] = meshgrid(theta_cap, phi_L);
    XL = b_um * cos(PH_L) - (a_um - b_um);
    YL = b_um * sin(PH_L) .* cos(TH_L);
    ZL = b_um * sin(PH_L) .* sin(TH_L);

    % --- Cylinder section ---
    x_cyl_nodes = linspace(-(a_um - b_um), (a_um - b_um), n_cyl);
    [TH_C, XC_grid] = meshgrid(theta_cap, x_cyl_nodes);
    XC = XC_grid;
    YC = b_um * cos(TH_C);
    ZC = b_um * sin(TH_C);

    % --- Right hemisphere (equator -> pole) ---
    phi_R = linspace(pi/2, 0, n_cap);
    [TH_R, PH_R] = meshgrid(theta_cap, phi_R);
    XR = b_um * cos(PH_R) + (a_um - b_um);
    YR = b_um * sin(PH_R) .* cos(TH_R);
    ZR = b_um * sin(PH_R) .* sin(TH_R);

    % --- Concatenate into single continuous surface ---
    capsule_X0 = [XL; XC; XR];
    capsule_Y0 = [YL; YC; YR];
    capsule_Z0 = [ZL; ZC; ZR];
    mesh_size  = size(capsule_X0);

    % Pre-flatten for fast rotation: [3 x N_mesh]
    capsule_flat = [capsule_X0(:)'; capsule_Y0(:)'; capsule_Z0(:)'];
end

function safe_close_figure(fig)
%SAFE_CLOSE_FIGURE  Close a figure if it is still valid (onCleanup target)
    if ~isempty(fig) && isgraphics(fig)
        close(fig);
    end
end

function safe_close_video(vid)
%SAFE_CLOSE_VIDEO  Close a VideoWriter, ignoring errors (onCleanup target)
    try
        close(vid);
    catch
        % Already closed or never opened; nothing to release.
    end
end

function cdata = normalize_anim_frame_cdata(cdata, nh, nw)
%NORMALIZE_ANIM_FRAME_CDATA  Resize RGB(A) frame to exact [nh x nw] for VideoWriter.
%  On Windows Hi-DPI, getframe often returns physical pixels larger than figure Position;
%  without this, writeVideo fails ("frames must be WxH").
    if isempty(cdata)
        error('render_motion_animation:EmptyAnimFrame', 'getframe returned empty cdata');
    end
    [h0, w0, ~] = size(cdata);
    if h0 == nh && w0 == nw
        return;
    end
    if exist('imresize', 'file') == 2
        cdata = imresize(cdata, [nh, nw]);
        return;
    end
    nc = size(cdata, 3);
    x = single(1:w0);
    y = single(1:h0);
    [Xi, Yi] = meshgrid(linspace(1, w0, nw), linspace(1, h0, nh));
    out = zeros(nh, nw, nc, 'uint8');
    for k = 1:nc
        Zi = interp2(x, y, single(cdata(:, :, k)), Xi, Yi, 'linear', 0);
        out(:, :, k) = im2uint8_mat(Zi);
    end
    cdata = out;
end

function [out_h, out_w] = fit_anim_frame_size(h0, w0, max_h, max_w)
%FIT_ANIM_FRAME_SIZE  Return an even-sized frame within codec limits.
    if h0 <= 0 || w0 <= 0
        error('render_motion_animation:InvalidAnimFrameSize', ...
              'Invalid animation frame size [%d x %d].', w0, h0);
    end

    scale = min([1, max_w / w0, max_h / h0]);
    out_w = max(2, 2 * floor((w0 * scale) / 2));
    out_h = max(2, 2 * floor((h0 * scale) / 2));

    % Enforce bounds again after even rounding.
    out_w = min(out_w, 2 * floor(max_w / 2));
    out_h = min(out_h, 2 * floor(max_h / 2));
end

function u8 = im2uint8_mat(Z)
    Z = max(0, min(255, round(Z)));
    u8 = uint8(Z);
end
