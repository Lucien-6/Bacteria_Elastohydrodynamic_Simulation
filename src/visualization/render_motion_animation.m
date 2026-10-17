function info = render_motion_animation(results, params, output_dir, opts)
%RENDER_MOTION_ANIMATION  Render the 4-panel motion animation (overview + co-moving windows)
%
%  Description:
%    Renders the bacterium motion animation from a finished simulation and
%    writes it to a video file in output_dir.  Fixed 16:9 canvas
%    (2560 x 1440 logical px) with four panels:
%      a (left)         : 3-D overview (orthographic, azimuth 45 deg, elevation
%                         25 deg), fixed limits, gravity arrow
%      b (centre)       : XZ projection (view along Y), fixed limits,
%                         2-um scale bar
%      c (right, top)   : co-moving 3-D window (orthographic, same view
%                         direction as panel a), centred on the centre of mass,
%                         with ticks, tick labels and grid (lab-frame coordinates)
%      d (right, bottom): co-moving XZ window (orthographic), centred on the
%                         centre of mass, with ticks, tick labels and grid
%                         (lab-frame coordinates)
%    The axes of panels a and b are sized to the rendered data box and every
%    panel title sits flush on the upper-left corner of its image.  Panel d has a
%    fixed half-width of 1.15 x the largest distance from the centre of mass to
%    any body / pilus point over all frames (at least 1.1 x the scale-bar
%    length); panel c shows a cube of half-width 1.05 x that distance, zoomed to
%    fill its axes.  Each window keeps constant magnification and carries a
%    2-um scale bar.  Ticks and grid of c and d follow the lab frame, so the grid
%    slides past the bacterium.  Gravity arrow: panel a only.
%
%    Every panel shows the semi-transparent body surface, every pilus, the
%    centre-of-mass marker and (panels a and b only) the full-resolution
%    centre-of-mass trajectory up to the current frame.  A global time stamp
%    (simulated time and playback speed relative to real time) and a legend,
%    both horizontally centred on the canvas, and the panel titles are drawn on
%    the canvas.
%
%    Rendering strategy: every graphics object (body, pili, markers,
%    trajectory, border lines, scale bars, lights) is created ONCE; each frame
%    only updates the XData/YData/ZData (or Vertices) of the existing objects,
%    the limits / camera of the two co-moving windows and the time stamp
%    before drawnow / getframe / writeVideo.  This replaces the former
%    per-frame cla + re-creation of all objects, which dominated the render
%    time.
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
%                   .scale_bar_um   Length of the scale bars [um]
%                                   (default: 2)
%
%  Outputs:
%    info - Structure with fields:
%             .file           Full path of the video file
%             .filename       File name (with extension)
%             .n_frames       Number of rendered frames
%             .elapsed_s      Total wall time of this function [s]
%             .stage_s        Structure of per-stage wall times [s]:
%                               prepass, setup, update, drawnow, getframe,
%                               resize, write
%             .window_half_um Half-width of the co-moving windows [um]
%             .playback_speed Simulated seconds per video second (NaN when
%                             only one frame is rendered)
%
%  Notes:
%    The capsule mesh resolution and the per-pilus colour palette below
%    intentionally mirror the STATIC VISUALISATION section of run_simulation.m
%    and run_simulation_HPC_Linux.m (snapshots of the summary figure); keep them
%    in sync when either changes.  The body and trajectory colours are defined
%    only here.
%
%    Font sizes are chosen for the 2560 x 1440 canvas (they remain legible when
%    the video is viewed at 1080p): ticks 16 pt, axis labels 18 pt, legend
%    16 pt, time stamp 24 pt, panel letters 24 pt bold.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-17
%  Updated: 2026-10-18 — 4-panel layout (YZ panel replaced by two co-moving
%           windows), legend, scale bars, global time stamp with playback
%           speed, 16:9 canvas, larger fonts and line widths;
%           co-moving windows with axes boxes and no trajectory in panel d,
%           orthographic panel a filling its axes, titles flush on the image
%           corners, centred time stamp / legend, manually drawn legend;
%           co-moving windows with ticks, tick labels, grid and axis labels,
%           no gravity arrows in panels c and d

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
scale_bar_um   = get_opt(opts, 'scale_bar_um', 2);

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
    r_body_um = max(vecnorm(capsule_flat, 2, 1));   % Max distance body point - CM [um]
else
    custom_verts_um = params.body_verts * 1e6;   % [V x 3] in um
    capsule_flat = [];
    mesh_size    = [];
    a_um         = [];
    r_body_um = max(vecnorm(custom_verts_um, 2, 2));
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
ink_color   = [0.15 0.15 0.15];          % Scale bars, gravity arrows, frames

% --- Typography and line widths (tuned for the 2560 x 1440 canvas) ---
fs_tick   = 16;    % Tick labels [pt]
fs_label  = 18;    % Axis labels [pt]
fs_legend = 16;    % Legend [pt]
fs_note   = 16;    % Scale-bar labels, playback note [pt]
fs_time   = 24;    % Time stamp [pt]
fs_letter = 24;    % Panel letters [pt]
lw_axes   = 1.2;   % Axes and frames [pt]
lw_pilus  = 2.0;   % Pilus lines [pt]
lw_traj   = 1.5;   % Trajectory line [pt]
lw_bar    = 4.0;   % Scale bars [pt]

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

% Playback speed: simulated seconds per video second
if n_frames > 1
    playback_speed = (results.t(frame_indices(end)) - results.t(frame_indices(1))) ...
                     / (n_frames - 1) * animation_fps;
else
    playback_speed = NaN;
end

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
r_pili_um = 0;    % Max distance pilus node - CM over all frames [um]
for k = 1:M_pili
    Xp_all = results.X_pili_um{k}(:, :, frame_indices);               % [Nn x 3 x n_frames]
    mn = min(Xp_all, [], [1 3]);
    mx = max(Xp_all, [], [1 3]);
    x_min = min(x_min, mn(1));  x_max = max(x_max, mx(1));
    y_min = min(y_min, mn(2));  y_max = max(y_max, mx(2));
    z_min = min(z_min, mn(3));  z_max = max(z_max, mx(3));
    rel_k = Xp_all - permute(cm_traj_um, [3 1 2]);                   % Relative to CM [Nn x 3 x n_frames]
    r_pili_um = max(r_pili_um, sqrt(max(sum(rel_k.^2, 2), [], 'all')));
end
clear Xp_all rel_k

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

% Co-moving windows: fixed half-width around the centre of mass
% (minimum: 1.1 x scale-bar length, so the scale bar always fits inside)
win_half_um = max(1.15 * max(r_body_um, r_pili_um), 1.1 * scale_bar_um);
logf('[ANIM] Co-moving windows: half-width %.2f um (max CM distance %.2f um)\n', ...
    win_half_um, max(r_body_um, r_pili_um));
stage_s.prepass = toc(t_stage);

%% ---- Canvas layout (logical pixels, origin bottom-left) ----
t_stage = tic;
anim_fig_w = 2560;   % 16:9 canvas
anim_fig_h = 1440;

anim_max_w = 3840;   % safe upper bound for H.264 width  (even)
anim_max_h = 2160;   % safe upper bound for H.264 height (even)

%% ---- Figure (created once, invisible) ----
fig_anim = figure('Units', 'pixels', 'Position', [50, 50, anim_fig_w, anim_fig_h], ...
                  'Color', 'w', 'Visible', 'off');
cleanup_fig = onCleanup(@() safe_close_figure(fig_anim));
set(fig_anim, 'DefaultAxesFontName',  'Arial', ...
              'DefaultAxesFontSize',  fs_tick, ...
              'DefaultAxesLineWidth', lw_axes);

% --- Legend geometry (drawn manually on the canvas for full spacing control) ---
leg_labels = [{'Body'}, arrayfun(@(k) sprintf('Pilus %d', k), 1:M_pili, ...
              'UniformOutput', false), {'Centre of mass', 'Trajectory'}];
n_leg         = numel(leg_labels);
leg_key_w     = 46;      % Width of the key symbol [px]
leg_key_gap   = 12;      % Gap between key and label [px]
leg_entry_gap = 90;      % Gap between neighbouring entries [px]
leg_row_h     = 40;      % Row pitch [px]
leg_max_row_w = 2300;    % Maximum width of one legend row [px]
ax_tmp = axes(fig_anim, 'Units', 'pixels', 'Position', [0 0 anim_fig_w anim_fig_h], ...
              'XLim', [0 anim_fig_w], 'YLim', [0 anim_fig_h], 'Visible', 'off');
leg_w = zeros(1, n_leg);                   % Entry widths (key + gap + label) [px]
for ie = 1:n_leg
    h_tmp = text(ax_tmp, 0, 0, leg_labels{ie}, 'FontName', 'Arial', ...
                 'FontSize', fs_legend, 'Interpreter', 'none');
    ext = get(h_tmp, 'Extent');
    leg_w(ie) = leg_key_w + leg_key_gap + ext(3);
end
delete(ax_tmp);
leg_row_of = ones(1, n_leg);               % Row index of every entry (greedy wrap)
cur_w = 0;
for ie = 1:n_leg
    add_w = leg_w(ie) + (cur_w > 0) * leg_entry_gap;
    if cur_w > 0 && cur_w + add_w > leg_max_row_w
        leg_row_of(ie:end) = leg_row_of(ie) + 1;
        cur_w = 0;
        add_w = leg_w(ie);
    end
    cur_w = cur_w + add_w;
end
leg_rows = max(leg_row_of);
leg_h_px = leg_rows * leg_row_h + 24;

top_strip  = 190;                          % Time stamp, playback note, panel letters
ax_bot     = 110 + leg_h_px;               % Tick labels + axis label + legend strip
ax_top     = anim_fig_h - top_strip;
ax_h       = ax_top - ax_bot;
gap_win    = 120;                          % Gap between the two co-moving windows
win_side   = floor((ax_h - gap_win) / 2);  % Square co-moving windows
win_x      = anim_fig_w - 40 - win_side;   % Left edge of the right column
win_y_top  = ax_top - win_side;            % Bottom of the upper window
win_y_bot  = ax_bot;                       % Bottom of the lower window
% Panels a and b: the axes rectangle is sized to the rendered data box
% (equal aspect, orthographic projection), top-aligned and centred in its
% slot, so the box fills the axes exactly and the panel title can sit flush on
% its upper-left corner.
slot_w     = 810;                          % Width of the slots of panels a and b [px]
slot_a_x   = 100;                          % Left edge of the slot of panel a
slot_b_x   = 1090;                         % Left edge of the slot of panel b
az_v = 45;  el_v = 25;                     % View of panel a and c
dx_l = diff(ax_lim_x);  dy_l = diff(ax_lim_y);  dz_l = diff(ax_lim_z);
proj_w = dx_l * abs(cosd(az_v)) + dy_l * abs(sind(az_v));
proj_h = (dx_l * abs(sind(az_v)) + dy_l * abs(cosd(az_v))) * abs(sind(el_v)) ...
         + dz_l * abs(cosd(el_v));
ppu_a  = min(slot_w / proj_w, ax_h / proj_h);   % Pixels per um, panel a
w_a    = proj_w * ppu_a;   h_a = proj_h * ppu_a;
ppu_b  = min(slot_w / dx_l, ax_h / dz_l);       % Pixels per um, panel b
w_b    = dx_l * ppu_b;     h_b = dz_l * ppu_b;
pos_3d_px  = [slot_a_x + (slot_w - w_a) / 2, ax_top - h_a, w_a, h_a];   % Panel a
pos_xz_px  = [slot_b_x + (slot_w - w_b) / 2, ax_top - h_b, w_b, h_b];   % Panel b
pos_3df_px = [win_x, win_y_top, win_side, win_side];   % Panel c
pos_xzf_px = [win_x, win_y_bot, win_side, win_side];   % Panel d
px2n = @(p) [p(1)/anim_fig_w, p(2)/anim_fig_h, p(3)/anim_fig_w, p(4)/anim_fig_h];

logf('[ANIM] Figure size (logical): %d x %d px  (a: 3D, b: XZ, c/d: co-moving 3D/XZ windows)\n', ...
    anim_fig_w, anim_fig_h);

%% ---- Axes (created once) ----
ax_3d  = axes(fig_anim, 'Units', 'normalized', 'PositionConstraint', 'innerposition', ...
              'Position', px2n(pos_3d_px));
ax_xz  = axes(fig_anim, 'Units', 'normalized', 'PositionConstraint', 'innerposition', ...
              'Position', px2n(pos_xz_px));
ax_3df = axes(fig_anim, 'Units', 'normalized', 'PositionConstraint', 'innerposition', ...
              'Position', px2n(pos_3df_px));
ax_xzf = axes(fig_anim, 'Units', 'normalized', 'PositionConstraint', 'innerposition', ...
              'Position', px2n(pos_xzf_px));
ax_all = [ax_3d, ax_xz, ax_3df, ax_xzf];     % 1: 3D, 2: XZ, 3: co-moving 3D, 4: co-moving XZ

cm_1 = results.X_cm_um(:, 1);

% Panel a: 3-D overview (orthographic, manual camera so that the box fills the axes)
configure_overview_3d(ax_3d, ax_lim_x, ax_lim_y, ax_lim_z, az_v, el_v, min(proj_w, proj_h));
set(ax_3d, 'FontSize', fs_tick);
xlabel(ax_3d, 'X (μm)', 'FontName', 'Arial', 'FontSize', fs_label);
ylabel(ax_3d, 'Y (μm)', 'FontName', 'Arial', 'FontSize', fs_label);
zlabel(ax_3d, 'Z (μm)', 'FontName', 'Arial', 'FontSize', fs_label);

% Panel b: XZ overview (view along Y axis)
configure_axes_view(ax_xz, ax_lim_x, ax_lim_y, ax_lim_z, 0, 0, 'off', true);
set(ax_xz, 'FontSize', fs_tick);
xlabel(ax_xz, 'X (μm)', 'FontName', 'Arial', 'FontSize', fs_label);
zlabel(ax_xz, 'Z (μm)', 'FontName', 'Arial', 'FontSize', fs_label);

% Panel c: co-moving 3-D window (orthographic; manual camera follows the CM;
% axes box without ticks).  The cube half-width cube_half_um contains every
% body / pilus point; the camera zoom makes the projected cube fill the axes.
cube_half_um = max(1.05 * max(r_body_um, r_pili_um), 0.85 * scale_bar_um);
cam_dist = 10 * cube_half_um;
[cam_dir, ~, cam_right, cam_up_s, view_half_um] = ...
    setup_follow_3d(ax_3df, cm_1(:).', cube_half_um, az_v, el_v, cam_dist, lw_axes, ...
                    fs_tick, fs_label);

% Panel d: co-moving XZ window (orthographic, view along Y; axes box without ticks)
setup_follow_xz(ax_xzf, cm_1(:).', win_half_um, ax_lim_y, lw_axes, fs_tick, fs_label);

drawnow;
add_gravity_indicator_outside(fig_anim, ax_3d, 20, 3.5, 1.5);

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
% Per axis: body, pili, centre-of-mass marker, trajectory (panels a, b only),
% [XZ border lines / scale bars], then lights.  Only their data are updated
% inside the frame loop.
n_ax   = numel(ax_all);
h_body = gobjects(1, n_ax);
h_pili = gobjects(M_pili, n_ax);
h_mark = gobjects(1, n_ax);
h_traj = gobjects(1, n_ax);

R_1  = quat_to_rotmat(results.q(:, 1));

for ai = 1:n_ax
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
                              '-', 'Color', pili_colors(k, :), 'LineWidth', lw_pilus);
    end

    h_mark(ai) = plot3(ax, cm_1(1), cm_1(2), cm_1(3), ...
                       'r+', 'MarkerSize', 9, 'LineWidth', 1.5);

    % Trajectory: nothing visible on the first frame (NaN point; plot3 with
    % empty data returns no handle), grows with the frame index afterwards.
    % Drawn in the overview panels a and b only (not in the co-moving
    % windows c and d).
    if ai <= 2
        h_traj(ai) = plot3(ax, NaN, NaN, NaN, ...
                           '-', 'Color', traj_color, 'LineWidth', lw_traj);
    end

    hold(ax, 'off');
end

% --- Panel a: lights, lighting, material, final view ---
camlight(ax_3d, 'headlight');
camlight(ax_3d, 'left');
lighting(ax_3d, 'gouraud');
material(ax_3d, 'dull');
configure_overview_3d(ax_3d, ax_lim_x, ax_lim_y, ax_lim_z, az_v, el_v, min(proj_w, proj_h));

% --- Panel b: light, material, view, border lines, scale bar ---
camlight(ax_xz, 'headlight');
lighting(ax_xz, 'gouraud');
material(ax_xz, 'dull');
configure_axes_view(ax_xz, ax_lim_x, ax_lim_y, ax_lim_z, 0, 0, 'off', true);
hold(ax_xz, 'on');
plot3(ax_xz, [ax_lim_x(2) ax_lim_x(2)], ax_lim_y([1 1]), ax_lim_z, ...
    'k-', 'LineWidth', lw_axes, 'HandleVisibility', 'off');
plot3(ax_xz, ax_lim_x, ax_lim_y([1 1]), [ax_lim_z(2) ax_lim_z(2)], ...
    'k-', 'LineWidth', lw_axes, 'HandleVisibility', 'off');
sb_x0 = ax_lim_x(1) + 0.05 * diff(ax_lim_x);
sb_z0 = ax_lim_z(1) + 0.06 * diff(ax_lim_z);
plot3(ax_xz, sb_x0 + [0 scale_bar_um], ax_lim_y([1 1]), sb_z0 * [1 1], ...
    '-', 'Color', ink_color, 'LineWidth', lw_bar, 'Clipping', 'off', ...
    'HandleVisibility', 'off');
text(ax_xz, sb_x0 + scale_bar_um / 2, ax_lim_y(1), sb_z0 + 0.012 * diff(ax_lim_z), ...
    sprintf('%g μm', scale_bar_um), 'FontName', 'Arial', 'FontSize', fs_note, ...
    'Color', ink_color, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'Clipping', 'off', 'HandleVisibility', 'off', ...
    'BackgroundColor', 'w', 'Margin', 1);
hold(ax_xz, 'off');

% --- Panel c: light, material, scale bar (co-moving 3-D window) ---
camlight(ax_3df, 'headlight');
camlight(ax_3df, 'left');
lighting(ax_3df, 'gouraud');
material(ax_3df, 'dull');
hold(ax_3df, 'on');
h_sb_3df = plot3(ax_3df, NaN, NaN, NaN, '-', 'Color', ink_color, 'LineWidth', lw_bar, ...
                 'Clipping', 'off', 'HandleVisibility', 'off');
h_sbt_3df = text(ax_3df, 0, 0, 0, sprintf('%g μm', scale_bar_um), ...
                 'FontName', 'Arial', 'FontSize', fs_note, 'Color', ink_color, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                 'Clipping', 'off', 'HandleVisibility', 'off', ...
                 'BackgroundColor', 'w', 'Margin', 1);
hold(ax_3df, 'off');

% --- Panel d: light, material, scale bar (co-moving XZ window) ---
camlight(ax_xzf, 'headlight');
lighting(ax_xzf, 'gouraud');
material(ax_xzf, 'dull');
hold(ax_xzf, 'on');
h_sb_xzf = plot3(ax_xzf, NaN, NaN, NaN, '-', 'Color', ink_color, 'LineWidth', lw_bar, ...
                 'Clipping', 'off', 'HandleVisibility', 'off');
h_sbt_xzf = text(ax_xzf, 0, 0, 0, sprintf('%g μm', scale_bar_um), ...
                 'FontName', 'Arial', 'FontSize', fs_note, 'Color', ink_color, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                 'Clipping', 'off', 'HandleVisibility', 'off', ...
                 'BackgroundColor', 'w', 'Margin', 1);
hold(ax_xzf, 'off');

%% ---- Canvas overlay (pixel coordinates): letters, time, legend ----
ax_ov = axes(fig_anim, 'Units', 'normalized', 'Position', [0 0 1 1], ...
             'XLim', [0 anim_fig_w], 'YLim', [0 anim_fig_h], 'Visible', 'off', ...
             'HitTest', 'off', 'Clipping', 'off');
hold(ax_ov, 'on');

% Panel titles: flush on the upper-left corner of the corresponding image
title_gap = 10;                            % Gap above the image [px]
panel_pos = {pos_3d_px, pos_xz_px, pos_3df_px, pos_xzf_px};
panel_str = {'{\bf a}  3-D view', '{\bf b}  XZ view', ...
             '{\bf c}  Co-moving 3-D view', '{\bf d}  Co-moving XZ view'};
for ip = 1:4
    pp = panel_pos{ip};
    text(ax_ov, pp(1), pp(2) + pp(4) + title_gap, panel_str{ip}, ...
         'Interpreter', 'tex', 'FontName', 'Arial', 'FontSize', fs_letter, ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% Global time stamp and playback speed (horizontally centred on the canvas)
x_mid = anim_fig_w / 2;
h_time = text(ax_ov, x_mid, anim_fig_h - 45, '', ...
              'Interpreter', 'tex', 'FontName', 'Arial', 'FontSize', fs_time, ...
              'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
              'VerticalAlignment', 'middle');
if isnan(playback_speed)
    speed_str = '';
elseif abs(playback_speed - 1) < 0.005
    speed_str = 'Playback: real time';
elseif playback_speed < 1
    speed_str = sprintf('Playback: %.3g× real time (slow motion)', playback_speed);
else
    speed_str = sprintf('Playback: %.3g× real time (accelerated)', playback_speed);
end
text(ax_ov, x_mid, anim_fig_h - 100, speed_str, ...
     'Interpreter', 'none', 'FontName', 'Arial', 'FontSize', fs_note, ...
     'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'middle');

% Legend: drawn manually (key + label per entry), every row centred on the canvas
traj_key_color = traj_color(1:3) * traj_color(4) + (1 - traj_color(4));  % Trajectory colour on white
y_leg_top = 14 + leg_rows * leg_row_h;     % Top of the legend block [px]
for ir = 1:leg_rows
    idx_r = find(leg_row_of == ir);
    row_w = sum(leg_w(idx_r)) + leg_entry_gap * (numel(idx_r) - 1);
    x_c   = x_mid - row_w / 2;             % Left edge of the first entry
    y_c   = y_leg_top - (ir - 0.5) * leg_row_h;
    for ie = idx_r
        if ie == 1                                     % Body
            patch(ax_ov, x_c + [0 leg_key_w leg_key_w 0], y_c + [-8 -8 8 8], body_color, ...
                  'FaceAlpha', 0.55, 'EdgeColor', 'none');
        elseif ie <= 1 + M_pili                        % Pilus k
            line(ax_ov, x_c + [0 leg_key_w], y_c * [1 1], ...
                 'Color', pili_colors(ie - 1, :), 'LineWidth', lw_pilus);
        elseif ie == 2 + M_pili                        % Centre of mass
            line(ax_ov, x_c + leg_key_w / 2, y_c, 'LineStyle', 'none', 'Marker', '+', ...
                 'Color', 'r', 'MarkerSize', 9, 'LineWidth', 1.5);
        else                                           % Trajectory
            line(ax_ov, x_c + [0 leg_key_w], y_c * [1 1], ...
                 'Color', traj_key_color, 'LineWidth', lw_traj);
        end
        text(ax_ov, x_c + leg_key_w + leg_key_gap, y_c, leg_labels{ie}, ...
             'FontName', 'Arial', 'FontSize', fs_legend, 'Interpreter', 'none', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        x_c = x_c + leg_w(ie) + leg_entry_gap;
    end
end
hold(ax_ov, 'off');

stage_s.setup = toc(t_stage);

%% ---- Frame loop: update data only ----
sb_off    = 0.86 * win_half_um;    % Scale-bar corner offset from window centre, panel d [um]
sb_off_3d = 0.88 * view_half_um;   % Same for panel c (screen coordinates) [um]
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

    for ai = 1:n_ax
        if is_capsule
            set(h_body(ai), 'XData', Xs, 'YData', Ys, 'ZData', Zs);
        else
            set(h_body(ai), 'Vertices', verts);
        end
    end

    for k = 1:M_pili
        Xp = results.X_pili_um{k}(:, :, ti);
        for ai = 1:n_ax
            set(h_pili(k, ai), 'XData', Xp(:,1), 'YData', Xp(:,2), 'ZData', Xp(:,3));
        end
    end

    for ai = 1:n_ax
        set(h_mark(ai), 'XData', cm_um(1), 'YData', cm_um(2), 'ZData', cm_um(3));
    end

    if fi > 1
        for ai = [1 2]
            set(h_traj(ai), 'XData', cm_traj_full_um(1, 1:ti), ...
                            'YData', cm_traj_full_um(2, 1:ti), ...
                            'ZData', cm_traj_full_um(3, 1:ti));
        end
    end

    % --- Co-moving windows follow the centre of mass ---
    cm_row = cm_um(:).';
    % Panel c: limits + manual orthographic camera (target = CM)
    set(ax_3df, 'XLim', cm_row(1) + cube_half_um * [-1 1], ...
                'YLim', cm_row(2) + cube_half_um * [-1 1], ...
                'ZLim', cm_row(3) + cube_half_um * [-1 1], ...
                'CameraTarget', cm_row, ...
                'CameraPosition', cm_row + cam_dist * cam_dir);
    p0 = cm_row - sb_off_3d * cam_right + sb_off_3d * cam_up_s;   % Scale-bar start (screen upper-left)
    p1 = p0 + scale_bar_um * cam_right;
    set(h_sb_3df, 'XData', [p0(1) p1(1)], 'YData', [p0(2) p1(2)], 'ZData', [p0(3) p1(3)]);
    pt = (p0 + p1) / 2 - 0.03 * view_half_um * cam_up_s;          % Label below the bar
    set(h_sbt_3df, 'Position', pt);
    % Panel d: limits only (view along Y, orthographic)
    set(ax_xzf, 'XLim', cm_row(1) + win_half_um * [-1 1], ...
                'ZLim', cm_row(3) + win_half_um * [-1 1]);
    set(h_sb_xzf, 'XData', cm_row(1) - sb_off + [0 scale_bar_um], ...
                  'YData', ax_lim_y([1 1]), ...
                  'ZData', (cm_row(3) - sb_off) * [1 1]);
    set(h_sbt_xzf, 'Position', [cm_row(1) - sb_off + scale_bar_um / 2, ax_lim_y(1), ...
                                cm_row(3) - sb_off + 0.03 * win_half_um]);

    set(h_time, 'String', sprintf('{\\itt} = %.4f s', results.t(ti)));
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
if ~isnan(playback_speed)
    logf('[ANIM] Playback speed: %.4g simulated s per video s (%.3g x real time)\n', ...
        playback_speed, playback_speed);
end

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
              'stage_s', stage_s, 'window_half_um', win_half_um, ...
              'playback_speed', playback_speed);

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

function configure_overview_3d(ax, lim_x, lim_y, lim_z, az, el, proj_short)
%CONFIGURE_OVERVIEW_3D  Orthographic 3-D overview whose projected box fills the axes
%  Applies the fixed limits, equal aspect, grid, box and view, then places the
%  camera manually.  MATLAB maps CameraViewAngle to the SHORTER side of the
%  axes, so the field of view is set to proj_short [um], the shorter side of the
%  projected data box; a box of that aspect ratio then fills the axes exactly.
    configure_axes_view(ax, lim_x, lim_y, lim_z, az, el, 'on', true);
    ctr  = [mean(lim_x), mean(lim_y), mean(lim_z)];
    pos0 = get(ax, 'CameraPosition');
    tgt0 = get(ax, 'CameraTarget');
    cam_dir = (pos0 - tgt0) / norm(pos0 - tgt0);
    cam_up  = get(ax, 'CameraUpVector');
    cam_up  = cam_up / norm(cam_up);
    cam_dist = 10 * norm([diff(lim_x), diff(lim_y), diff(lim_z)]);
    set(ax, 'CameraTarget', ctr, ...
            'CameraPosition', ctr + cam_dist * cam_dir, ...
            'CameraUpVector', cam_up, ...
            'CameraViewAngle', 2 * atand(proj_short / 2 / cam_dist));
end

function [cam_dir, cam_up, cam_right, cam_up_s, view_half] = ...
        setup_follow_3d(ax, cm_row, cube_half, az, el, cam_dist, lw, fs_tick, fs_label)
%SETUP_FOLLOW_3D  Orthographic co-moving 3-D window centred on cm_row
%  The axes box is a cube of half-width cube_half [um] around cm_row, with
%  ticks, tick labels, grid and axis labels (lab-frame coordinates).  The
%  camera is placed manually and zoomed so that the projected cube just fills
%  the square axes.
%  Outputs: unit vector target -> camera and camera up vector (as defined by
%  view(ax, az, el)), the screen right / up unit vectors in data space, and the
%  half-width of the visible square [um] (orthographic scale).
    xlim(ax, cm_row(1) + cube_half * [-1 1]);
    ylim(ax, cm_row(2) + cube_half * [-1 1]);
    zlim(ax, cm_row(3) + cube_half * [-1 1]);
    set(ax, 'DataAspectRatio', [1 1 1], 'Projection', 'orthographic');
    view(ax, az, el);
    pos0 = get(ax, 'CameraPosition');
    tgt0 = get(ax, 'CameraTarget');
    cam_dir = (pos0 - tgt0) / norm(pos0 - tgt0);
    cam_up  = get(ax, 'CameraUpVector');
    cam_up  = cam_up / norm(cam_up);
    cam_fwd   = -cam_dir;
    cam_right = cross(cam_fwd, cam_up);    cam_right = cam_right / norm(cam_right);
    cam_up_s  = cross(cam_right, cam_fwd); cam_up_s  = cam_up_s / norm(cam_up_s);
    % Half-extent of the projected cube along screen right / up (+2 % margin)
    view_half = 1.02 * cube_half * max(sum(abs(cam_right)), sum(abs(cam_up_s)));
    set(ax, 'CameraTarget', cm_row, ...
            'CameraPosition', cm_row + cam_dist * cam_dir, ...
            'CameraUpVector', cam_up, ...
            'CameraViewAngle', 2 * atand(view_half / cam_dist));
    set_window_axes(ax, lw, fs_tick, fs_label, true);
end

function setup_follow_xz(ax, cm_row, half_um, lim_y, lw, fs_tick, fs_label)
%SETUP_FOLLOW_XZ  Orthographic co-moving XZ window (view along Y) centred on cm_row
%  Axes box with ticks, tick labels, grid and axis labels (lab-frame
%  coordinates).  The Y range stays fixed to the global range so that nothing is
%  clipped along the viewing direction.
    xlim(ax, cm_row(1) + half_um * [-1 1]);
    ylim(ax, lim_y);
    zlim(ax, cm_row(3) + half_um * [-1 1]);
    set(ax, 'DataAspectRatio', [1 1 1], 'Projection', 'orthographic');
    view(ax, 0, 0);
    set_window_axes(ax, lw, fs_tick, fs_label, false);
end

function set_window_axes(ax, lw, fs_tick, fs_label, with_y_label)
%SET_WINDOW_AXES  Axes box, ticks (lab-frame coordinates), tick labels, grid, axis labels
%  The tick positions are automatic, so they are recomputed when the limits
%  follow the centre of mass; the grid therefore slides past the bacterium.
    set(ax, 'Box', 'on', 'TickDir', 'out', 'LineWidth', lw, 'FontSize', fs_tick, ...
            'Layer', 'top');
    grid(ax, 'on');
    xlabel(ax, 'X (μm)', 'FontName', 'Arial', 'FontSize', fs_label);
    if with_y_label
        ylabel(ax, 'Y (μm)', 'FontName', 'Arial', 'FontSize', fs_label);
    end
    zlabel(ax, 'Z (μm)', 'FontName', 'Arial', 'FontSize', fs_label);
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
