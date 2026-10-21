function info = render_motion_animation(results, params, output_dir, opts)
%RENDER_MOTION_ANIMATION  Render the 4-panel motion animation (overview + co-moving windows)
%
%  Description:
%    Renders the bacterium motion animation from a finished simulation and
%    writes it to a video file in output_dir.  Fixed 16:9 canvas
%    (2560 x 1440 logical px) with four panels:
%      a : 3-D overview (orthographic, azimuth 45 deg, elevation 25 deg),
%          fixed limits, gravity arrow
%      b : XZ projection (view along Y), fixed limits, 2-um scale bar
%      c : co-moving 3-D window (orthographic, same view direction as panel a),
%          centred on the centre of mass, with ticks, tick labels and grid
%          (lab-frame coordinates)
%      d : co-moving XZ window (orthographic), centred on the centre of mass,
%          with ticks, tick labels and grid (lab-frame coordinates)
%    The arrangement of the four panels on the canvas is chosen automatically
%    (see below), so the canvas is used as completely as possible.
%    Adaptive layout: every axes rectangle has the aspect ratio of its data box
%    (equal scale on all axes) and a natural size (panels a, b: 810 x 1076 px,
%    panels c, d: 440 x 478 px).  All guillotine arrangements (nested
%    side-by-side / stacked splits) of the panels in the reading orders a-b-c-d,
%    a-c-b-d and a-c-d-b are tried; the one that fits the largest COMMON scale
%    factor of the four rectangles into the free canvas region wins, so the
%    relative sizes of the panels never change.  The slack is spread over the
%    gaps (up to 120 px) and the group is centred.  The chosen arrangement is
%    logged, e.g. H(a,H(b,V(c,d))).  Every panel title sits flush on the
%    upper-left corner of its image.  The
%    co-moving windows c and d are NOT cubes / squares: along every lab axis
%    their extent is the reach of the body and all pili over all frames
%    (relative to the centre of mass) plus a margin of 8 % (at least 0.1 um),
%    so the bacterium fills the window and its shape stays visible; the scale
%    is equal on the X, Y and Z axes (nothing is stretched), the axes rectangle
%    of each window has the aspect ratio of its box, and the box is centred on
%    the middle of the occupied range (constant offset from the centre of mass).
%    Each window keeps constant magnification (no scale bar; panel b carries
%    the 2-um scale bar).  Ticks and grid
%    of c and d follow the lab frame, so the grid slides past the bacterium.
%    Tick spacing on c and d is never finer than 0.5 μm.
%    Gravity arrow: panel a only.
%
%    The body keeps a uniform light green.  A short dark-green patch sits on one
%    side of the waist and rotates with the body, so a roll about the long
%    axis stays visible without a line or a half-and-half split (custom
%    bodies: long axis = principal axis of the vertex cloud).
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
%                   .scale_bar_um   Length of the scale bar of panel b [um]
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
%             .window_half_um Half-extents [hx hy hz] of the co-moving box [um]
%             .window_offset_um Offset of the box centre from the centre of
%                             mass [1 x 3, um]
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
%  Updated: 2026-10-19 — v8.6.3: co-moving windows sized per axis from the
%           reach of the bacterium (non-cubic box, equal scale on all axes);
%           no scale bar in panels c and d; two dark body stripes make the
%           roll visible; adaptive panel layout (best of all guillotine
%           arrangements, common scale factor)
%  Updated: 2026-10-21 — co-moving windows c and d use a minimum tick
%           spacing of 0.5 μm; body roll is shown by a short dark-green patch
%           on one side of the waist

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
ink_color   = [0.15 0.15 0.15];          % Scale bars, gravity arrows, frames

% --- Roll cue: a short patch on one side of the waist ---
% A uniform body hides rotation about the long axis.  One dark-green patch,
% shorter than the body and only a small fraction of the circumference, is
% fixed in the body frame and therefore rotates with it.
body_color_dark = [0.03 0.25 0.10];
waist_alpha     = 0.95;
if is_capsule
    e_long  = [1 0 0];
    e_side  = [0 1 0];
    waist_origin = [0 0 0];
    waist_pts = capsule_flat.';
else
    waist_origin = mean(custom_verts_um, 1);
    waist_pts = custom_verts_um - waist_origin;
    [Vp, Dp] = eig(cov(waist_pts));
    [~, i_long] = max(diag(Dp));
    e_long = Vp(:, i_long).';
    e_tmp = [0 1 0];
    if abs(dot(e_tmp, e_long)) > 0.9
        e_tmp = [0 0 1];
    end
    e_side = e_tmp - dot(e_tmp, e_long) * e_long;
    e_side = e_side / norm(e_side);
end
[waist_flat, waist_size] = waist_mark_surface(waist_pts, e_long, e_side);
waist_flat = waist_flat + waist_origin.';
clear waist_pts e_long e_side e_tmp waist_origin i_long Vp Dp

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
rel_min = [0 0 0];    % Reach of body + pili relative to the CM per lab axis over all frames [um]
rel_max = [0 0 0];    % (the CM itself is included, so the CM marker is always inside)
for idx = 1:n_frames
    ti   = frame_indices(idx);
    R_fi = quat_to_rotmat(results.q(:, ti));
    if is_capsule
        pts_rel = R_fi * capsule_flat;                               % [3 x N_mesh], relative to CM
    else
        pts_rel = R_fi * custom_verts_um';                           % [3 x N_verts], relative to CM
    end
    rel_min = min(rel_min, min(pts_rel, [], 2).');
    rel_max = max(rel_max, max(pts_rel, [], 2).');
    pts = pts_rel + results.X_cm_um(:, ti);
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
    rel_k = Xp_all - permute(cm_traj_um, [3 1 2]);                   % Relative to CM [Nn x 3 x n_frames]
    rel_min = min(rel_min, min(rel_k, [], [1 3]));
    rel_max = max(rel_max, max(rel_k, [], [1 3]));
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

% Co-moving windows: a fixed, generally NON-cubic box that follows the centre
% of mass.  Its extent along every lab axis is the reach of the body + pili
% over all frames (rel_min / rel_max) plus a small margin, so the bacterium
% fills the box and its shape stays visible; the unit scale is equal on all
% three axes (no stretching).  The box is centred on win_off_um (constant
% offset from the CM) so that it is symmetric about the occupied range.
rel_span = rel_max - rel_min;
win_pad  = max(0.08 * rel_span, 0.10);                  % Margin per axis [um]
win_half_um = rel_span / 2 + win_pad;                   % [hx hy hz] half-extents [um]
win_off_um  = (rel_min + rel_max) / 2;                  % Box centre relative to the CM [um]
logf(['[ANIM] Co-moving windows (box, equal scale): half-extents X=%.2f, Y=%.2f, Z=%.2f um; ', ...
      'centre offset from CM (%.2f, %.2f, %.2f) um\n'], win_half_um, win_off_um);
win_tick_min_um = 0.5;   % Finest tick spacing on panels c and d [μm]
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

% Panel titles (also needed here to measure their widths for the layout)
panel_str = {'{\bf a}  3-D view', '{\bf b}  XZ view', ...
             '{\bf c}  Co-moving 3-D view', '{\bf d}  Co-moving XZ view'};

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
title_w = zeros(1, 4);                     % Rendered widths of the panel titles [px]
for ip = 1:4
    h_tmp = text(ax_tmp, 0, 0, panel_str{ip}, 'Interpreter', 'tex', 'FontName', 'Arial', ...
                 'FontSize', fs_letter);
    ext = get(h_tmp, 'Extent');
    title_w(ip) = ext(3);
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

% ---- Panel geometry (natural sizes) ----
az_v = 45;  el_v = 25;                     % View of panel a and c
dx_l = diff(ax_lim_x);  dy_l = diff(ax_lim_y);  dz_l = diff(ax_lim_z);
proj_w = dx_l * abs(cosd(az_v)) + dy_l * abs(sind(az_v));            % Projected box, panel a [um]
proj_h = (dx_l * abs(sind(az_v)) + dy_l * abs(cosd(az_v))) * abs(sind(el_v)) ...
         + dz_l * abs(cosd(el_v));
dx_w = 2 * win_half_um(1);  dy_w = 2 * win_half_um(2);  dz_w = 2 * win_half_um(3);
proj_w_c = dx_w * abs(cosd(az_v)) + dy_w * abs(sind(az_v));          % Projected box, panel c [um]
proj_h_c = (dx_w * abs(sind(az_v)) + dy_w * abs(cosd(az_v))) * abs(sind(el_v)) ...
           + dz_w * abs(cosd(el_v));
% Natural size of every axes rectangle (equal scale on all axes): panels a and b
% fitted into 810 x 1076 px, panels c and d into 440 x 478 px.  Only the
% aspect ratios and these relative sizes matter: the layout below multiplies
% all four rectangles by ONE common factor, so their relative sizes are kept.
ref_ab = [810, 1076];
ref_cd = [440, 478];
ppu_a = min(ref_ab(1) / proj_w,   ref_ab(2) / proj_h);     % Pixels per um, panel a
ppu_b = min(ref_ab(1) / dx_l,     ref_ab(2) / dz_l);       % Pixels per um, panel b
ppu_c = min(ref_cd(1) / proj_w_c, ref_cd(2) / proj_h_c);   % Pixels per um, panel c
ppu_d = min(ref_cd(1) / dx_w,     ref_cd(2) / dz_w);       % Pixels per um, panel d
lay.rect = [proj_w * ppu_a,   proj_h * ppu_a;
            dx_l * ppu_b,     dz_l * ppu_b;
            proj_w_c * ppu_c, proj_h_c * ppu_c;
            dx_w * ppu_d,     dz_w * ppu_d];                % [w h] of the rectangles a-d [px]
% Space reserved around each rectangle for tick labels, axis labels and the
% title: [left right bottom top] [px].  Panels a and c are 3-D (Y labels on the
% right); panel a also keeps room for the gravity arrow.
lay.marg = [110  90 100 60;
            110  30  90 60;
            110 110 100 60;
            110  30  90 60];
lay.title_w = title_w(:);

% ---- Adaptive layout: choose the arrangement that allows the largest panels ----
% Candidates are all guillotine layouts (nested horizontal / vertical splits)
% of the four panels in the reading orders a-b-c-d, a-c-b-d and a-c-d-b.  For
% every candidate the largest common scale factor s of the four rectangles that
% still fits into the free canvas region is found by bisection; the candidate
% with the largest s wins (ties: smallest bounding box).  The remaining slack is
% spread over the gaps between the panels (up to 120 px) and the whole group is
% centred in the region.
region_x0  = 30;                                 % Free canvas region [px]
region_w   = anim_fig_w - 60;
region_y0  = leg_h_px + 4;                       % Above the legend
region_top = anim_fig_h - 125;                   % Below the time stamp and the playback note
region_h   = region_top - region_y0;
gap_min    = 10;                                 % Gap between neighbouring boxes [px]
gap_max    = 120;
best_s = -Inf;  best_area = Inf;  best_tree = [];
orders = [1 2 3 4; 1 3 2 4; 1 3 4 2];
for io = 1:size(orders, 1)
    cand = enum_layout_trees(orders(io, :));
    for ic = 1:numel(cand)
        s_c = round(max_layout_scale(cand{ic}, gap_min, region_w, region_h, lay), 3);
        [w_c0, h_c0] = layout_size(cand{ic}, s_c, gap_min, lay);
        if s_c > best_s || (s_c == best_s && w_c0 * h_c0 < best_area)
            best_s = s_c;  best_area = w_c0 * h_c0;  best_tree = cand{ic};
        end
    end
end
lay_gap = gap_min;                               % Spread the slack over the gaps
for it = 1:20
    g_mid = (lay_gap + gap_max) / 2;
    [w_t, h_t] = layout_size(best_tree, best_s, g_mid, lay);
    if w_t <= region_w && h_t <= region_h, lay_gap = g_mid; else, gap_max = g_mid; end
end
[w_t, h_t] = layout_size(best_tree, best_s, lay_gap, lay);
rects = place_layout(best_tree, region_x0 + (region_w - w_t) / 2, ...
                     region_top - (region_h - h_t) / 2, w_t, best_s, lay_gap, lay, zeros(4, 4));
pos_3d_px  = rects(1, :);   % Panel a
pos_xz_px  = rects(2, :);   % Panel b
pos_3df_px = rects(3, :);   % Panel c
pos_xzf_px = rects(4, :);   % Panel d
logf('[ANIM] Layout %s: scale %.2f, gap %.0f px (panel a %.0f x %.0f px, c %.0f x %.0f px)\n', ...
    layout_desc(best_tree), best_s, lay_gap, rects(1, 3), rects(1, 4), rects(3, 3), rects(3, 4));
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

% Panel c: co-moving 3-D window (orthographic; manual camera follows the box
% centre).  The box has half-extents win_half_um (per axis) around CM + win_off_um;
% the camera zoom makes the projected box fill the axes rectangle exactly.
cam_dist = 10 * norm(2 * win_half_um);
[cam_dir, ~, ~, ~] = ...
    setup_follow_3d(ax_3df, cm_1(:).' + win_off_um, win_half_um, az_v, el_v, cam_dist, ...
                    min(proj_w_c, proj_h_c), lw_axes, fs_tick, fs_label);

% Panel d: co-moving XZ window (orthographic, view along Y)
setup_follow_xz(ax_xzf, cm_1(:).' + win_off_um, win_half_um, ax_lim_y, lw_axes, fs_tick, fs_label);

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
% Per axis: body, waist mark, pili, centre-of-mass marker, trajectory
% (panels a, b only), [panel b: XZ border lines / scale bar], then lights.
% Only their data are updated inside the frame loop.
n_ax   = numel(ax_all);
h_body = gobjects(1, n_ax);
h_waist = gobjects(1, n_ax);
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
                             'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
                             'Parent', ax);
    end
    waist_1 = R_1 * waist_flat + cm_1;
    h_waist(ai) = surf(ax, reshape(waist_1(1, :), waist_size), ...
                           reshape(waist_1(2, :), waist_size), ...
                           reshape(waist_1(3, :), waist_size), ...
                       'FaceColor', body_color_dark, 'FaceAlpha', waist_alpha, ...
                       'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
                       'HandleVisibility', 'off');

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

% --- Panel c: lights and material (co-moving 3-D window; no scale bar) ---
camlight(ax_3df, 'headlight');
camlight(ax_3df, 'left');
lighting(ax_3df, 'gouraud');
material(ax_3df, 'dull');

% --- Panel d: light and material (co-moving XZ window; no scale bar) ---
camlight(ax_xzf, 'headlight');
lighting(ax_xzf, 'gouraud');
material(ax_xzf, 'dull');

%% ---- Canvas overlay (pixel coordinates): letters, time, legend ----
ax_ov = axes(fig_anim, 'Units', 'normalized', 'Position', [0 0 1 1], ...
             'XLim', [0 anim_fig_w], 'YLim', [0 anim_fig_h], 'Visible', 'off', ...
             'HitTest', 'off', 'Clipping', 'off');
hold(ax_ov, 'on');

% Panel titles: flush on the upper-left corner of the corresponding image
title_gap = 10;                            % Gap above the image [px]
panel_pos = {pos_3d_px, pos_xz_px, pos_3df_px, pos_xzf_px};
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
            patch(ax_ov, x_c + leg_key_w * [0.38 0.62 0.62 0.38], y_c + [-4 -4 4 4], ...
                  body_color_dark, 'FaceAlpha', 0.95, 'EdgeColor', 'none');
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

    waist_lab = R_fi * waist_flat + cm_um;
    for ai = 1:n_ax
        if is_capsule
            set(h_body(ai), 'XData', Xs, 'YData', Ys, 'ZData', Zs);
        else
            set(h_body(ai), 'Vertices', verts);
        end
        set(h_waist(ai), 'XData', reshape(waist_lab(1, :), waist_size), ...
                         'YData', reshape(waist_lab(2, :), waist_size), ...
                         'ZData', reshape(waist_lab(3, :), waist_size));
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
    cm_row  = cm_um(:).';
    box_ctr = cm_row + win_off_um;                % Centre of the co-moving box
    % Panel c: limits + manual orthographic camera (target = box centre)
    lim_c = [box_ctr(1) + win_half_um(1) * [-1 1]; ...
             box_ctr(2) + win_half_um(2) * [-1 1]; ...
             box_ctr(3) + win_half_um(3) * [-1 1]];
    set(ax_3df, 'XLim', lim_c(1, :), 'YLim', lim_c(2, :), 'ZLim', lim_c(3, :), ...
                'XTick', window_ticks(lim_c(1, :), win_tick_min_um), ...
                'YTick', window_ticks(lim_c(2, :), win_tick_min_um), ...
                'ZTick', window_ticks(lim_c(3, :), win_tick_min_um), ...
                'CameraTarget', box_ctr, ...
                'CameraPosition', box_ctr + cam_dist * cam_dir);
    % Panel d: limits only (view along Y, orthographic).  Y is the viewing
    % direction and keeps the global range; X and Z follow the box.
    lim_dx = box_ctr(1) + win_half_um(1) * [-1 1];
    lim_dz = box_ctr(3) + win_half_um(3) * [-1 1];
    set(ax_xzf, 'XLim', lim_dx, 'ZLim', lim_dz, ...
                'XTick', window_ticks(lim_dx, win_tick_min_um), ...
                'ZTick', window_ticks(lim_dz, win_tick_min_um));

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
              'window_offset_um', win_off_um, ...
              'playback_speed', playback_speed);

end


%% ==================== Local Functions ====================

function trees = enum_layout_trees(order)
%ENUM_LAYOUT_TREES  All binary guillotine layouts of the panels in the given order
%  A leaf is a panel index; an inner node is a struct with the split
%  orientation o ('H': side by side, 'V': stacked) and the two sub-layouts.
    n = numel(order);
    if n == 1
        trees = {order};
        return
    end
    trees = {};
    for k = 1:n - 1
        L = enum_layout_trees(order(1:k));
        R = enum_layout_trees(order(k+1:end));
        for il = 1:numel(L)
            for ir = 1:numel(R)
                for o = 'HV'
                    trees{end+1} = struct('o', o, 'l', L(il), 'r', R(ir)); %#ok<AGROW>
                end
            end
        end
    end
end

function [w, h] = layout_size(t, s, g, lay)
%LAYOUT_SIZE  Size [px] of the box of a (sub-)layout at scale s and gap g
%  The box of a panel is its rectangle (scaled by s) plus the label margins;
%  it is at least as wide as the panel title.
    if ~isstruct(t)
        w = max(s * lay.rect(t, 1) + lay.marg(t, 1) + lay.marg(t, 2), ...
                lay.marg(t, 1) + lay.title_w(t) + 10);
        h = s * lay.rect(t, 2) + lay.marg(t, 3) + lay.marg(t, 4);
        return
    end
    [w1, h1] = layout_size(t.l, s, g, lay);
    [w2, h2] = layout_size(t.r, s, g, lay);
    if t.o == 'H'
        w = w1 + w2 + g;   h = max(h1, h2);
    else
        w = max(w1, w2);   h = h1 + h2 + g;
    end
end

function s = max_layout_scale(t, g, reg_w, reg_h, lay)
%MAX_LAYOUT_SCALE  Largest common scale (<= 2) at which the layout fits the region
    s_lo = 0.2;  s_hi = 2;
    [w, h] = layout_size(t, s_hi, g, lay);
    if w <= reg_w && h <= reg_h
        s = s_hi;
        return
    end
    for it = 1:22
        s_mid = (s_lo + s_hi) / 2;
        [w, h] = layout_size(t, s_mid, g, lay);
        if w <= reg_w && h <= reg_h, s_lo = s_mid; else, s_hi = s_mid; end
    end
    s = s_lo;
end

function rects = place_layout(t, x, ytop, wa, s, g, lay, rects)
%PLACE_LAYOUT  Axes rectangles [x y w h] (origin bottom-left) of all panels
%  (x, ytop) is the top-left corner and wa the allocated width of the box of
%  the sub-layout.  Side by side boxes are top-aligned, stacked boxes centred.
    if ~isstruct(t)
        [wo, ~] = layout_size(t, s, g, lay);
        bx = x + (wa - wo) / 2;
        rw = s * lay.rect(t, 1);   rh = s * lay.rect(t, 2);
        rects(t, :) = [bx + lay.marg(t, 1), ytop - lay.marg(t, 4) - rh, rw, rh];
        return
    end
    [w1, h1] = layout_size(t.l, s, g, lay);
    [w2, ~]  = layout_size(t.r, s, g, lay);
    if t.o == 'H'
        rects = place_layout(t.l, x, ytop, w1, s, g, lay, rects);
        rects = place_layout(t.r, x + w1 + g, ytop, w2, s, g, lay, rects);
    else
        wn = max(w1, w2);
        rects = place_layout(t.l, x, ytop, wn, s, g, lay, rects);
        rects = place_layout(t.r, x, ytop - h1 - g, wn, s, g, lay, rects);
    end
end

function str = layout_desc(t)
%LAYOUT_DESC  Text form of a layout, e.g. H(a,V(b,H(c,d)))
    if ~isstruct(t)
        str = char('a' + t - 1);
    else
        str = sprintf('%s(%s,%s)', t.o, layout_desc(t.l), layout_desc(t.r));
    end
end

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

function [cam_dir, cam_up, cam_right, cam_up_s] = ...
        setup_follow_3d(ax, ctr, half_um, az, el, cam_dist, proj_short, lw, fs_tick, fs_label)
%SETUP_FOLLOW_3D  Orthographic co-moving 3-D window around the box centre ctr
%  The axes box has the half-extents half_um = [hx hy hz] [um] around ctr (not
%  necessarily a cube; equal scale on all axes), with ticks, tick labels, grid
%  and axis labels (lab-frame coordinates).  The camera is placed manually and
%  zoomed so that the projected box fills the axes rectangle, whose aspect ratio
%  equals that of the projected box: MATLAB maps CameraViewAngle to the SHORTER
%  side of the axes, so the field of view is set to proj_short [um], the shorter
%  side of the projected box.
%  Outputs: unit vector target -> camera, camera up vector (as defined by
%  view(ax, az, el)), and the screen right / up unit vectors in data space.
    xlim(ax, ctr(1) + half_um(1) * [-1 1]);
    ylim(ax, ctr(2) + half_um(2) * [-1 1]);
    zlim(ax, ctr(3) + half_um(3) * [-1 1]);
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
    set(ax, 'CameraTarget', ctr, ...
            'CameraPosition', ctr + cam_dist * cam_dir, ...
            'CameraUpVector', cam_up, ...
            'CameraViewAngle', 2 * atand(proj_short / 2 / cam_dist));
    set_window_axes(ax, lw, fs_tick, fs_label, true);
    apply_window_ticks(ax, 0.5);
end

function setup_follow_xz(ax, ctr, half_um, lim_y, lw, fs_tick, fs_label)
%SETUP_FOLLOW_XZ  Orthographic co-moving XZ window (view along Y) around the box centre ctr
%  Axes box with half-extents half_um(1) (X) and half_um(3) (Z), equal scale on
%  both axes, with ticks, tick labels, grid and axis labels (lab-frame
%  coordinates).  The Y range stays fixed to the global range so that nothing is
%  clipped along the viewing direction.
    xlim(ax, ctr(1) + half_um(1) * [-1 1]);
    ylim(ax, lim_y);
    zlim(ax, ctr(3) + half_um(3) * [-1 1]);
    set(ax, 'DataAspectRatio', [1 1 1], 'Projection', 'orthographic');
    view(ax, 0, 0);
    set_window_axes(ax, lw, fs_tick, fs_label, false);
    apply_window_ticks(ax, 0.5);
end

function set_window_axes(ax, lw, fs_tick, fs_label, with_y_label)
%SET_WINDOW_AXES  Axes box, ticks (lab-frame coordinates), tick labels, grid, axis labels
%  Tick positions are set explicitly (see window_ticks) and refreshed whenever
%  the limits follow the centre of mass, so the grid slides past the bacterium.
%  Spacing is at least 0.5 μm.
    set(ax, 'Box', 'on', 'TickDir', 'out', 'LineWidth', lw, 'FontSize', fs_tick, ...
            'Layer', 'top');
    grid(ax, 'on');
    xlabel(ax, 'X (μm)', 'FontName', 'Arial', 'FontSize', fs_label);
    if with_y_label
        ylabel(ax, 'Y (μm)', 'FontName', 'Arial', 'FontSize', fs_label);
    end
    zlabel(ax, 'Z (μm)', 'FontName', 'Arial', 'FontSize', fs_label);
end

function apply_window_ticks(ax, min_step)
%APPLY_WINDOW_TICKS  Place X/Y/Z ticks on the current limits, spacing >= min_step
    set(ax, 'XTick', window_ticks(get(ax, 'XLim'), min_step), ...
            'YTick', window_ticks(get(ax, 'YLim'), min_step), ...
            'ZTick', window_ticks(get(ax, 'ZLim'), min_step));
end

function ticks = window_ticks(lim, min_step)
%WINDOW_TICKS  Lab-frame ticks inside lim with spacing at least min_step
%  The step is a 1-2-5 value.  Short axes, whose automatic step would be
%  0.2 μm, therefore use 0.5 μm; longer axes keep a coarser step.
    span = lim(2) - lim(1);
    if ~(isfinite(span) && span > 0)
        ticks = mean(lim);
        return;
    end
    step = max(min_step, nice_step_at_least(span / 6));
    tol = 1e-9 * max(1, abs(step));
    first = ceil(lim(1) / step - tol) * step;
    last  = floor(lim(2) / step + tol) * step;
    if last + tol < first
        ticks = round(mean(lim) / step) * step;
        return;
    end
    n = round((last - first) / step);
    ticks = first + (0:n) * step;
end

function step = nice_step_at_least(raw)
%NICE_STEP_AT_LEAST  Smallest {1, 2, 5} x 10^k that is at least raw
    if ~(isfinite(raw) && raw > 0)
        step = 1;
        return;
    end
    exp10 = floor(log10(raw));
    base = 10^exp10;
    mant = raw / base;
    if mant <= 1
        factor = 1;
    elseif mant <= 2
        factor = 2;
    elseif mant <= 5
        factor = 5;
    else
        factor = 10;
    end
    step = factor * base;
end

function [mark_flat, mark_size] = waist_mark_surface(pts, e_long, e_side)
%WAIST_MARK_SURFACE  Short patch on one side of the body waist
%  pts is [N x 3] and centred on the body.  The patch covers the central
%  18 % of the body length and ±32° of azimuth about e_side, and sits 2 %
%  outside the local radius.  mark_flat is [3 x N] and mark_size is the
%  surf grid.
    e_long = e_long(:) / norm(e_long);
    e_side = e_side(:) / norm(e_side);
    e_bin  = cross(e_long, e_side);
    s = pts * e_long;
    radial = pts - s * e_long.';
    rad = sqrt(sum(radial.^2, 2));
    s_half = max(abs(s));
    near_waist = abs(s) <= 0.35 * s_half;
    if ~any(near_waist)
        near_waist = true(size(s));
    end
    radius = median(rad(near_waist));
    half_len = 0.18 * s_half;
    azimuth = 32 * pi / 180;
    n_s = 16;
    n_a = 10;
    ss = linspace(-half_len, half_len, n_s);
    aa = linspace(-azimuth, azimuth, n_a);
    [SS, AA] = meshgrid(ss, aa);
    lifted = 1.02 * radius;
    P = e_long * SS(:).' + lifted * (e_side * cos(AA(:)).' + e_bin * sin(AA(:)).');
    mark_flat = P;
    mark_size = size(SS);
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
