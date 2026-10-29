function info = render_flowfield_from_results(varargin)
%RENDER_FLOWFIELD_FROM_RESULTS  Flow fields from a saved simulation
%
%  Description:
%    Draws final, selected-time, and full-sequence flow fields from a
%    version-1 results structure.  This is the computational entry used by
%    post_processing/render_flowfield_from_results.m.  Calling it with no
%    arguments, or with a results-file path, opens that same file browser.
%
%    Colour scales:
%      final      - 99th percentile of the final frame
%      times      - 99th percentile of the requested frames only
%      animation  - 99th percentile of every saved output
%    Every saved output, including the last, uses its recorded
%    closing-interval Brownian force.  Velocity grids are not retained
%    between the sampling pass and the drawing pass.
%
%  Syntax:
%    info = render_flowfield_from_results(results, params, output_dir, opts)
%    info = render_flowfield_from_results()
%    info = render_flowfield_from_results(mat_file, mode, requested_times_s, fps, skip, opts)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28

if nargin == 0 || ischar(varargin{1}) || isstring(varargin{1})
    info = open_saved_flowfield_ui(varargin{:});
    return
end
if nargin ~= 4 || ~isstruct(varargin{1}) || ~isstruct(varargin{2}) || ~isstruct(varargin{4})
    error('render_flowfield_from_results:MissingOptions', ...
        'Pass loaded results, params, an output directory, and an options struct.');
end
results = varargin{1};
params = varargin{2};
output_dir = varargin{3};
opts = varargin{4};
if ~isfield(opts, 'mode') || isempty(opts.mode)
    error('render_flowfield_from_results:MissingMode', 'opts.mode is required.');
end

logf = opt_of(opts, 'logf', @(varargin) fprintf(varargin{:}));
silent = logical(opt_of(opts, 'silent', false));
if ~isfolder(output_dir)
    mkdir(output_dir);
end

mode = lower(char(opts.mode));
flow_opts = flow_options(opts);
t_start = tic;
info = blank_info(mode, opt_of(opts, 'input_file', ''));

switch mode
    case 'final'
        info = render_final(info, results, params, output_dir, opts, flow_opts, logf, silent);
    case 'times'
        info = render_times(info, results, params, output_dir, opts, flow_opts, logf, silent);
    case 'animation'
        info = render_animation(info, results, params, output_dir, opts, flow_opts, logf, silent);
    case 'all'
        info = render_all(info, results, params, output_dir, opts, flow_opts, logf, silent);
    otherwise
        error('render_flowfield_from_results:BadMode', ...
            'mode must be final, times, animation, or all.');
end
info.elapsed_s = toc(t_start);
logf('[FLOW] Post-processing completed in %.1f s\n', info.elapsed_s);
end


function info = render_final(info, results, params, output_dir, opts, flow_opts, logf, silent)
n_times = results.flow_state.n_times;
logf('[FLOW] Final frame: one evaluation of output %d.\n', n_times);
held = hold_frames(results, params, n_times, n_times, flow_opts, logf, silent, opts);
cmax = colour_limit(samples_for(held, n_times));
png = fullfile(output_dir, 'flowfield_panels_rerender_final.png');
[info, held] = draw_static_held(info, held, n_times, cmax, ...
    'Frame colour scale', {png}, flow_opts, logf, silent, opts);
info.final_output_index = n_times;
info.output_index = [info.output_index, n_times];
info.matched_times_s = [info.matched_times_s, results.t(n_times)];
info.convention{end + 1} = 'closing_interval';
info.colour_max.final = cmax;
info.grid_recompute_count = held.grid_recompute_count;
logf('[FLOW] Grid cache %.1f MiB | grid recomputes = %d\n', ...
    held.retained_bytes / 1024^2, held.grid_recompute_count);
end


function info = render_times(info, results, params, output_dir, opts, flow_opts, logf, silent)
if ~isfield(opts, 'requested_times_s') || isempty(opts.requested_times_s)
    error('render_flowfield_from_results:MissingTimes', ...
        'requested_times_s is required for times output.');
end
[idx_unique, requested_index, requested_times] = match_output_times( ...
    results.t, opts.requested_times_s);
logf('[FLOW] Requested-time colour scale uses %d matched output(s).\n', numel(idx_unique));
held = hold_frames(results, params, idx_unique, idx_unique, flow_opts, logf, silent, opts);
cmax = colour_limit(samples_for(held, idx_unique));
names = cell(1, numel(idx_unique));
for n = 1:numel(idx_unique)
    names{n} = fullfile(output_dir, sprintf('flowfield_panels_t%06d.png', n));
end
[info, held] = draw_static_held(info, held, idx_unique, cmax, ...
    'Requested-time colour scale', names, flow_opts, logf, silent, opts);
info.requested_times_s = reshape(requested_times, 1, []);
info.requested_output_index = reshape(requested_index, 1, []);
info.times_output_index = reshape(idx_unique, 1, []);
info.output_index = [info.output_index, reshape(idx_unique, 1, [])];
info.matched_times_s = [info.matched_times_s, reshape(results.t(idx_unique), 1, [])];
for n = 1:numel(idx_unique)
    info.convention{end + 1} = convention_for(idx_unique(n), results.flow_state.n_times);
end
info.colour_max.times = cmax;
info.grid_recompute_count = held.grid_recompute_count;
logf('[FLOW] Grid cache %.1f MiB | grid recomputes = %d\n', ...
    held.retained_bytes / 1024^2, held.grid_recompute_count);
end


function info = render_animation(info, results, params, output_dir, opts, flow_opts, logf, silent)
[fps, skip, frame_index] = animation_frame_index(results, opts);
n_times = results.flow_state.n_times;
logf('[FLOW] Animation colour scale uses all %d saved outputs.\n', n_times);
held = hold_frames(results, params, 1:n_times, frame_index, flow_opts, logf, silent, opts);
cmax = colour_limit(samples_for(held, 1:n_times));
[info, held] = write_animation_held(info, held, frame_index, cmax, fps, skip, ...
    results, output_dir, opts, flow_opts, logf, silent);
info.colour_max.animation = cmax;
info.grid_recompute_count = held.grid_recompute_count;
logf('[FLOW] Grid cache %.1f MiB | grid recomputes = %d\n', ...
    held.retained_bytes / 1024^2, held.grid_recompute_count);
end


function info = render_all(info, results, params, output_dir, opts, flow_opts, logf, silent)
if ~isfield(opts, 'requested_times_s') || isempty(opts.requested_times_s)
    error('render_flowfield_from_results:MissingTimes', ...
        'requested_times_s is required for times output.');
end
[fps, skip, frame_index] = animation_frame_index(results, opts);
[idx_unique, requested_index, requested_times] = match_output_times( ...
    results.t, opts.requested_times_s);
n_times = results.flow_state.n_times;
retain = unique([n_times, reshape(idx_unique, 1, []), reshape(frame_index, 1, [])]);
logf('[FLOW] Combined output evaluates each saved frame once.\n');
held = hold_frames(results, params, 1:n_times, retain, flow_opts, logf, silent, opts);

cmax_final = colour_limit(samples_for(held, n_times));
png = fullfile(output_dir, 'flowfield_panels_rerender_final.png');
[info, held] = draw_static_held(info, held, n_times, cmax_final, ...
    'Frame colour scale', {png}, flow_opts, logf, silent, opts);
info.final_output_index = n_times;
info.output_index = [info.output_index, n_times];
info.matched_times_s = [info.matched_times_s, results.t(n_times)];
info.convention{end + 1} = 'closing_interval';
info.colour_max.final = cmax_final;

cmax_times = colour_limit(samples_for(held, idx_unique));
names = cell(1, numel(idx_unique));
for n = 1:numel(idx_unique)
    names{n} = fullfile(output_dir, sprintf('flowfield_panels_t%06d.png', n));
end
[info, held] = draw_static_held(info, held, idx_unique, cmax_times, ...
    'Requested-time colour scale', names, flow_opts, logf, silent, opts);
info.requested_times_s = reshape(requested_times, 1, []);
info.requested_output_index = reshape(requested_index, 1, []);
info.times_output_index = reshape(idx_unique, 1, []);
info.output_index = [info.output_index, reshape(idx_unique, 1, [])];
info.matched_times_s = [info.matched_times_s, reshape(results.t(idx_unique), 1, [])];
for n = 1:numel(idx_unique)
    info.convention{end + 1} = convention_for(idx_unique(n), n_times);
end
info.colour_max.times = cmax_times;

cmax_anim = colour_limit(samples_for(held, 1:n_times));
[info, held] = write_animation_held(info, held, frame_index, cmax_anim, fps, skip, ...
    results, output_dir, opts, flow_opts, logf, silent);
info.colour_max.animation = cmax_anim;
info.grid_recompute_count = held.grid_recompute_count;
logf('[FLOW] Grid cache %.1f MiB | grid recomputes = %d\n', ...
    held.retained_bytes / 1024^2, held.grid_recompute_count);
end


function [fps, skip, frame_index] = animation_frame_index(results, opts)
fps = opt_of(opts, 'fps', []);
skip = opt_of(opts, 'skip', []);
if ~(isnumeric(fps) && isscalar(fps) && isfinite(fps) && fps > 0)
    error('render_flowfield_from_results:InvalidFps', 'fps must be a positive scalar.');
end
if ~(isnumeric(skip) && isscalar(skip) && isfinite(skip) && skip >= 1 && skip == round(skip))
    error('render_flowfield_from_results:InvalidSkip', 'skip must be a positive integer.');
end
n_times = results.flow_state.n_times;
frame_index = 1:skip:n_times;
if frame_index(end) ~= n_times
    frame_index(end + 1) = n_times;
end
end


function held = hold_frames(results, params, visit, retain, flow_opts, logf, silent, opts)
%HOLD_FRAMES  Evaluate each visited output once and retain render grids.
if nargin < 8
    opts = struct();
end
cap = 512 * 1024^2;
if isfield(opts, 'grid_cache_bytes') && ~isempty(opts.grid_cache_bytes)
    cap = opts.grid_cache_bytes;
end
n_times = results.flow_state.n_times;
held.pack = cell(n_times, 1);
held.sample = cell(n_times, 1);
held.retained_bytes = 0;
held.grid_recompute_count = 0;
retain = unique(reshape(retain, 1, []));
visit = reshape(visit, 1, []);
for n = 1:numel(visit)
    idx = visit(n);
    packed = evaluate_saved_index(results, params, idx, flow_opts, logf, silent);
    held.sample{idx} = packed.planes.samples;
    if any(retain == idx)
        nbytes = struct_bytes(packed.planes);
        if held.retained_bytes + nbytes <= cap
            held.pack{idx} = packed;
            held.retained_bytes = held.retained_bytes + nbytes;
        else
            packed.planes = [];
            held.pack{idx} = packed;
            logf(['[FLOW] Grid cache full; output %d will rebuild its grid ', ...
                  'from saved coefficients.\n'], idx);
        end
    end
    logf('[FLOW] Evaluated output %d at t = %.6g s (%s).\n', ...
        idx, packed.t_s, packed.convention);
end
end


function samples = samples_for(held, indices)
parts = held.sample(reshape(indices, 1, []));
samples = vertcat(parts{:});
end


function [info, held] = draw_static_held(info, held, indices, cmax, scale_note, names, flow_opts, logf, silent, opts)
visible = 'on';
if silent
    visible = 'off';
end
run_id = opt_of(opts, 'run_id', 'flowfield');
dpi = opt_of(opts, 'dpi', 600);
indices = reshape(indices, 1, []);
for n = 1:numel(indices)
    [packed, held] = materialize_pack(held, indices(n), flow_opts, logf, silent);
    frame = pack_frame(packed, run_id, cmax, scale_note);
    render_flowfield_panels(frame, struct('action', 'static', ...
        'visible', visible, 'bundle_export', false, 'close_after', true, ...
        'png_path', names{n}, 'dpi', dpi));
    info.static_files{end + 1} = names{n};
    logf('[FLOW] Static panel saved: %s (t = %.6g s, %s)\n', ...
        names{n}, packed.t_s, packed.convention);
end
end


function [info, held] = write_animation_held(info, held, frame_index, cmax, fps, skip, results, output_dir, opts, flow_opts, logf, silent)
n_frames = numel(frame_index);
n_times = results.flow_state.n_times;
if n_frames > 1
    playback = (results.t(frame_index(end)) - results.t(frame_index(1))) ...
        / (n_frames - 1) * fps;
else
    playback = NaN;
end
[vid, anim_file, cleanup_vid] = open_flow_video(output_dir, opts, fps);
handles = [];
write_png = logical(opt_of(opts, 'write_frame_png', false));
dpi = opt_of(opts, 'dpi', 600);
frame_dir = fullfile(output_dir, 'flowfield_frames');
if write_png && ~isfolder(frame_dir)
    mkdir(frame_dir);
end
vid_h = [];
vid_w = [];
run_id = opt_of(opts, 'run_id', 'flowfield');
if write_png
    frame_files = cell(1, n_frames);
else
    frame_files = {};
end
for n = 1:n_frames
    idx = frame_index(n);
    [packed, held] = materialize_pack(held, idx, flow_opts, logf, silent);
    frame = pack_frame(packed, run_id, cmax, 'Global colour scale');
    if n == 1
        handles = render_flowfield_panels(frame, ...
            struct('action', 'create', 'visible', 'off'));
        cleanup_fig = onCleanup(@() close_if_open(handles.fig));
    else
        handles = render_flowfield_panels(frame, ...
            struct('action', 'update', 'handles', handles));
    end
    shot = getframe(handles.fig);
    [shot.cdata, vid_h, vid_w] = fit_frame(shot.cdata, vid_h, vid_w);
    writeVideo(vid, shot.cdata);
    if write_png
        png = fullfile(frame_dir, sprintf('flowfield_frame_%06d.png', n));
        render_flowfield_panels(frame, struct('action', 'export_png', ...
            'figure', handles.fig, 'png_path', png, 'dpi', dpi));
        frame_files{n} = png;
    end
    logf('[FLOW] Animation frame %d/%d at t = %.6g s (%s).\n', ...
        n, n_frames, packed.t_s, packed.convention);
end
cleanup_vid();
if exist('cleanup_fig', 'var')
    cleanup_fig();
end
info.animation_file = anim_file;
info.frame_files = frame_files;
info.fps = fps;
info.skip = skip;
info.n_frames = n_frames;
info.playback_speed = playback;
info.animation_output_index = reshape(frame_index, 1, []);
info.output_index = [info.output_index, reshape(frame_index, 1, [])];
info.matched_times_s = [info.matched_times_s, reshape(results.t(frame_index), 1, [])];
for n = 1:n_frames
    info.convention{end + 1} = convention_for(frame_index(n), n_times);
end
logf('[FLOW] Animation saved: %s (%d frames)\n', anim_file, n_frames);
end


function [packed, held] = materialize_pack(held, index, flow_opts, logf, silent)
packed = held.pack{index};
if isempty(packed.planes)
    packed.planes = flowfield_grid_velocity(packed.hydro, packed.params, ...
        flow_opts, logf, silent);
    held.pack{index} = packed;
    held.grid_recompute_count = held.grid_recompute_count + 1;
    logf('[FLOW] Rebuilt velocity grid for output %d.\n', index);
end
end


function n = struct_bytes(value)
n = 0;
if isnumeric(value)
    n = numel(value) * 8;
elseif islogical(value)
    n = numel(value);
elseif isstruct(value)
    names = fieldnames(value);
    for i = 1:numel(names)
        n = n + struct_bytes(value.(names{i}));
    end
elseif iscell(value)
    for i = 1:numel(value)
        n = n + struct_bytes(value{i});
    end
end
end

function packed = evaluate_saved_index(results, params, index, flow_opts, logf, silent)
spec = flowfield_state_data(results, params, index, ...
    kind_for(index, results.flow_state.n_times));
try
    hydro = flowfield_from_saved_state(spec.Z, spec.params, spec.t_dimless, logf);
    planes = flowfield_grid_velocity(hydro, spec.params, flow_opts, logf, silent);
catch ME
    err = MException('render_flowfield_from_results:FrameFailed', ...
        'Flow-field evaluation failed at output %d (t = %.6g s).', ...
        index, spec.t_s);
    err = addCause(err, ME);
    throw(err);
end
packed.planes = planes;
packed.hydro = hydro;
packed.geom = hydro.geom;
packed.params = spec.params;
packed.L_to_um = hydro.L_to_um;
packed.t_s = spec.t_s;
packed.convention = spec.convention;
packed.index = index;
end


function frame = pack_frame(packed, run_id, cmax, scale_note)
frame.planes = packed.planes;
frame.geom = packed.geom;
frame.params = packed.params;
frame.L_to_um = packed.L_to_um;
frame.cmax = cmax;
frame.time_s = packed.t_s;
frame.run_id = run_id;
frame.scale_note = scale_note;
end


function [idx_unique, requested_index, requested_times] = match_output_times(t_saved, t_req)
requested_times = double(t_req(:));
if isempty(requested_times) || any(~isfinite(requested_times)) || any(requested_times < 0)
    error('render_flowfield_from_results:BadRequestedTime', ...
        'Requested times must be finite and nonnegative, in seconds.');
end
t_saved = double(t_saved(:));
if any(requested_times < t_saved(1) | requested_times > t_saved(end))
    error('render_flowfield_from_results:TimeOutOfRange', ...
        'Requested times must lie within the saved interval [%.6g, %.6g] s.', ...
        t_saved(1), t_saved(end));
end
requested_index = zeros(numel(requested_times), 1);
for n = 1:numel(requested_times)
    [~, requested_index(n)] = min(abs(t_saved - requested_times(n)));
end
idx_unique = unique(requested_index, 'stable').';
end


function kind = kind_for(index, n_times)
% Every output uses the single recorded interval force. n_times is accepted
% so callers can keep the output index visible at the call site.
if index < 1 || index > n_times
    error('render_flowfield_from_results:BadIndex', ...
        'Output index %d is outside 1:%d.', index, n_times);
end
kind = 'history';
end


function convention = convention_for(index, n_times)
kind_for(index, n_times);
convention = 'closing_interval';
end


function cmax = colour_limit(samples)
samples = samples(isfinite(samples));
if isempty(samples)
    cmax = 1;
    return
end
cmax = prctile(samples, 99);
if ~(isfinite(cmax) && cmax > 0)
    cmax = max(samples);
end
if ~(isfinite(cmax) && cmax > 0)
    cmax = 1;
end
end


function flow_opts = flow_options(opts)
flow_opts = struct();
names = {'grid_spacing_um', 'grid_range_um', 'quiver_skip', 'chunk_size'};
for n = 1:numel(names)
    if isfield(opts, names{n}) && ~isempty(opts.(names{n}))
        flow_opts.(names{n}) = opts.(names{n});
    end
end
if isfield(flow_opts, 'grid_spacing_um') && ~(flow_opts.grid_spacing_um > 0)
    error('render_flowfield_from_results:BadGrid', 'grid_spacing_um must be positive.');
end
if isfield(flow_opts, 'grid_range_um') && ~(flow_opts.grid_range_um > 0)
    error('render_flowfield_from_results:BadGrid', 'grid_range_um must be positive.');
end
if isfield(flow_opts, 'quiver_skip') && ~(flow_opts.quiver_skip >= 1 && ...
        flow_opts.quiver_skip == round(flow_opts.quiver_skip))
    error('render_flowfield_from_results:BadGrid', 'quiver_skip must be a positive integer.');
end
if isfield(flow_opts, 'chunk_size') && ~(flow_opts.chunk_size >= 1 && ...
        flow_opts.chunk_size == round(flow_opts.chunk_size))
    error('render_flowfield_from_results:BadGrid', 'chunk_size must be a positive integer.');
end
end


function [vid, anim_file, cleanup_vid] = open_flow_video(output_dir, opts, fps)
basename = opt_of(opts, 'file_basename', 'flowfield_animation');
if ispc || ismac
    anim_file = fullfile(output_dir, [basename, '.mp4']);
    profile = 'MPEG-4';
else
    anim_file = fullfile(output_dir, [basename, '.avi']);
    profile = 'Motion JPEG AVI';
end
vid = VideoWriter(anim_file, profile);
vid.FrameRate = fps;
vid.Quality = opt_of(opts, 'animation_quality', 100);
open(vid);
cleanup_vid = onCleanup(@() close_video(vid));
end


function [cdata, vid_h, vid_w] = fit_frame(cdata, vid_h, vid_w)
if isempty(vid_h)
    vid_h = max(2, 2 * floor(size(cdata, 1) / 2));
    vid_w = max(2, 2 * floor(size(cdata, 2) / 2));
end
if size(cdata, 1) == vid_h && size(cdata, 2) == vid_w
    cdata = cdata(:, :, 1:min(3, size(cdata, 3)));
    return
end
if exist('imresize', 'file') == 2
    cdata = imresize(cdata(:, :, 1:min(3, size(cdata, 3))), [vid_h, vid_w]);
else
    error('render_flowfield_from_results:FrameSize', ...
        'Animation frames changed size and imresize is unavailable.');
end
end


function info = blank_info(mode, input_file)
info.input_file = input_file;
info.mode = mode;
info.requested_times_s = [];
info.matched_times_s = [];
info.requested_output_index = [];
info.output_index = [];
info.convention = {};
info.static_files = {};
info.animation_file = '';
info.frame_files = {};
info.fps = [];
info.skip = [];
info.n_frames = 0;
info.playback_speed = NaN;
info.colour_max = struct();
info.elapsed_s = 0;
end


function value = opt_of(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end


function close_video(vid)
try
    close(vid);
catch
end
end


function close_if_open(fig)
if ~isempty(fig) && isgraphics(fig)
    close(fig);
end
end
