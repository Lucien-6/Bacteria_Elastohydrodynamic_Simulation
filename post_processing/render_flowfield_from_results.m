function info = render_flowfield_from_results(mat_file, mode, requested_times_s, fps, skip, opts)
%RENDER_FLOWFIELD_FROM_RESULTS  Browse a saved run and draw its flow field
%
%  Description:
%    Independent entry for flow-field figures and animations.  When mat_file
%    is omitted, a file dialog is shown.  Frame rate and stride fall back to
%    config.json in the same folder, then to 30 fps and a stride of 1.
%    Grid settings fall back to that file's flowfield section.
%
%    Optional opts.drop_first_frame (default false): true omits saved output
%    index 1 from times / animation / all (drawing and colour scale), so a
%    startup transient does not dominate the figure.  final mode is unchanged.
%
%  Syntax:
%    info = render_flowfield_from_results()
%    info = render_flowfield_from_results(mat_file)
%    info = render_flowfield_from_results(mat_file, mode)
%    info = render_flowfield_from_results(mat_file, mode, requested_times_s)
%    info = render_flowfield_from_results(mat_file, mode, requested_times_s, fps)
%    info = render_flowfield_from_results(mat_file, mode, requested_times_s, fps, skip)
%    info = render_flowfield_from_results(mat_file, mode, requested_times_s, fps, skip, opts)
%
%  Outputs:
%    info - Summary structure, or [] when the file dialog is cancelled.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28
%  Updated: 2026-11-11 — opts.drop_first_frame
%  Updated: 2026-11-13 — optional stream_spacing_um and stream_length_um

info = [];
this_file = mfilename('fullpath');
proj_root = fileparts(fileparts(this_file));
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));
chlamy_dir = fullfile(proj_root, 'chlamydomonas');
if isfolder(chlamy_dir)
    addpath(chlamy_dir);
    chlamy_core = fullfile(chlamy_dir, 'core');
    if isfolder(chlamy_core)
        addpath(chlamy_core);
    end
end

if nargin < 1 || isempty(mat_file)
    [f_name, f_dir] = uigetfile({'*.mat', 'MAT-files (*.mat)'}, ...
        'Select simulation_results.mat');
    if isequal(f_name, 0)
        fprintf('No file selected - exiting.\n');
        return
    end
    mat_file = fullfile(f_dir, f_name);
end
mat_file = char(mat_file);
if ~isfile(mat_file)
    error('render_flowfield_from_results:FileNotFound', ...
        'Results file not found: %s', mat_file);
end
case_dir = fileparts(mat_file);
if isempty(case_dir)
    case_dir = pwd;
end

if nargin < 2 || isempty(mode), mode = 'final'; end
if nargin < 3, requested_times_s = []; end
if nargin < 4, fps = []; end
if nargin < 5, skip = []; end
if nargin < 6 || isempty(opts), opts = struct(); end

S = load_saved_results(mat_file);
if ~isfield(S, 'results') || ~isfield(S, 'params')
    error('render_flowfield_from_results:MissingVariables', ...
        'The file must contain the variables "results" and "params": %s', mat_file);
end

[cfg_fps, cfg_skip, cfg_flow] = read_case_config(fullfile(case_dir, 'config.json'));
if isempty(fps), fps = cfg_fps; end
if isempty(skip), skip = cfg_skip; end
opts = apply_flow_defaults(opts, cfg_flow);
opts.mode = mode;
opts.requested_times_s = requested_times_s;
opts.fps = fps;
opts.skip = skip;
opts.input_file = mat_file;
if ~isfield(opts, 'run_id') || isempty(opts.run_id)
    [~, opts.run_id] = fileparts(case_dir);
end
if ~isfield(opts, 'output_dir') || isempty(opts.output_dir)
    out_dir = case_dir;
else
    out_dir = opts.output_dir;
end
if ~isfield(opts, 'silent') || isempty(opts.silent)
    opts.silent = false;
end
if ~isfield(opts, 'drop_first_frame') || isempty(opts.drop_first_frame)
    opts.drop_first_frame = false;
end

fprintf('Results file : %s\n', mat_file);
fprintf('Mode         : %s\n', char(mode));
fprintf('Drop first   : %d\n', logical(opts.drop_first_frame));
info = call_saved_flowfield(S.results, S.params, out_dir, opts);
end


function S = load_saved_results(mat_file)
%LOAD_SAVED_RESULTS  Load results & params without benign legacy-MAT warnings.
%  Historical .mat files may contain a non-serialisable body_precomp.dM_BB
%  decomposition object (rebuilt automatically by ensure_body_precomp) or
%  anonymous active-moment handles recorded before folder reorganisation
%  (repaired automatically by restore_saved_moment_handles).
w_decomp = warning('query', 'MATLAB:decomposition:LoadNotSupported');
w_handle = warning('query', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
restore_w = onCleanup(@() restore_warnings(w_decomp, w_handle));
warning('off', 'MATLAB:decomposition:LoadNotSupported');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
S = load(mat_file, 'results', 'params');
end


function restore_warnings(w_decomp, w_handle)
warning(w_decomp.state, 'MATLAB:decomposition:LoadNotSupported');
warning(w_handle.state, 'MATLAB:dispatcher:UnresolvedFunctionHandle');
end


function [fps, skip, flow] = read_case_config(cfg_json)
fps = 30;
skip = 1;
flow = struct();
if ~isfile(cfg_json)
    return
end
try
    jc = jsondecode(fileread(cfg_json));
catch ME
    warning('render_flowfield_from_results:BadConfigJson', ...
        'Could not read %s (%s). Using 30 fps and a stride of 1.', ...
        cfg_json, ME.message);
    return
end
if isfield(jc, 'output')
    if isfield(jc.output, 'animation_fps') && ~isempty(jc.output.animation_fps)
        fps = jc.output.animation_fps;
    end
    if isfield(jc.output, 'animation_skip') && ~isempty(jc.output.animation_skip)
        skip = jc.output.animation_skip;
    end
end
if isfield(jc, 'flowfield') && isstruct(jc.flowfield)
    flow = jc.flowfield;
end
end


function opts = apply_flow_defaults(opts, flow)
pairs = {'grid_spacing_um', 'grid_spacing_um'; ...
         'grid_range_um', 'grid_range_um'; ...
         'quiver_skip', 'quiver_skip'; ...
         'chunk_size', 'chunk_size'; ...
         'stream_spacing_um', 'stream_spacing_um'; ...
         'stream_length_um', 'stream_length_um'; ...
         'use_gpu', 'use_gpu'; ...
         'use_parallel', 'use_parallel'; ...
         'scale_from_video_frames', 'scale_from_video_frames'};
for n = 1:size(pairs, 1)
    target = pairs{n, 1};
    source = pairs{n, 2};
    if (~isfield(opts, target) || isempty(opts.(target))) && isfield(flow, source)
        opts.(target) = flow.(source);
    end
end
end
