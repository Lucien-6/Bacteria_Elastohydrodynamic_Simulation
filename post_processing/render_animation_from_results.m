function info = render_animation_from_results(mat_file, fps, skip)
%RENDER_ANIMATION_FROM_RESULTS  Re-render the motion animation from saved results
%
%  Description:
%    Re-creates the 4-panel motion animation (3D and XZ overviews plus
%    co-moving 3D and XZ windows) of a finished simulation from its
%    simulation_results.mat, without re-running the ODE
%    integration.  Useful for changing the frame stride or frame rate, or for
%    regenerating a lost / corrupted video.
%
%    The rendering itself is done by src/visualization/render_motion_animation.m,
%    the same function used by run_simulation.m and run_simulation_HPC_Linux.m.
%    The output is written next to the .mat file as
%    bacteria_3d_motion_rerender.mp4 (Windows / macOS) or
%    bacteria_3d_motion_rerender.avi (Linux), so the original video is never
%    overwritten.
%
%  Syntax:
%    info = render_animation_from_results()
%    info = render_animation_from_results(mat_file)
%    info = render_animation_from_results(mat_file, fps)
%    info = render_animation_from_results(mat_file, fps, skip)
%
%  Inputs:
%    mat_file - (Optional) Path of simulation_results.mat.  When omitted or
%               empty, a file-selection dialog is shown.
%    fps      - (Optional) Video frame rate [frames/s].
%    skip     - (Optional) Frame stride over the saved time points.
%               When fps or skip is omitted, the values are read from
%               config.json (output.animation_fps / output.animation_skip) in
%               the same folder; if that file is missing, fps = 30 and
%               skip = 1.
%
%  Outputs:
%    info - Structure returned by render_motion_animation (fields: file,
%           filename, n_frames, elapsed_s, stage_s, window_half_um,
%           window_offset_um, playback_speed); empty when the file selection is cancelled.
%
%  Platform and requirements:
%    MATLAB R2018b+ (min/max over multiple dimensions); the MPEG-4 profile
%    requires Windows or macOS, Linux falls back to Motion JPEG AVI.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-17
%  Updated: 2026-10-18 — header follows the 4-panel layout of the renderer
%  Updated: 2026-10-19 — info gains window_offset_um (v8.6.3)

info = [];

% Standalone entry: add project root and src (same policy as postprocess_cases.m)
this_file = mfilename('fullpath');
proj_root = fileparts(fileparts(this_file));
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));

%% ---- Select the results file ----
if nargin < 1 || isempty(mat_file)
    [f_name, f_dir] = uigetfile({'*.mat', 'MAT-files (*.mat)'}, ...
        'Select simulation_results.mat');
    if isequal(f_name, 0)
        fprintf('No file selected - exiting.\n');
        return;
    end
    mat_file = fullfile(f_dir, f_name);
end
mat_file = char(mat_file);
if ~isfile(mat_file)
    error('render_animation_from_results:FileNotFound', ...
          'Results file not found: %s', mat_file);
end
case_dir = fileparts(mat_file);
if isempty(case_dir)
    case_dir = pwd;
end

%% ---- Load and validate ----
S = load(mat_file, 'results', 'params');
if ~isfield(S, 'results') || ~isfield(S, 'params')
    error('render_animation_from_results:MissingVariables', ...
          'The file must contain the variables "results" and "params": %s', mat_file);
end
require_fields(S.results, {'t', 'q', 'X_cm_um', 'X_pili_um'}, 'results');
require_fields(S.params,  {'M_pili', 'geometry_mode', 'body_type', 'dim'}, 'params');
require_fields(S.params.dim, {'major_axis', 'minor_axis'}, 'params.dim');
if ~(strcmp(S.params.geometry_mode, 'parametric') && strcmp(S.params.body_type, 'capsule'))
    require_fields(S.params, {'body_verts', 'body_faces'}, 'params');
end

%% ---- Frame rate and stride ----
cfg_fps  = 30;
cfg_skip = 1;
cfg_json = fullfile(case_dir, 'config.json');
if (nargin < 2 || isempty(fps) || nargin < 3 || isempty(skip)) && isfile(cfg_json)
    try
        jc = jsondecode(fileread(cfg_json));
        if isfield(jc, 'output')
            if isfield(jc.output, 'animation_fps') && ~isempty(jc.output.animation_fps)
                cfg_fps = jc.output.animation_fps;
            end
            if isfield(jc.output, 'animation_skip') && ~isempty(jc.output.animation_skip)
                cfg_skip = jc.output.animation_skip;
            end
        end
    catch
        % Unreadable config.json: keep the defaults.
    end
end
if nargin < 2 || isempty(fps),  fps  = cfg_fps;  end
if nargin < 3 || isempty(skip), skip = cfg_skip; end
if ~(isnumeric(fps) && isscalar(fps) && isfinite(fps) && fps > 0)
    error('render_animation_from_results:InvalidFps', 'fps must be a positive scalar.');
end
if ~(isnumeric(skip) && isscalar(skip) && isfinite(skip) && skip >= 1 && skip == round(skip))
    error('render_animation_from_results:InvalidSkip', 'skip must be a positive integer.');
end

%% ---- Render ----
fprintf('Results file : %s\n', mat_file);
fprintf('Frame rate   : %g fps | frame stride: %d\n', fps, skip);

opts = struct('fps', fps, 'skip', skip, 'silent', false, ...
              'file_basename', 'bacteria_3d_motion_rerender');
info = render_motion_animation(S.results, S.params, case_dir, opts);

fprintf('Done: %s (%d frames, %.1f s)\n', info.file, info.n_frames, info.elapsed_s);

end


function require_fields(S, names, label)
%REQUIRE_FIELDS  Error out if any of the named fields is missing from struct S
    missing = names(~isfield(S, names));
    if ~isempty(missing)
        error('render_animation_from_results:MissingFields', ...
              'Missing field(s) in %s: %s', label, strjoin(missing, ', '));
    end
end
