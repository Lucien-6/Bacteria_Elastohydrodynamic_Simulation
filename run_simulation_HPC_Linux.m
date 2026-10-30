function [results, params, output_dir] = run_simulation_HPC_Linux(cfg_or_path)
%RUN_SIMULATION_HPC_LINUX  Execute the HPC/Linux coupled body-pili simulation.
%
%  Description:
%    Maintained, headless-cluster-oriented fork of run_simulation.m.  This
%    file is NOT auto-synchronised: physics, algorithm, and logging fixes in
%    run_simulation.m must be ported here manually.  Intentional differences:
%      1. VideoWriter uses MPEG-4 (.mp4) on Windows/macOS and Motion JPEG AVI
%         (.avi) on Linux, where MATLAB does not support the MPEG-4 profile.
%         This choice is implemented in the shared render_motion_animation.m.
%      2. Existing-output detection recognises .avi products.
%      3. Optional progress_queue/progress_tag fields provide throttled
%         (at least 30 s apart) ODE heartbeat messages for silent batch jobs.
%    Other simulation, post-processing, logging, and flow-field behaviour is
%    intended to match run_simulation.m.
%
%  Syntax:
%    [results, params, output_dir] = run_simulation_HPC_Linux(cfg)
%    [results, params, output_dir] = run_simulation_HPC_Linux(json_path)
%
%  Inputs:
%    cfg_or_path - Configuration struct (from load_sim_config or inline),
%                  OR a char/string path to a JSON config file.
%                  Important optional fields:
%                    cfg.compute_flowfield - logical; default true if omitted.
%                    cfg.brown_rng_seed - positive numeric scalar; rounded and
%                      used with a run-local mt19937ar Brownian stream.  If
%                      omitted, an MMddHHmmss seed is generated and saved.
%                    cfg.silent - runtime-only logical (default false);
%                      suppresses console/progress/warnings but not file logs
%                      or errors.  silent + preview_initial requires
%                      auto_confirm = true; otherwise the function throws
%                      run_simulation_HPC_Linux:SilentPreviewDeadlock.
%                    cfg.jacobian_reuse - runtime-only logical (default true).
%                      With Brownian motion enabled, ode15s re-uses the last
%                      exact finite-difference Jacobian across the dt_B
%                      restarts (LaggedJacobian.m; v8.5.1).  Only the Newton
%                      iteration matrix is affected, not the discrete
%                      equations or the error control.  false reproduces
%                      v8.5.0 bit-for-bit.  Ignored for deterministic runs.
%                    cfg.jacobian_max_age - runtime-only positive number of
%                      Brownian nodes (default 50); a fresh exact Jacobian is
%                      computed when the cached one is this many nodes old.
%                    cfg.hydro_radius_factor - positive number (default sqrt(2);
%                      also JSON material.hydro_radius_factor).  Pilus
%                      hydrodynamic regularisation length eps_pili = factor *
%                      r_pili (v8.5.2).  sqrt(2) gives the exact no-slip
%                      cylinder twist drag and translation within ~3 %;
%                      sqrt(e) gives exact translation; 1 reproduces <= 8.5.1.
%                    cfg.progress_queue - optional parallel.pool.DataQueue
%                      receiving heartbeat structs; empty means no heartbeat.
%                    cfg.progress_tag - caller-defined value copied into each
%                      heartbeat as its tag field; empty by default.
%
%  Errors raised while running are captured by a top-level try/catch and
%  logged to simulation.log as a [FATAL] block (identifier, message, full
%  stack trace) before being rethrown, so simulation.log always records
%  why a run failed instead of truncating silently.
%
%  Outputs:
%    results    - Results struct.  A segment timeout returns valid partial
%                 results rather than []; principal fields include
%                 t [1 x NT, s], X_cm_um [3 x NT, um],
%                 X_pili_um {M}[N_k+1 x 3 x NT, um],
%                 V_sink [signed dZ/dt, um/s], V_cm_um [3 x NT, um/s],
%                 tilt_angle [rad], F_sink_pN [pN], units, and the
%                 corresponding dimensionless state histories, completed,
%                 termination_reason, target_time_s, and achieved_time_s.
%                 jacobian_reuse is a struct of lagged-Jacobian counters
%                 (enabled, n_fresh, n_lagged, n_on_demand, n_age_refresh,
%                 n_invalidated, n_fallback, n_fd_rhs); only .enabled = false
%                 when the cache is not used.
%    params     - Model parameters and dimensional scale factors; [] if aborted.
%    output_dir - Created run directory; '' if aborted before allocation.
%
%  Side effects and files:
%    Adds the project root and project-root/src recursively to the MATLAB path.
%    Writes the same MAT/CSV/JSON/log/plot products as run_simulation.m.
%    Animation output is bacteria_3d_motion.mp4 on Windows/macOS and
%    bacteria_3d_motion.avi on Linux.  Brownian runs use an isolated stream
%    and save its seed without modifying the global RNG.  Since v8.5.0 the
%    body-only Brownian increment of each interval dt_B is realised as a
%    constant thermal generalised force acting through the coupled ODE
%    (sample_brownian_body_force) rather than as a displacement kick.
%    A configured progress_queue receives
%    structs with fields type='heartbeat', tag, pct, and elapsed.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17
%  Updated: 2026-07-21
%  Updated: 2026-09-04 — v8.5.0 Brownian force-realisation scheme
%  Updated: 2026-10-10 — v8.5.1 lagged Jacobian for Brownian ode15s restarts
%           and static-analysis clean-up (no behavioural change)
%  Updated: 2026-10-11 — v8.5.2 pilus hydrodynamic radius factor (hydro_radius_factor)
%  Updated: 2026-10-17 — v8.6.0 animation rendering moved to the shared
%           src/visualization/render_motion_animation.m (graphics objects created once,
%           only data updated per frame; video content unchanged)
%  Updated: 2026-10-18 — v8.6.2 animation layout: 4 panels (overviews + co-moving
%           windows), legend, scale bar, time stamp with playback speed
%  Updated: 2026-10-19 — v8.6.3 animation: windows sized per axis from the reach of the
%           bacterium, roll-marker stripes on the body, adaptive panel layout
%  Updated: 2026-10-20 — v8.6.5 realistic single-flagellum E. coli example case
%           (no solver or animation change)
%  Updated: 2026-10-21 — v8.6.6 custom-body surface reconstruction; animation
%           tick floor 0.5 um and a short waist roll mark
%  Updated: 2026-10-30 — v8.7.1 flow-field figure: blank interior, cell-sized
%           axes, smooth outline, vertical colour bar beside panel d
%  Updated: 2026-10-29 — v8.7.0 saved flow-field history and post-processing;
%           one mobility assembly per flow frame
%  Log banner version: 8.7.1 (2026-10-30) - emitted via project_version(); keep README / CHANGELOG in sync
%
%  Path policy: defensive addpath calls here allow
%  run_simulation_HPC_Linux(cfg) without main_coupled_simulation.m.
%  Standalone drivers in batch_sweep/ and post_processing/ resolve and add
%  the project path independently.

proj_root = fileparts(mfilename('fullpath'));
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));

%% ---- Input parsing ----
if ischar(cfg_or_path) || isstring(cfg_or_path)
    cfg = load_sim_config(cfg_or_path);
else
    cfg = cfg_or_path;
end

results    = [];
params     = [];
output_dir = '';   %#ok<NASGU> documented default if aborted before allocation

%% ---- Unpack cfg to local variables ----
case_name              = cfg.case_name;
geometry_mode          = get_cfg_field(cfg, 'geometry_mode', 'parametric');

% Silent mode: when true, logf writes to the log file only (no console echo).
% Used by batch parallel runs to prevent N_workers stdout interleaving.
% Not a physical parameter — omit from JSON; inject at runtime only.
silent                 = logical(get_cfg_field(cfg, 'silent', false));

% Optional live-progress heartbeat (HPC_Linux addition, 2026-07-03): when a
% batch driver (e.g. batch_repeat_brownian_HPC_Linux) needs real-time
% sub-case progress under cfg.silent = true -- where console echo is
% otherwise fully suppressed and the only visibility is "this repeat
% finished" -- it may pass a parallel.pool.DataQueue handle here plus a
% tag identifying this run (e.g. the repeat index k).  When present, the
% ODE segment loop below periodically sends a throttled progress message
% through the queue.  Not a physical parameter — omit from JSON; inject at
% runtime only.  Absent/empty by default: fully backward compatible with
% standalone/interactive use and with the original run_simulation.m.
progress_queue         = get_cfg_field(cfg, 'progress_queue', []);
progress_tag           = get_cfg_field(cfg, 'progress_tag', []);

% Layered silence (v7.9.4): in addition to the dual_log console-echo gate above,
% suppress (a) direct fprintf calls in this function, (b) odetpbar/render_progress_bar
% progress bars, and (c) warnings emitted by src/ sub-functions (~24 distinct IDs
% across 16 files).  Enumerating every warning ID is fragile, so we toggle 'all'
% and restore via onCleanup so the caller's warning settings are always preserved
% (normal return, error, or Ctrl-C).  Errors are NOT suppressed.
if silent
    warn_state_orig = warning('off', 'all');
    cleanup_warn = onCleanup(@() warning(warn_state_orig));
end

% Body geometry (parametric mode)
if strcmp(geometry_mode, 'parametric')
    minor_axis             = cfg.minor_axis;
    major_axis             = cfg.major_axis;
    N_body                 = cfg.N_body;
else
    minor_axis             = NaN;
    major_axis             = NaN;
    N_body                 = NaN;
end
alpha_0                = cfg.alpha_0;
Pili_Matrix            = cfg.Pili_Matrix;
seg_per_um             = cfg.seg_per_um;
E                      = cfg.E;
r_pili                 = cfg.r_pili;
hydro_radius_factor    = get_cfg_field(cfg, 'hydro_radius_factor', sqrt(2));
sigma_p                = cfg.sigma_p;
mu                     = cfg.mu;
rho_bac                = cfg.rho_bac;
rho_fluid              = cfg.rho_fluid;
T_obs                  = cfg.T_obs;
AbsTol                 = cfg.AbsTol;
RelTol                 = cfg.RelTol;
delta                  = cfg.delta;
tlim                   = cfg.tlim;
dt                     = cfg.dt;
enable_brown           = cfg.enable_brown;
Temper                 = cfg.Temper;
dt_B                   = cfg.dt_B;
preview_initial        = cfg.preview_initial;
generate_animation     = cfg.generate_animation;
animation_fps          = cfg.animation_fps;
animation_skip         = cfg.animation_skip;
auto_confirm           = cfg.auto_confirm;
flowfield_grid_spacing = cfg.flowfield_grid_spacing;
flowfield_grid_range   = cfg.flowfield_grid_range;
flowfield_quiver_skip  = cfg.flowfield_quiver_skip;
flowfield_chunk_size   = cfg.flowfield_chunk_size;

if isfield(cfg, 'compute_flowfield')
    compute_flowfield = logical(cfg.compute_flowfield);
else
    compute_flowfield = true;   % backward-compatible default
end

% New feature fields (defaults for backward compatibility)
if isfield(cfg, 'pili_ext') && ~isempty(cfg.pili_ext)
    intrinsic_curvature_cfg = {cfg.pili_ext.intrinsic_curvature};
    active_moment_cfg       = {cfg.pili_ext.active_moment};
    if isfield(cfg.pili_ext, 'root_coupling')
        root_coupling_cfg = {cfg.pili_ext.root_coupling};
    else
        root_coupling_cfg = {};
    end
    if isfield(cfg.pili_ext, 'motor_torque')
        motor_torque_cfg = {cfg.pili_ext.motor_torque};
    else
        motor_torque_cfg = {};
    end
else
    intrinsic_curvature_cfg = {};
    active_moment_cfg       = {};
    root_coupling_cfg       = {};
    motor_torque_cfg        = {};
end
if isfield(cfg, 'background_flow')
    background_flow_cfg = cfg.background_flow;
else
    background_flow_cfg = struct('type', 'off');
end

if isfield(cfg, 'eps_body') && ~isempty(cfg.eps_body)
    eps_body_cfg = cfg.eps_body;
else
    eps_body_cfg = [];
end

% Custom geometry config
if strcmp(geometry_mode, 'custom') && isfield(cfg, 'custom_geometry')
    custom_geometry_cfg = cfg.custom_geometry;
else
    custom_geometry_cfg = struct();
end

% Initial orientation (C2): custom-mode-only general quaternion/axis/Euler override
if isfield(cfg, 'initial_orientation') && ~isempty(cfg.initial_orientation)
    initial_orientation_cfg = cfg.initial_orientation;
else
    initial_orientation_cfg = [];
end

% ---- Brownian RNG: isolated stream, independent of geometry/cache RNG ----
brown_rng_seed_used = NaN;
brown_rng_stream    = [];
if enable_brown
    if isfield(cfg, 'brown_rng_seed') && isscalar(cfg.brown_rng_seed) ...
            && isnumeric(cfg.brown_rng_seed) && isfinite(cfg.brown_rng_seed) ...
            && cfg.brown_rng_seed > 0
        brown_rng_seed_used = round(double(cfg.brown_rng_seed));
    else
        brown_rng_seed_used = brown_rng_seed_from_clock();
    end
    brown_rng_stream = RandStream('mt19937ar', 'Seed', brown_rng_seed_used);
    cfg.brown_rng_seed = brown_rng_seed_used;
end

if strcmp(geometry_mode, 'parametric') && size(Pili_Matrix, 2) > 0
    pili_lengths_um = Pili_Matrix(2, :);
    N_segments      = round(pili_lengths_um * seg_per_um);
    N_segments      = max(N_segments, 10);
else
    N_segments      = zeros(1, 0);
end
NT = round(T_obs / dt) + 1;

tic_sim_start = tic;

%% =====================================================================
%  OUTPUT DIRECTORY AND LOGGING SETUP
%  =====================================================================

run_timestamp = char(datetime('now', 'Format', 'uuuuMMdd''-''HHmmss'));
if isempty(strtrim(case_name))
    run_id     = ['run-' run_timestamp];
    output_dir = fullfile(pwd, 'output', run_id);
else
    run_id     = strtrim(case_name);
    output_dir = fullfile(pwd, 'output', run_id);
end

if exist(output_dir, 'dir')
    existing_files = dir(fullfile(output_dir, '*.*'));
    existing_files = existing_files(~[existing_files.isdir]);
    data_exts = {'.mat', '.csv', '.png', '.pdf', '.html', '.mp4', '.avi', '.log', '.json', '.txt'};
    has_data  = any(cellfun(@(ext) any(endsWith({existing_files.name}, ext)), data_exts));

    if has_data
        if auto_confirm
            output_dir = [output_dir '-' run_timestamp];
            run_id     = [run_id '-' run_timestamp];
            mkdir(output_dir);
            if ~silent
                fprintf('[OUTPUT] Auto-renamed to avoid overwrite: %s\n', output_dir);
            end
        else
            overwrite_choice = questdlg( ...
                sprintf(['Output directory already contains simulation data:\n' ...
                         '  %s\n\n' ...
                         'How would you like to proceed?'], output_dir), ...
                'Output Directory Exists', ...
                'Overwrite', 'Rename (append timestamp)', 'Abort', ...
                'Rename (append timestamp)');
            switch overwrite_choice
                case 'Overwrite'
                    for ef = 1:numel(existing_files)
                        delete(fullfile(output_dir, existing_files(ef).name));
                    end
                    if ~silent
                        fprintf('[OUTPUT] Cleared %d existing file(s) and proceeding.\n', ...
                            numel(existing_files));
                    end
                case 'Rename (append timestamp)'
                    output_dir = [output_dir '-' run_timestamp];
                    run_id     = [run_id '-' run_timestamp];
                    mkdir(output_dir);
                    if ~silent
                        fprintf('[OUTPUT] Renamed to: %s\n', output_dir);
                    end
                otherwise
                    if ~silent
                        fprintf('Simulation aborted: output directory conflict.\n');
                    end
                    return;
            end
        end
    end
else
    mkdir(output_dir);
end

% --- Initialise dual-channel logger (console + file) ---
log_file = fullfile(output_dir, 'simulation.log');
log_fid  = fopen(log_file, 'w', 'n', 'UTF-8');
if log_fid == -1
    error('Cannot create log file: %s', log_file);
end

% Guarantee the log file handle is always released -- even on error or
% Ctrl-C interruption -- regardless of where execution stops below.
% onCleanup fires on every exit path (normal return, error, Ctrl-C),
% unlike try/catch which cannot intercept Ctrl-C.  safe_fclose() is a
% no-op if the fid has already been closed by the normal-path fclose()
% near the end of this function, so registering it here is harmless.
cleanup_log_fid = onCleanup(@() safe_fclose(log_fid));

% Bugfix: the trailing flag forwarded to dual_log must be the
% "echo to console" flag, i.e. ~silent (echo unless silent), not silent
% itself.  The previous `..., silent)` passed the flag un-negated, which
% meant normal (non-silent, silent=false) runs echoed NOTHING through
% logf to the console -- only simulation.log received the ~100 [INIT]/
% [PARAMS]/[ODE]/... messages -- while cfg.silent=true (batch/parallel)
% runs echoed everything, the exact opposite of the documented contract.
logf = @(varargin) dual_log(log_fid, varargin{:}, ~silent);  % echo unless silent

% Log-file-only variant (no console echo) for high-frequency diagnostics
logf_file = @(varargin) log_file_only(log_fid, varargin{:});

% Wrap the remainder of the run in try/catch so any runtime error
% (numerical divergence, out-of-memory, invalid geometry, etc.) is
% recorded in simulation.log with a full stack trace before propagating
% to the caller, instead of silently truncating the log file.  Ctrl-C is
% not caught here (MATLAB does not route it through try/catch) but is
% still handled by cleanup_log_fid above.  Early `return`s further below
% (output-directory conflict, user-aborted preview) are unaffected: a
% `return` simply exits this try block normally without touching catch.
try

logf('\n');
logf('================================================================\n');
logf('  Bacteria Elastohydrodynamic Simulation\n');
logf(['  Version ', project_version(), '\n']);
logf('================================================================\n');
logf('  Run ID    : %s\n', run_id);
logf('  Timestamp : %s\n', run_timestamp);
logf('  Output dir: %s\n', output_dir);
logf('  Platform  : MATLAB %s on %s\n', version, computer);
logf('================================================================\n\n');

% --- Save configuration snapshot for reproducibility ---
config_json_path = fullfile(output_dir, 'config.json');
save_sim_config(cfg, config_json_path);
logf('[CONFIG] Configuration saved: %s\n', config_json_path);

if enable_brown
    seed_txt = fullfile(output_dir, 'brownian_rng_seed.txt');
    fid_seed = fopen(seed_txt, 'w', 'n', 'UTF-8');
    if fid_seed ~= -1
        fprintf(fid_seed, '%d\n', brown_rng_seed_used);
        fprintf(fid_seed, ...
            'Seed for the isolated mt19937ar Brownian stream. Re-run with the same seed via cfg.brown_rng_seed or JSON brownian.rng_seed.\n');
        fclose(fid_seed);
        logf('[CONFIG] Brownian RNG seed file: %s\n', seed_txt);
    end
end

% --- Log all parameters ---
if strcmp(geometry_mode, 'parametric')
    logf('[PARAMS] Body: a = %.2f um, b = %.2f um, N_body = %d, alpha_0 = %.4f rad (%.1f deg)\n', ...
        major_axis*1e6, minor_axis*1e6, N_body, alpha_0, rad2deg(alpha_0));
else
    logf('[PARAMS] Body: N_body from file, alpha_0 = %.4f rad (%.1f deg)\n', ...
        alpha_0, rad2deg(alpha_0));
    logf('[PARAMS] Custom body: %s\n', custom_geometry_cfg.body_points_file);
end
logf('[PARAMS] Geometry mode: %s\n', geometry_mode);
logf('[PARAMS] Pili: M = %d, seg_per_um = %d, method = RSS+Rotlet\n', ...
    size(Pili_Matrix, 2), seg_per_um);
for kk = 1:size(Pili_Matrix, 2)
    logf('[PARAMS]   Pili %d: type=%d, L=%.1f um, site=%d, N_seg=%d\n', ...
        kk, Pili_Matrix(1,kk), Pili_Matrix(2,kk), ...
        Pili_Matrix(8,kk), N_segments(kk));
    % v5.0+: per-pilus intrinsic curvature & active moment summary
    if kk <= numel(intrinsic_curvature_cfg) && ~isempty(intrinsic_curvature_cfg{kk})
        ic_kk = intrinsic_curvature_cfg{kk};
        if isstruct(ic_kk) && isfield(ic_kk, 'type')
            logf('[PARAMS]     intrinsic_curvature: type = %s\n', char(ic_kk.type));
        end
    end
    if kk <= numel(active_moment_cfg) && ~isempty(active_moment_cfg{kk})
        am_kk = active_moment_cfg{kk};
        if isstruct(am_kk) && isfield(am_kk, 'type')
            am_type = char(am_kk.type);
            switch lower(am_type)
                case 'traveling_wave'
                    amp_str  = fmt_field_value(am_kk, 'amplitude_Nm');
                    freq_str = fmt_field_value(am_kk, 'omega_Hz');
                    wvl_str  = fmt_field_value(am_kk, 'wavenumber');
                    logf(['[PARAMS]     active_moment: traveling_wave | ', ...
                          'amp = %s Nm, omega = %s Hz, wavenumber = %s\n'], ...
                          amp_str, freq_str, wvl_str);
                case {'off', 'custom'}
                    logf('[PARAMS]     active_moment: %s\n', am_type);
                otherwise
                    logf('[PARAMS]     active_moment: %s\n', am_type);
            end
        end
    end
end
% v5.0+: global background flow summary
bg_type = 'off';
if isstruct(background_flow_cfg) && isfield(background_flow_cfg, 'type')
    bg_type = char(background_flow_cfg.type);
end
switch lower(bg_type)
    case 'uniform'
        umax_str = fmt_field_value(background_flow_cfg, 'U_max');
        dir_vec  = get_field_default(background_flow_cfg, 'direction', [1, 0, 0]);
        logf('[PARAMS] Background flow: uniform | U_max = %s m/s | dir = [%.3f, %.3f, %.3f]\n', ...
            umax_str, dir_vec(1), dir_vec(2), dir_vec(3));
    case 'shear_z'
        umax_str = fmt_field_value(background_flow_cfg, 'U_max');
        dir_vec  = get_field_default(background_flow_cfg, 'direction',  [1, 0, 0]);
        ax_vec   = get_field_default(background_flow_cfg, 'shear_axis', [0, 0, 1]);
        h_str    = fmt_field_value(background_flow_cfg, 'shear_height_m');
        logf(['[PARAMS] Background flow: shear_z | U_max = %s m/s | ', ...
              'dir = [%.3f, %.3f, %.3f] | axis = [%.3f, %.3f, %.3f] | h = %s m\n'], ...
              umax_str, dir_vec(1), dir_vec(2), dir_vec(3), ...
              ax_vec(1), ax_vec(2), ax_vec(3), h_str);
    otherwise
        logf('[PARAMS] Background flow: off\n');
end
logf('[PARAMS] Material: E = %.2e Pa, r_pili = %.1f nm, sigma = %.1f\n', ...
    E, r_pili*1e9, sigma_p);
logf('[PARAMS] Pilus hydrodynamic radius: eps_pili = %.4f * r_pili = %.2f nm\n', ...
    hydro_radius_factor, hydro_radius_factor*r_pili*1e9);
logf('[PARAMS] Fluid: mu = %.2e Pa*s, rho_bac = %.1f kg/m^3, rho_fluid = %.3f kg/m^3\n', ...
    mu, rho_bac, rho_fluid);
logf('[PARAMS] Simulation: T_obs = %.1f s, dt = %.2e s, NT = %d, AbsTol = %.0e, RelTol = %.0e\n', ...
    T_obs, dt, NT, AbsTol, RelTol);
logf('[PARAMS] Numerics: delta = pi/%.0f rad, timeout = %d s\n', pi/delta, tlim);
logf('[PARAMS] Preview: %s\n', mat2str(preview_initial));
logf('[PARAMS] Animation: fps = %d, skip = %d\n', animation_fps, animation_skip);
logf('[PARAMS] Flowfield: grid = %.2f um, range = %.1f um, quiver_skip = %d\n', ...
    flowfield_grid_spacing, flowfield_grid_range, flowfield_quiver_skip);
logf('[PARAMS] Brownian: enable = %s, Temper = %.2f K, dt_B = %.4f s, scheme = constant thermal force per interval (implicit Ito drift)\n', ...
    mat2str(enable_brown), Temper, dt_B);
if enable_brown
    logf('[PARAMS] Brownian RNG: seed = %d (isolated mt19937ar), saved for reproducibility\n', ...
        brown_rng_seed_used);
end
logf('[PARAMS] Auto-confirm: %s\n', mat2str(auto_confirm));
if isempty(strtrim(case_name))
    logf('[PARAMS] Case name: %s (auto-generated from timestamp)\n', run_id);
else
    logf('[PARAMS] Case name: %s\n', strtrim(case_name));
end

%% =====================================================================
%  PARAMETER VALIDATION
%  =====================================================================

logf('\n[VALIDATE] Running parameter sanity checks...\n');
% Pass the full cfg struct directly to the validator.  All field accesses
% inside validate_simulation_params are guarded by isfield(), so optional
% or mode-specific fields (custom_geometry, initial_orientation, delta,
% alpha_0, etc.) are handled gracefully without a duplicated field list
% here.  Adding a new parameter therefore only requires updating the
% validator, not this call site.
validate_simulation_params(cfg, logf);
logf('[VALIDATE] All parameter checks passed.\n');

%% =====================================================================
%  INITIALISATION
%  =====================================================================

logf('\n[INIT] Initialising simulation parameters...\n');
t_init_start = tic;

params = init_coupled_params( ...
    'geometry_mode',    geometry_mode, ...
    'custom_geometry',  custom_geometry_cfg, ...
    'N_body',       N_body, ...
    'minor_axis',   minor_axis, ...
    'major_axis',   major_axis, ...
    'alpha_0',      alpha_0, ...
    'initial_orientation', initial_orientation_cfg, ...
    'mu',           mu, ...
    'E',            E, ...
    'r_pili',       r_pili, ...
    'hydro_radius_factor', hydro_radius_factor, ...
    'sigma',        sigma_p, ...
    'rho_bac',      rho_bac, ...
    'rho_fluid',    rho_fluid, ...
    'N_segments',   N_segments, ...
    'Pili_Matrix',  Pili_Matrix, ...
    'seg_per_um',   seg_per_um, ...
    'delta',        delta, ...
    'Temper',       Temper, ...
    'enable_brown', enable_brown, ...
    'intrinsic_curvature', intrinsic_curvature_cfg, ...
    'active_moment',       active_moment_cfg, ...
    'root_coupling',       root_coupling_cfg, ...
    'motor_torque',        motor_torque_cfg, ...
    'background_flow',     background_flow_cfg, ...
    'eps_body',            eps_body_cfg, ...
    'logf',                logf, ...
    'logf_file',           logf_file);

if enable_brown
    params.brown_rng_seed = brown_rng_seed_used;
end

% Logger handles are now stored inside params by init_coupled_params.
% They will be removed before saving to .mat (see [SAVE] section).

logf('[INIT] Building initial state vector...\n');
Z0 = init_coupled_state(params, logf);
numerical_diagnostics_state('reset');
initial_health = check_coupled_system_health(Z0, params, 'initial state');
logf(['[DIAG] Initial rcond = %.3e | reciprocity = %.3e | ', ...
      'min dissipation eigenvalue = %.3e\n'], ...
    initial_health.rcond, initial_health.weighted_reciprocity_residual, ...
    initial_health.min_dissipation_eigenvalue);

t_init_elapsed = toc(t_init_start);

% --- Derived quantities from parameter structure ---
L_0_dim = params.scales.L_0;
F_sink_N = params.Eg * params.scales.F_0;
V_body_m3 = params.dim.volume;

% Log derived scales
logf('[INIT] Characteristic scales:\n');
logf('[INIT]   L_0 = %.4e m  (length)\n',           params.scales.L_0);
logf('[INIT]   T_0 = %.4e s  (time)\n',             params.scales.T_0);
logf('[INIT]   U_0 = %.4e m/s  (velocity)\n',       params.scales.U_0);
logf('[INIT]   F_0 = %.4e N  (force)\n',            params.scales.F_0);
logf('[INIT]   EI  = %.4e N*m^2  (bending stiffness)\n', params.scales.EI);
logf('[INIT] Dimensionless parameters:\n');
logf('[INIT]   Eg      = %.6f  (elastogravitational number)\n', params.Eg);
logf('[INIT]   epsilon = %.4e  (pilus regularisation length = %.4f * r_pili / L_0)\n', ...
    params.epsilon, params.hydro_radius_factor);
if ~isempty(eps_body_cfg)
    logf('[INIT]   eps_body (RSM override, L0): %.6f  (h_mean = %.6f)\n', ...
        params.eps_body, params.body_h_mean);
end
logf('[INIT]   alpha_0 = %.4f rad (%.1f deg) -initial body tilt\n', ...
    params.alpha_0, rad2deg(params.alpha_0));
logf('[INIT]   F_sink = %.2e N, F_0 = %.2e N\n', F_sink_N, params.scales.F_0);
logf('[INIT]   V_body = %.4e m^3\n', V_body_m3);
logf('[INIT]   q_0 = [%.4f; %.4f; %.4f; %.4f]\n', Z0(4), Z0(5), Z0(6), Z0(7));
logf('[INIT] State vector : %d DOFs\n', length(Z0));
logf('[INIT] Total points : %d (body + pili)\n', ...
    params.N_body + sum(params.N_segments + 1));
logf('[INIT] -------------------------------------------------------------\n');
logf('[INIT] Mobility Method : RSM (body) + RSS (pili) + Rotlet  [Cortez 2018]\n');
logf('[INIT] Regularisation Parameters\n');
logf('[INIT]   Body spacing : h_min = %.4f L0, h_mean = %.4f L0\n', ...
    params.body_h_min, params.body_h_mean);
logf('[INIT]   Body  : eps_body = %.4f L0 = %.1f nm (eps/h_mean = %.4f)\n', ...
    params.eps_body, params.eps_body * L_0_dim * 1e9, params.eps_body / params.body_h_mean);
if params.M_pili > 0
    for kk = 1:params.M_pili
        ds_k    = params.ell(kk) / params.N_segments(kk);
        ds_nm   = ds_k * L_0_dim * 1e9;
        eps_nm  = params.epsilon * L_0_dim * 1e9;
        logf('[INIT]   Pili %d: ell = %.4f L_0, N_seg = %d, ds = %.4f L_0 (%.2f nm)\n', ...
            kk, params.ell(kk), params.N_segments(kk), ds_k, ds_nm);
        logf('[INIT]          eps_pili = %.4f L0 (%.1f nm = %.4f r_pili), ds/eps = %.1f\n', ...
            params.epsilon, eps_nm, params.hydro_radius_factor, ds_nm / eps_nm);
    end
    dis_dim_nm = params.dim.dis * 1e9;
    logf('[INIT]   Gap : dis = 1.1*(eps_body + r_pili/L0) = %.1f nm\n', dis_dim_nm);
    logf('[INIT]   Note: eps_ij = eps_j (determined by force point j)\n');
    logf('[INIT]     Body force -> any point: eps = %.1f nm\n', ...
        params.eps_body * L_0_dim * 1e9);
    logf('[INIT]     Pili force -> any point: eps = %.1f nm (hydrodynamic radius)\n', ...
        params.epsilon * L_0_dim * 1e9);
end
logf('[INIT] -------------------------------------------------------------\n');
logf('[INIT] Background flow : %s', params.bg_flow.type);
if ~strcmp(params.bg_flow.type, 'off')
    logf('  (U_max_dimless = %.2e)', params.bg_flow.U_max_dimless);
end
logf('\n');
if params.M_pili > 0
    for kk = 1:params.M_pili
        ic_kk_type = 'frozen';
        if kk <= numel(intrinsic_curvature_cfg) && ~isempty(intrinsic_curvature_cfg{kk})
            if isstruct(intrinsic_curvature_cfg{kk}), ic_kk_type = intrinsic_curvature_cfg{kk}.type;
            elseif ischar(intrinsic_curvature_cfg{kk}) || isstring(intrinsic_curvature_cfg{kk}), ic_kk_type = char(intrinsic_curvature_cfg{kk});
            end
        end
        am_kk_type = params.am_cfg{kk}.type;
        rc_kk = params.root_coupling{kk};
        mt_kk_type = params.motor_cfg{kk}.type;
        if strcmp(rc_kk, 'clamped')
            logf('[INIT]   Pili %d: intrinsic_curvature = %s, active_moment = %s, root_coupling = clamped\n', ...
                kk, ic_kk_type, am_kk_type);
        else
            logf('[INIT]   Pili %d: intrinsic_curvature = %s, active_moment = %s, root_coupling = free_torsion, motor_torque = %s\n', ...
                kk, ic_kk_type, am_kk_type, mt_kk_type);
            if strcmp(mt_kk_type, 'constant')
                logf('[INIT]          motor: torque = %.2f pN*nm, max_omega = %.1f Hz (dimless: tau=%.4e, w_max=%.4e)\n', ...
                    params.motor_cfg{kk}.torque_pN_nm, params.motor_cfg{kk}.max_omega_Hz, ...
                    params.motor_cfg{kk}.torque_dimless, params.motor_cfg{kk}.max_omega_dimless);
            end
        end
    end
end
logf('[INIT] Initialisation completed in %.2f s\n\n', t_init_elapsed);

%% =====================================================================
%  INITIAL POSE PREVIEW (optional, interactive)
%  =====================================================================

if preview_initial
    % Guard against a headless deadlock: under cfg.silent=true (batch /
    % parallel workers), a modal questdlg cannot be answered by anyone,
    % so a worker would hang forever with no console/log indication of
    % why.  Fail fast with a clear error instead.
    if silent && ~auto_confirm
        error('run_simulation_HPC_Linux:SilentPreviewDeadlock', ...
            ['cfg.silent = true with preview_initial = true requires ', ...
             'auto_confirm = true (an interactive questdlg cannot be ', ...
             'answered on a headless/parallel worker). Set auto_confirm ', ...
             '= true or preview_initial = false.']);
    end

    logf('[PREVIEW] Generating initial pose scatter plot...\n');

    % Reconstruct initial geometry and convert to um
    geom_preview = reconstruct_geometry(Z0, params);
    L_to_um_prev = params.scales.L_0 * 1e6;

    X_body_prev  = geom_preview.X_body * L_to_um_prev;     % [N_body x 3]
    X_pili_prev  = cell(params.M_pili, 1);
    for kp = 1:params.M_pili
        X_pili_prev{kp} = geom_preview.X_pili{kp} * L_to_um_prev;
    end

    % Per-pili distinguishable colours (same palette as animation section)
    prev_pili_palette = [0.20 0.40 0.80;   % Medium blue
                         0.85 0.33 0.10;   % Burnt orange
                         0.47 0.67 0.19;   % Olive green
                         0.64 0.08 0.18;   % Crimson
                         0.49 0.18 0.56;   % Purple
                         0.93 0.69 0.13;   % Goldenrod
                         0.30 0.75 0.93;   % Sky blue
                         0.60 0.60 0.60;   % Grey
                         0.00 0.45 0.74;   % Dark blue
                         0.85 0.47 0.66];  % Rose pink
    prev_pili_colors = prev_pili_palette( ...
        mod((1:params.M_pili)-1, size(prev_pili_palette,1)) + 1, :);

    body_preview_color = [0.20 0.65 0.33];   % Forest green

    % Compute per-axis bounding box
    all_pts = X_body_prev;
    for kp = 1:params.M_pili
        all_pts = [all_pts; X_pili_prev{kp}]; %#ok<AGROW>
    end
    pv_min = min(all_pts);
    pv_max = max(all_pts);
    pv_span = pv_max - pv_min;
    pv_margin = max(0.10 * pv_span, 0.15);   % [1x3] per-axis margin

    % Create figure with 3 views.  Visibility follows `silent` (belt-and-
    % braces alongside the deadlock guard above) so a batch/parallel
    % worker never pops a real window even if auto_confirm was set.
    fig_preview_vis = 'on';
    if silent
        fig_preview_vis = 'off';
    end
    fig_preview = figure('Name', 'Initial Pose Preview', ...
                         'NumberTitle', 'off', ...
                         'Position', [80, 80, 1500, 500], ...
                         'Color', 'w', 'Visible', fig_preview_vis);
    set(fig_preview, 'DefaultAxesFontName',  'Arial', ...
                     'DefaultAxesFontSize',  11, ...
                     'DefaultAxesLineWidth', 0.8);

    view_configs = {45, 30, '3D View'; ...
                     0,  0, 'XZ View (side)'; ...
                    90,  0, 'YZ View (front)'};

    for vi = 1:3
        ax_pv = subplot(1, 3, vi);
        hold(ax_pv, 'on');

        % Body scatter
        scatter3(ax_pv, X_body_prev(:,1), X_body_prev(:,2), X_body_prev(:,3), ...
                 6, body_preview_color, 'filled', 'MarkerFaceAlpha', 0.5);

        % Pili scatter (per-pili colour)
        legend_entries = {'Body'};
        for kp = 1:params.M_pili
            Xp = X_pili_prev{kp};
            scatter3(ax_pv, Xp(:,1), Xp(:,2), Xp(:,3), ...
                     12, prev_pili_colors(kp,:), 'filled');
            legend_entries{end+1} = sprintf('Pili %d', kp); %#ok<AGROW>
        end

        % Formatting
        pv_lim_x = [pv_min(1)-pv_margin(1), pv_max(1)+pv_margin(1)];
        pv_lim_y = [pv_min(2)-pv_margin(2), pv_max(2)+pv_margin(2)];
        pv_lim_z = [pv_min(3)-pv_margin(3), pv_max(3)+pv_margin(3)];

        % Enforce minimum axis range
        if strcmp(geometry_mode, 'parametric')
            pv_min_range = 4 * 2 * (major_axis * 1e6);
        else
            pv_min_range = 4 * max(max(X_body_prev) - min(X_body_prev));
        end
        if diff(pv_lim_x) < pv_min_range
            cx = mean(pv_lim_x); pv_lim_x = cx + [-1, 1] * pv_min_range / 2;
        end
        if diff(pv_lim_y) < pv_min_range
            cy = mean(pv_lim_y); pv_lim_y = cy + [-1, 1] * pv_min_range / 2;
        end
        if diff(pv_lim_z) < pv_min_range
            cz = mean(pv_lim_z); pv_lim_z = cz + [-1, 1] * pv_min_range / 2;
        end

        xlim(ax_pv, pv_lim_x);
        ylim(ax_pv, pv_lim_y);
        zlim(ax_pv, pv_lim_z);
        set(ax_pv, 'DataAspectRatio', [1 1 1], 'TickDir', 'out');
        grid(ax_pv, 'on');
        box(ax_pv,  'on');
        xlabel(ax_pv, 'X (μm)', 'FontName', 'Arial', 'FontSize', 12);
        ylabel(ax_pv, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 12);
        zlabel(ax_pv, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 12);
        view(ax_pv, view_configs{vi, 1}, view_configs{vi, 2});
        title(ax_pv, view_configs{vi, 3}, ...
              'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');

        if vi == 1
            lg_pv = legend(ax_pv, legend_entries, ...
                           'Location', 'northeast', 'FontSize', 9, 'Box', 'off');
            set(lg_pv, 'FontName', 'Arial');
        end

        hold(ax_pv, 'off');
    end

    sgtitle(sprintf('Initial Pose Preview  [\\alpha_0 = %.1f°,  %s]', ...
        rad2deg(alpha_0), run_id), ...
        'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');

    drawnow;

    % Add gravity indicators
    for vi = 1:3
        ax_pv_i = subplot(1, 3, vi);
        add_gravity_indicator_outside(fig_preview, ax_pv_i, 11, 1.8);
    end

    logf('[PREVIEW] Initial pose figure displayed.\n');
    logf('[PREVIEW] Body: %d points (green)  |  Pili: %d root(s)\n', ...
        params.N_body, params.M_pili);
    % Interactive confirmation dialog (skipped in headless mode)
    if auto_confirm
        user_choice = 'Proceed';
        logf('[PREVIEW] Auto-confirm enabled -skipping interactive dialog.\n');
    else
        logf('[PREVIEW] Waiting for user confirmation...\n');
        user_choice = questdlg( ...
            sprintf(['Initial pose preview is shown.\n\n' ...
                     'alpha_0 = %.4f rad (%.1f deg)\n' ...
                     'Body points: %d\n' ...
                     'Pili: %d root(s)\n\n' ...
                     'Proceed with simulation?'], ...
                     alpha_0, rad2deg(alpha_0), params.N_body, params.M_pili), ...
            'Initial Pose Confirmation', ...
            'Proceed', 'Abort', 'Proceed');
    end

    if strcmp(user_choice, 'Abort') || isempty(user_choice)
        logf('[PREVIEW] User aborted simulation.\n');
        close(fig_preview);
        fclose(log_fid);
        if ~silent
            fprintf('Simulation aborted by user at initial pose preview.\n');
        end
        return;
    end

    % Save preview figure to disk before closing (raster + vector + HTML on R2026a+)
    preview_png  = fullfile(output_dir, 'initial_pose.png');
    preview_pdf  = fullfile(output_dir, 'initial_pose.pdf');
    preview_html = fullfile(output_dir, 'initial_pose.html');
    export_figure_png_pdf_html(fig_preview, preview_png, preview_pdf, preview_html);
    logf('[PREVIEW] Initial pose saved : %s  |  %s\n', preview_png, preview_pdf);
    if isfile(preview_html)
        logf('[PREVIEW] Interactive HTML : %s\n', preview_html);
    end
    logf('[PREVIEW] User confirmed. Closing preview and continuing...\n\n');
    close(fig_preview);
end

%% =====================================================================
%  MAIN SIMULATION LOOP (with adaptive basis rotation)
%  =====================================================================

logf('================================================================\n');
logf('  STARTING ODE INTEGRATION\n');
logf('================================================================\n\n');

T_0            = params.scales.T_0;
T_obs_dimless  = T_obs / T_0;
ts_dimless     = linspace(0, T_obs_dimless, NT);

logf('[ODE] Dimensionless time span : [0, %.4e]\n', T_obs_dimless);
logf('[ODE] Physical time span      : [0, %.2f] s\n', T_obs);
logf('[ODE] Output time points      : %d (dt = %.4e s)\n\n', NT, T_obs/(NT-1));

% --- Pre-allocate result storage ---
M_pili = params.M_pili;
N_segs = params.N_segments;

results.t_dimless = ts_dimless;
results.X_cm      = zeros(3, NT);
results.q         = zeros(4, NT);

[X_cm0, q0, ~]       = unpack_state(Z0, N_segs);
results.X_cm(:, 1)   = X_cm0;
results.q(:, 1)      = q0;

results.X_pili = cell(M_pili, 1);
geom0 = reconstruct_geometry(Z0, params);
for k = 1:M_pili
    results.X_pili{k} = zeros(N_segs(k)+1, 3, NT);
    results.X_pili{k}(:, :, 1) = geom0.X_pili{k};
end

if enable_brown
    results.brown_rng_seed = brown_rng_seed_used;
end

% Flow-field reconstruction history. Saved on every run, whether or not a
% flow field is drawn. brown_force(:,i) is the one thermal force that
% produced output i, including the last output. It is zero at the initial
% sample. The saved parameter value is cleared after integration; the
% final flow field uses brown_force(:,end), not that cleared value.
nZ_flow = numel(Z0);
results.flow_state.version = 1;
results.flow_state.Z = zeros(nZ_flow, NT);
results.flow_state.Z(:, 1) = Z0(:);
results.flow_state.brown_force = zeros(6, NT);
results.flow_state.rot = cell(M_pili, 1);
for k = 1:M_pili
    results.flow_state.rot{k} = zeros(3, 3, NT);
    results.flow_state.rot{k}(:, :, 1) = params.rot{k};
end
results.flow_state.units.Z = 'dimensionless ODE state (lengths in L_0, angles in rad)';
results.flow_state.units.rot = 'pilus local basis; rows map lab-frame vectors into the basis';
results.flow_state.units.brown_force = ['dimensionless generalised thermal force, ', ...
    'same units as params.brown_force_dimless'];
results.flow_state.brownian_endpoint_convention = 'closing_interval';
results.flow_state.state_dimension = nZ_flow;
results.flow_state.n_times = NT;

% --- Integration loop ---
clear odetpbar;                    % Reset persistent vars from prior runs
clear ode_coupled_rhs;             % Reset DPSI1_PREV motor-speed cache
tic_total        = tic;
t_last_heartbeat = tic;   % throttle timer for the optional progress_queue heartbeat
T_achieved       = 0;
last_Z           = Z0;
current_time_ind = 1;
basis_changes    = 0;
segment_count    = 0;
seg_walltimes    = [];          % Wall time per segment [s]
seg_dtstar       = [];          % Dimensionless time advanced per segment
termination_reason = 'completed';
timed_out          = false;

% ---- Lagged Jacobian for ode15s (v8.5.1; Brownian runs only) ----
%  Brownian runs restart ode15s at every dt_B node and each restart pays a
%  dense finite-difference Jacobian (n_ode extra RHS calls).  The lagged
%  Jacobian re-uses the last exact FD Jacobian (see LaggedJacobian.m); it
%  affects only the Newton iteration matrix, not the discrete equations or
%  the error control.  Runtime option cfg.jacobian_reuse (default true);
%  set false to reproduce the previous behaviour bit-for-bit.
%  cfg.jacobian_max_age (default 50) caps the cache age in Brownian nodes.
jac_reuse_enabled = params.brown_enabled && ...
    logical(get_cfg_field(cfg, 'jacobian_reuse', true));
jac_cache = [];
if jac_reuse_enabled
    jac_cache = LaggedJacobian(get_cfg_field(cfg, 'jacobian_max_age', 50));
    logf('[JAC] Lagged Jacobian enabled (max age = %d Brownian nodes)\n', ...
        jac_cache.max_age);
end

% Segment-level [ODE] announcements are throttled to avoid flooding
% simulation.log / console when Brownian dynamics splits the run into
% many short segments (one per dt_B node -- e.g. T_obs/dt_B ~= 4000 for
% a typical long Brownian run).  Per-segment timing/ETA statistics
% (seg_walltimes, seg_dtstar) are still recorded on EVERY segment
% regardless of this throttle; only the "Segment N: ..." announcement
% text is skipped.  Non-Brownian runs have few segments (basis-rotation
% restarts only), so every segment is logged (interval = 1).
if enable_brown
    seg_log_interval = 50;
else
    seg_log_interval = 1;
end

% ---- Brownian dynamics node tracking ----
%  v8.5.0 scheme: at the start of every Brownian interval [T_brown_last,
%  T_brown_next] a constant generalised thermal force f_th is drawn from the
%  6x6 effective body mobility (sample_brownian_body_force) and stored in
%  params.brown_force_dimless; ode_coupled_rhs applies it to the body rows
%  for the whole interval.  The pili therefore follow the body's thermal
%  motion through the coupled hydrodynamic solve (no root-projection kink),
%  and the Ito drift kB*T*div(M) is generated implicitly by the mobility
%  gradient along the trajectory (no finite-difference correction).
if params.brown_enabled
    dt_B_dimless   = dt_B / params.scales.T_0;
    T_brown_last   = 0;
    T_brown_next   = min(dt_B_dimless, T_obs_dimless);
    logf(['[BROWN] Body-only Brownian dynamics enabled ', ...
          '(constant 6-DOF thermal generalised force per interval, ', ...
          'coupled effective mobility probe at each interval start)\n']);
    logf('[BROWN] dt_B = %.4f s = %.4e T*\n', dt_B, dt_B_dimless);
    logf(['[BROWN] Ito drift: IMPLICIT (constant-force realisation along the ', ...
          'coupled trajectory reproduces kB*T*div(M) to O(dt_B))\n']);
    brown_kick_count    = 0;      % completed Brownian intervals
    brown_log_interval  = 100;
    brown_step_sizes_s  = [];
    brown_force_active  = false;  % true while an interval's force is held
    brown_tol           = dt_B_dimless * 1e-8;
else
    dt_B_dimless = Inf;
    T_brown_last = 0;
    T_brown_next = Inf;
    brown_step_sizes_s = [];
    brown_force_active = false;
    brown_tol          = 0;
end
params.brown_force_dimless = zeros(6, 1);

while T_achieved < T_obs_dimless
    Z = last_Z;
    segment_count = segment_count + 1;
    log_this_segment = (segment_count <= 3) || (mod(segment_count, seg_log_interval) == 0);

    % Adaptive basis check for all pili
    [Z, params, rotated_flags] = adaptive_basis_check(Z, params);
    if any(rotated_flags)
        if jac_reuse_enabled
            jac_cache.invalidate();   % Euler angles re-parameterised
        end
        basis_changes = basis_changes + sum(rotated_flags);
        for kk = 1:length(rotated_flags)
            if rotated_flags(kk)
                logf('[ODE] Basis rotation triggered for pili %d at t* = %.4e\n', ...
                    kk, T_achieved);
            end
        end
        basis_health = check_coupled_system_health( ...
            Z, params, sprintf('post-basis-change segment %d', segment_count));
        logf('[DIAG] Post-basis rcond = %.3e | reciprocity = %.3e\n', ...
            basis_health.rcond, basis_health.weighted_reciprocity_residual);
    end

    % Determine integration target for this segment
    T_seg_target = min(T_brown_next, T_obs_dimless);

    % ---- Brownian: realise the thermal force for a NEW interval ----
    %  Executed once per Brownian interval (a basis-rotation or timeout
    %  restart inside an interval keeps the already-sampled force so that
    %  the realised increment stays consistent with the drawn deviates).
    if params.brown_enabled && ~brown_force_active && ...
            T_achieved < T_obs_dimless - brown_tol
        dt_B_step = (T_seg_target - T_achieved) * params.scales.T_0;   % [s]
        [f_th_dimless, brown_info] = sample_brownian_body_force( ...
            Z, params, dt_B_step, brown_rng_stream);
        params.brown_force_dimless = f_th_dimless;
        brown_force_active = true;
        if jac_reuse_enabled
            jac_cache.tick();         % one more Brownian node since last FD J
        end

        if brown_kick_count < 5 || mod(brown_kick_count + 1, brown_log_interval) == 0
            logf(['[BROWN] Interval #%d start at t* = %.4e | dt = %.4e s | ', ...
                  '|F_th| = %.2e N | |tau_th| = %.2e N*m | ', ...
                  'target |dX| = %.2e m | target |dTheta| = %.2e rad\n'], ...
                brown_kick_count + 1, T_achieved, dt_B_step, ...
                norm(brown_info.f_th_dim(1:3)), norm(brown_info.f_th_dim(4:6)), ...
                norm(brown_info.dU_target_dim(1:3)), norm(brown_info.dU_target_dim(4:6)));
        end
    end

    % Preserve the force acting in this segment. Interval-end bookkeeping
    % clears params.brown_force_dimless before output samples are written;
    % those samples keep the force that produced them.
    flow_interval_force = params.brown_force_dimless(:);

    % Setup ODE solver
    tstart_ode   = tic;
    eventFunc    = @(t_ev, Z_ev) singularity_event(t_ev, Z_ev, params, tstart_ode, tlim);
    % Always install the timeout OutputFcn.  In silent mode it omits only
    % the progress-bar rendering, not the wall-clock interrupt.
    progressFunc = @(t_p, y_p, fl_p, varargin) ...
        timeout_output(t_p, y_p, fl_p, tstart_ode, tlim, silent);

    ode_ops = odeset('OutputFcn',  progressFunc, ...
                     'Events',     eventFunc, ...
                     'Stats',      'off', ...
                     'AbsTol',     AbsTol, ...
                     'RelTol',     RelTol, ...
                     'Vectorized', 'off');

    dZ_func = @(t_ode, Z_ode) ode_coupled_rhs(t_ode, Z_ode, params);

    if log_this_segment
        logf('[ODE] Segment %d: integrating t* = [%.4e, %.4e]\n', ...
            segment_count, T_achieved, T_seg_target);
    end

    if jac_reuse_enabled
        % Lagged Jacobian.  Any solver error, or a premature stop that is
        % neither an event nor a timeout, re-runs the segment from the same
        % initial state with ode15s' default Jacobian (safety net).
        jac_cache.beginSegment(dZ_func);
        ode_ops_jac = odeset(ode_ops, 'Jacobian', ...
            @(t_j, Z_j) jac_cache.evaluate(t_j, Z_j));
        try
            sol = ode15s(dZ_func, [T_achieved, T_seg_target], Z, ode_ops_jac);
            jac_fallback = isempty(sol.xe) && ...
                sol.x(end) < T_seg_target - 1e-12 && toc(tstart_ode) < tlim;
        catch ME_jac
            jac_fallback = true;
            logf('[JAC] WARNING: ode15s with lagged Jacobian failed (%s); re-running the segment with the default Jacobian.\n', ...
                ME_jac.message);
        end
        if jac_fallback
            jac_cache.noteFallback();
            logf('[JAC] Segment %d re-run with default finite-difference Jacobian.\n', ...
                segment_count);
            sol = ode15s(dZ_func, [T_achieved, T_seg_target], Z, ode_ops);
        end
    else
        sol  = ode15s(dZ_func, [T_achieved, T_seg_target], Z, ode_ops);
    end
    seg_wall = toc(tstart_ode);
    if log_this_segment
        logf('[ODE] Segment %d completed (%.1f s, %d steps)\n', ...
            segment_count, seg_wall, length(sol.x));
    end

    % Check if event was triggered
    stopped_by_timeout_output = seg_wall >= tlim && ...
        sol.x(end) < T_seg_target - 1e-14;
    timeout_this_segment = stopped_by_timeout_output;
    if ~isempty(sol.xe)
        event_ids = unique(sol.ie(:)');
        timeout_this_segment = timeout_this_segment || any(event_ids == 1);
        event_labels = cell(1, numel(event_ids));
        for ei = 1:numel(event_ids)
            switch event_ids(ei)
                case 1
                    event_labels{ei} = 'timeout';
                case 2
                    event_labels{ei} = 'euler-singularity-margin';
                otherwise
                    event_labels{ei} = sprintf('unknown-event-%d', event_ids(ei));
            end
        end
        logf('[ODE]   Event triggered at t* = %.4e | reason = %s\n', ...
            sol.xe(end), strjoin(event_labels, ','));
    elseif stopped_by_timeout_output
        logf('[ODE]   OutputFcn stopped at t* = %.4e | reason = timeout\n', ...
            sol.x(end));
    end

    % Update state
    T_prev     = T_achieved;
    T_achieved = sol.x(end);
    Z_end      = sol.y(:, end);

    % Record segment timing for ETA estimation
    seg_walltimes(end+1) = seg_wall;            %#ok<AGROW>
    seg_dtstar(end+1)    = T_achieved - T_prev; %#ok<AGROW>

    % Renormalise quaternion
    q_end         = Z_end(4:7);
    q_norm_before = norm(q_end);
    Z_end(4:7)    = q_end / q_norm_before;

    % Enforce geometric clamped BC (correct constraint drift)
    Z_end = enforce_clamped_geometry(Z_end, params);

    % ---- Close the Brownian interval at a completed Brownian target ----
    %  The thermal force sampled at the interval start has now acted on the
    %  body through the coupled solve for the whole interval; the state
    %  Z_end already contains the realised (mobility-consistent) increment,
    %  so no displacement kick or root re-projection is applied here.  The
    %  force is cleared; a fresh one is drawn at the top of the next segment.
    if params.brown_enabled && ~timeout_this_segment && brown_force_active
        if T_achieved >= T_seg_target - brown_tol
            dt_B_step = (T_achieved - T_brown_last) * params.scales.T_0;
            params.brown_force_dimless = zeros(6, 1);
            brown_force_active = false;

            brown_kick_count   = brown_kick_count + 1;
            T_brown_last       = T_achieved;
            brown_step_sizes_s(end+1) = dt_B_step; %#ok<AGROW>
            if T_brown_last < T_obs_dimless - brown_tol
                T_brown_next = min(T_brown_last + dt_B_dimless, T_obs_dimless);
            else
                T_brown_next = Inf;
            end

            if brown_kick_count <= 5 || mod(brown_kick_count, brown_log_interval) == 0
                logf('[BROWN] Interval #%d closed at t* = %.4e | dt = %.4e s\n', ...
                    brown_kick_count, T_achieved, dt_B_step);
            end
        end
    end
    % ---- End Brownian interval bookkeeping ----

    last_Z = Z_end;

    % Log quaternion health
    if abs(q_norm_before - 1) > 1e-6
        logf('[ODE]   WARNING: quaternion norm drift = %.2e (renormalised)\n', ...
            abs(q_norm_before - 1));
    end

    % Evaluate at output time points
    mask     = ts_dimless <= T_achieved + 1e-14;
    valid_ts = ts_dimless(mask);
    valid_ts = valid_ts(current_time_ind+1 : end);

    if ~isempty(valid_ts)
        valid_ts = max(valid_ts, sol.x(1));
        valid_ts = min(valid_ts, sol.x(end));
        Z_at_ts = deval(sol, valid_ts);
        for i = 1:length(valid_ts)
            ti  = current_time_ind + i;
            Z_i = Z_at_ts(:, i);

            % Renormalise quaternion
            Z_i(4:7) = Z_i(4:7) / norm(Z_i(4:7));

            % Enforce geometric clamped BC
            Z_i = enforce_clamped_geometry(Z_i, params);

            [X_cm_i, q_i, ~] = unpack_state(Z_i, N_segs);
            results.X_cm(:, ti) = X_cm_i;
            results.q(:, ti)    = q_i;

            geom_i = reconstruct_geometry(Z_i, params);
            for k = 1:M_pili
                results.X_pili{k}(:, :, ti) = geom_i.X_pili{k};
                results.flow_state.rot{k}(:, :, ti) = params.rot{k};
            end
            results.flow_state.Z(:, ti) = Z_i(:);
            results.flow_state.brown_force(:, ti) = flow_interval_force;
        end
    end

    current_time_ind = find(ts_dimless <= T_achieved + 1e-14, 1, 'last');

    % A timeout can occur between requested output nodes.  Preserve the exact
    % terminal state as the final sample before truncating preallocated arrays.
    if timeout_this_segment
        if abs(results.t_dimless(current_time_ind) - T_achieved) > 1e-14
            ti = current_time_ind + 1;
            results.t_dimless(ti) = T_achieved;
            [X_cm_i, q_i, ~] = unpack_state(Z_end, N_segs);
            results.X_cm(:, ti) = X_cm_i;
            results.q(:, ti)    = q_i;
            geom_i = reconstruct_geometry(Z_end, params);
            for k = 1:M_pili
                results.X_pili{k}(:, :, ti) = geom_i.X_pili{k};
                results.flow_state.rot{k}(:, :, ti) = params.rot{k};
            end
            results.flow_state.Z(:, ti) = Z_end(:);
            results.flow_state.brown_force(:, ti) = flow_interval_force;
            current_time_ind = ti;
        end
        timed_out = true;
        termination_reason = 'timeout';
        logf('[ODE] Integration aborted after segment %d timeout at t* = %.6e.\n', ...
            segment_count, T_achieved);
    end

    % Log progress with ETA
    progress_pct = min(100, T_achieved / T_obs_dimless * 100);
    elapsed_ode  = toc(tic_total);

    if progress_pct > 0 && progress_pct < 100
        eta_ode_s = elapsed_ode * (100 - progress_pct) / progress_pct;
        logf('[ODE]   Progress: %.1f%% | Elapsed: %s | ETA: %s\n', ...
            progress_pct, format_duration(elapsed_ode), ...
            format_duration(eta_ode_s));
    else
        logf('[ODE]   Progress: %.1f%% | Elapsed: %s\n', ...
            progress_pct, format_duration(elapsed_ode));
    end

    % Optional live-progress heartbeat: reuses the progress_pct/elapsed_ode
    % just computed above, so this adds no extra computation cost.
    % Throttled to at most once per 30 s per run so a several-thousand-
    % segment Brownian integration does not flood the DataQueue with
    % thousands of messages.  No-op (near-zero overhead) when the caller
    % did not provide a progress_queue (e.g. standalone/interactive use).
    if ~isempty(progress_queue) && toc(t_last_heartbeat) >= 30
        send(progress_queue, struct('type', 'heartbeat', 'tag', progress_tag, ...
            'pct', progress_pct, 'elapsed', elapsed_ode));
        t_last_heartbeat = tic;
    end

    if timed_out
        break;
    end
end

elapsed_total = toc(tic_total);

% Clear the active parameter after the loop so a saved params struct does
% not look like an integration still in progress. The last output keeps
% its own recorded force in flow_state.brown_force(:,end).
params.brown_force_dimless = zeros(6, 1);

% Remove every uncomputed preallocated sample before diagnostics,
% post-processing, persistence, and visualisation.
results.t_dimless = results.t_dimless(1:current_time_ind);
results.X_cm      = results.X_cm(:, 1:current_time_ind);
results.q         = results.q(:, 1:current_time_ind);
results.flow_state.Z = results.flow_state.Z(:, 1:current_time_ind);
results.flow_state.brown_force = results.flow_state.brown_force(:, 1:current_time_ind);
for k = 1:M_pili
    results.X_pili{k} = results.X_pili{k}(:, :, 1:current_time_ind);
    results.flow_state.rot{k} = results.flow_state.rot{k}(:, :, 1:current_time_ind);
end
results.flow_state.n_times = current_time_ind;
NT = current_time_ind;
results.completed          = ~timed_out;
results.termination_reason = termination_reason;
results.target_time_s      = T_obs;
results.achieved_time_s    = T_achieved * params.scales.T_0;
results.segment_count      = segment_count;
results.brownian_step_sizes_s = brown_step_sizes_s;
results.diagnostics = numerical_diagnostics_state('snapshot');
if jac_reuse_enabled
    results.jacobian_reuse = jac_cache.stats();
    logf(['[JAC] Lagged Jacobian summary: fresh = %d (on-demand %d, age %d) | ', ...
          'lagged = %d | basis invalidations = %d | fallbacks = %d | ', ...
          'FD RHS calls = %d\n'], ...
        results.jacobian_reuse.n_fresh, results.jacobian_reuse.n_on_demand, ...
        results.jacobian_reuse.n_age_refresh, results.jacobian_reuse.n_lagged, ...
        results.jacobian_reuse.n_invalidated, results.jacobian_reuse.n_fallback, ...
        results.jacobian_reuse.n_fd_rhs);
else
    results.jacobian_reuse = struct('enabled', false);
end
diag = results.diagnostics;
logf(['[DIAG] Max backward residual = %.3e (warn %.1e, err %.1e) | ', ...
      'min rcond = %.3e (warn %.1e, near-sing %.1e; joint-abort with resid err) | ', ...
      'max reciprocity residual = %.3e (err %.1e) | ', ...
      'min dissipation eigenvalue = %.3e (rel err %.1e)\n'], ...
    diag.max_backward_residual, ...
    diag.backward_residual_warning_threshold, ...
    diag.backward_residual_error_threshold, ...
    diag.min_rcond, ...
    diag.rcond_warning_threshold, ...
    diag.rcond_near_singular_threshold, ...
    diag.max_weighted_reciprocity_residual, ...
    diag.reciprocity_error_threshold, ...
    diag.min_dissipation_eigenvalue, ...
    diag.dissipation_error_relative_threshold);
logf(['[DIAG] Check counts: linear_solve = %d | rcond = %d | ', ...
      'hydrodynamic = %d\n'], ...
    diag.linear_solve_count, diag.rcond_check_count, ...
    diag.hydrodynamic_check_count);

requested_generate_animation = generate_animation;
requested_compute_flowfield  = compute_flowfield;
if timed_out
    generate_animation = false;
    compute_flowfield  = false;
end

logf('\n================================================================\n');
if timed_out
    logf('  INTEGRATION ABORTED: TIMEOUT\n');
else
    logf('  INTEGRATION COMPLETED\n');
end
logf('================================================================\n');
logf('[ODE] Total wall time      : %s\n',  format_duration(elapsed_total));
logf('[ODE] Integration segments : %d\n',     segment_count);
logf('[ODE] Basis rotations      : %d\n',     basis_changes);
logf('[ODE] Final t*             : %.6e (target: %.6e)\n', ...
    T_achieved, T_obs_dimless);
logf('[ODE] Quaternion max drift : %.2e\n', ...
    max(abs(sqrt(sum(results.q.^2, 1)) - 1)));
if ~isempty(seg_walltimes)
    logf('[ODE] Segment timing stats : mean = %s, min = %s, max = %s\n', ...
        format_duration(mean(seg_walltimes)), ...
        format_duration(min(seg_walltimes)), ...
        format_duration(max(seg_walltimes)));
end
logf('================================================================\n\n');

%% =====================================================================
%  POST-PROCESSING: Convert to Physical Units
%  =====================================================================
%  Output unit convention (all physical fields in results):
%    Time [s] | Length [um] | Force [pN] | Angle [rad]
%    Velocity [um/s] | Angular velocity [rad/s]
%  =====================================================================

logf('[POST] Converting results to physical units...\n');

L_0 = params.scales.L_0;    % Length scale [m]
T_0 = params.scales.T_0;    % Time scale   [s]
U_0 = params.scales.U_0;    % Velocity scale [m/s]
F_0 = params.scales.F_0;    % Force scale  [N]

L_to_um = L_0 * 1e6;        % Conversion: dimensionless length -> um

% --- Time [s] ---
results.t = results.t_dimless * T_0;

% --- Body centre-of-mass position [um] ---
results.X_cm_um = results.X_cm * L_to_um;

% --- Pili node coordinates [um] ---
results.X_pili_um = cell(M_pili, 1);
for k = 1:M_pili
    results.X_pili_um{k} = results.X_pili{k} * L_to_um;
end

% --- Sinking force [pN] ---
results.F_sink_pN = params.Eg * F_0 * 1e12;

% --- Unit metadata ---
results.units.time     = 's';
results.units.length   = 'um';
results.units.force    = 'pN';
results.units.angle    = 'rad';
results.units.velocity = 'um/s';
results.units.ang_vel  = 'rad/s';
results.units.L_0_m    = L_0;
results.units.T_0_s    = T_0;
results.units.U_0_m_s  = U_0;
results.units.F_0_N    = F_0;

% --- Translational velocity [um/s] ---
% A timeout may occur before the first requested output interval.  Preserve
% that single state and report an undefined velocity instead of calling
% gradient with insufficient samples.
if NT >= 2
    V_sink_all = gradient(results.X_cm_um(3,:), results.t);
    results.V_cm_um = [ gradient(results.X_cm_um(1,:), results.t); ...
                        gradient(results.X_cm_um(2,:), results.t); ...
                        gradient(results.X_cm_um(3,:), results.t) ];
else
    V_sink_all = NaN(1, NT);
    results.V_cm_um = NaN(3, NT);
end
results.V_sink = V_sink_all;

% --- Body tilt angle [rad] ---
%  Elevation angle of the body long axis (body-local X) from horizontal,
%  consistent with the alpha_0 convention in init_coupled_params:
%    = +pi/2  body points +Z (vertical up)
%    =  0     body is horizontal
%    = -pi/2  body points -Z (vertical down)
tilt_angle_vec = zeros(1, NT);
for ti_pp = 1:NT
    R_ti = quat_to_rotmat(results.q(:, ti_pp));
    body_long_axis_z       = R_ti(3, 1);            % Z-component of body X-axis
    tilt_angle_vec(ti_pp)  = asin(max(-1, min(1, body_long_axis_z)));
end
results.tilt_angle = tilt_angle_vec;                % [rad], [1 x NT]

logf('[POST] Physical unit conversion complete.\n');
logf('[POST] Sinking force     : %.4f pN\n',   results.F_sink_pN);
logf('[POST] Final Z position  : %.4e um\n',    results.X_cm_um(3, end));
logf('[POST] Mean V_sink       : %.4e um/s\n',  mean(V_sink_all));
logf('[POST] Steady V_sink     : %.4e um/s (last 5%%)\n', ...
    mean(V_sink_all(max(1, round(0.95*NT)):end)));
logf('[POST] Initial tilt      : %.4f rad (%.2f deg)\n', ...
    tilt_angle_vec(1), rad2deg(tilt_angle_vec(1)));
logf('[POST] Final tilt        : %.4f rad (%.2f deg)\n', ...
    tilt_angle_vec(end), rad2deg(tilt_angle_vec(end)));
logf('\n');

%% =====================================================================
%  SAVE DATA (.mat)
%  =====================================================================

logf('[SAVE] Saving simulation data...\n');

n_flow = results.flow_state.n_times;
flow_len_ok = n_flow == numel(results.t) && n_flow == numel(results.t_dimless) && ...
    size(results.flow_state.Z, 2) == n_flow && ...
    size(results.flow_state.brown_force, 2) == n_flow && ...
    size(results.flow_state.Z, 1) == results.flow_state.state_dimension;
for k = 1:M_pili
    flow_len_ok = flow_len_ok && size(results.flow_state.rot{k}, 3) == n_flow;
end
if ~flow_len_ok
    error('run_simulation_HPC_Linux:FlowStateLength', ...
        'flow_state history is inconsistent with the saved time vector.');
end
logf('[SAVE] Flow-state history: version %d | %d states x %d outputs | force convention %s\n', ...
    results.flow_state.version, results.flow_state.state_dimension, n_flow, ...
    results.flow_state.brownian_endpoint_convention);
logf('[SAVE] Flow-state final sample vs last_Z max |diff| = %.3e\n', ...
    max(abs(results.flow_state.Z(:, end) - last_Z(:))));

% Remove function handles before saving (not serialisable to .mat)
params = rmfield(params, {'logf', 'logf_file'});

mat_file = fullfile(output_dir, 'simulation_results.mat');
save(mat_file, 'results', 'params', '-v7.3');
mat_info = dir(mat_file);
logf('[SAVE] Data saved  : %s (%.2f MB)\n', mat_file, mat_info.bytes/1e6);

% --- Export scalar time series as CSV ---
csv_file = fullfile(output_dir, 'timeseries.csv');
T_csv = table( ...
    results.t', ...
    results.X_cm_um(1,:)', results.X_cm_um(2,:)', results.X_cm_um(3,:)', ...
    results.V_cm_um(1,:)', results.V_cm_um(2,:)', results.V_cm_um(3,:)', ...
    results.V_sink', ...
    results.tilt_angle', rad2deg(results.tilt_angle)', ...
    'VariableNames', {'t_s', 'X_um', 'Y_um', 'Z_um', ...
                      'Vx_um_s', 'Vy_um_s', 'Vz_um_s', ...
                      'Vsink_um_s', 'tilt_rad', 'tilt_deg'});
writetable(T_csv, csv_file);
csv_info = dir(csv_file);
logf('[SAVE] CSV exported: %s (%.1f KB)\n', csv_file, csv_info.bytes/1e3);

%% =====================================================================
%  PRE-COMPUTE BODY SURFACE MESH (for visualisation)
%  =====================================================================
%  For parametric mode: Generate a smooth parametric capsule mesh.
%  For custom mode: Use the triangulated mesh from custom_body_building.
%  Units: micrometres [um].
%  =====================================================================

logf('[PLOT] Building body surface mesh for visualisation...\n');

if strcmp(params.geometry_mode, 'parametric') && strcmp(params.body_type, 'capsule')
    % --- Parametric capsule surface mesh ---
    a_um = major_axis * 1e6;       % Semi-major axis [um]
    b_um = minor_axis * 1e6;       % Semi-minor axis (cap radius) [um]

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
    mesh_size   = size(capsule_X0);

    % Pre-flatten for fast rotation: [3 x N_mesh]
    capsule_flat = [capsule_X0(:)'; capsule_Y0(:)'; capsule_Z0(:)'];

    body_mesh_type = 'parametric_surf';
    logf('[PLOT] Capsule mesh: %d x %d = %d vertices\n', ...
        mesh_size(1), mesh_size(2), numel(capsule_X0));

else
    % --- Custom body triangulated mesh ---
    body_mesh_type = 'custom_trisurf';

    % Convert vertices from [m] to [um]
    custom_verts_um = params.body_verts * 1e6;

    logf('[PLOT] Custom mesh: %d vertices, %d faces\n', ...
        size(params.body_verts, 1), size(params.body_faces, 1));
end

% --- Colour palette (journal style) ---

% Per-pili distinguishable colours (up to 10; cycles if M_pili > 10)
pili_palette = [0.20 0.40 0.80;      % Medium blue
                0.85 0.33 0.10;      % Burnt orange
                0.47 0.67 0.19;      % Olive green
                0.64 0.08 0.18;      % Crimson
                0.49 0.18 0.56;      % Purple
                0.93 0.69 0.13;      % Goldenrod
                0.30 0.75 0.93;      % Sky blue
                0.60 0.60 0.60;      % Grey
                0.00 0.45 0.74;      % Dark blue
                0.85 0.47 0.66];     % Rose pink
pili_colors = pili_palette(mod((1:M_pili)-1, size(pili_palette,1)) + 1, :);

%% =====================================================================
%  STATIC VISUALISATION (saved to PNG)
%  =====================================================================

logf('[PLOT] Generating summary figure...\n');

fig_summary_vis = 'on';
if silent
    fig_summary_vis = 'off';
end
fig_summary = figure('Position', [100, 100, 1600, 1000], 'Color', 'w', ...
                     'Visible', fig_summary_vis);
set(fig_summary, 'DefaultAxesFontName',  'Arial', ...
                 'DefaultAxesFontSize',  11, ...
                 'DefaultAxesLineWidth', 0.8, ...
                 'DefaultAxesTickDir',   'out', ...
                 'DefaultLineLineWidth', 1.8);

t_s  = results.t;
X_um = results.X_cm_um;

% --- Subplot 1: Body centre-of-mass trajectory ---
subplot(2, 3, 1);
plot(t_s, X_um(1,:), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8); hold on;
plot(t_s, X_um(2,:), '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.8);
plot(t_s, X_um(3,:), '-', 'Color', [0.20 0.40 0.80], 'LineWidth', 1.8);
xlabel('Time (s)',                 'FontName', 'Arial', 'FontSize', 12);
ylabel('Position (μm)',        'FontName', 'Arial', 'FontSize', 12);
title('Body Centre of Mass',      'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
lg = legend('X', 'Y', 'Z', 'Location', 'best');
set(lg, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off');
grid on; box on;

% --- Subplot 2: Sinking velocity ---
subplot(2, 3, 2);
plot(t_s, results.V_sink, 'k-', 'LineWidth', 1.8);
xlabel('Time (s)',                        'FontName', 'Arial', 'FontSize', 12);
ylabel('{\itV}_{sink} (μm/s)',        'FontName', 'Arial', 'FontSize', 12);
title('Sinking Velocity',                 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
grid on; box on;

% --- Subplot 3: Body tilt angle ---
subplot(2, 3, 3);
plot(t_s, rad2deg(results.tilt_angle), '-', ...
     'Color', [0.49 0.18 0.56], 'LineWidth', 1.8);
xlabel('Time (s)',                'FontName', 'Arial', 'FontSize', 12);
ylabel('Tilt angle \alpha (°)',   'FontName', 'Arial', 'FontSize', 12);
title('Body Tilt Angle',          'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
grid on; box on;

% --- Subplot 4-6: 3D/XZ/YZ snapshot views ---
ax_snap_3d = subplot(2, 3, 4);
ax_snap_xz = subplot(2, 3, 5);
ax_snap_yz = subplot(2, 3, 6);
snap_axes  = [ax_snap_3d, ax_snap_xz, ax_snap_yz];
for ai = 1:numel(snap_axes)
    hold(snap_axes(ai), 'on');
    grid(snap_axes(ai), 'on');
    box(snap_axes(ai),  'on');
end

n_snapshots  = min(10, NT);
snap_indices = round(linspace(1, NT, n_snapshots));

% Time-gradient colormaps: light (early) -dark (late)
body_cmap = [linspace(0.75, 0.10, n_snapshots)', ...   % R
             linspace(0.92, 0.55, n_snapshots)', ...   % G
             linspace(0.78, 0.25, n_snapshots)'];      % B  (light→dark green)

% Per-pili time-gradient colormaps: light (early) -dark (late)
pili_cmaps = cell(M_pili, 1);
for k = 1:M_pili
    base_c = pili_colors(k, :);
    % Blend from light (0.7*white + 0.3*base) to dark (base * 0.65)
    light_c = 0.70 + 0.30 * base_c;
    dark_c  = base_c * 0.65;
    pili_cmaps{k} = [linspace(light_c(1), dark_c(1), n_snapshots)', ...
                     linspace(light_c(2), dark_c(2), n_snapshots)', ...
                     linspace(light_c(3), dark_c(3), n_snapshots)'];
end

% Fixed axis limits across all static snapshot views
all_x = [];
all_y = [];
all_z = [];
for si = 1:n_snapshots
    ti = snap_indices(si);
    if ti <= size(results.X_cm, 2)
        R_snap = quat_to_rotmat(results.q(:, ti));
        cm_um  = results.X_cm_um(:, ti);
        if strcmp(body_mesh_type, 'parametric_surf')
            pts_snap = R_snap * capsule_flat + cm_um;    % [3 x N_mesh]
        else
            pts_snap = R_snap * custom_verts_um' + cm_um; % [3 x N_verts]
        end
        all_x = [all_x; pts_snap(1,:)']; %#ok<AGROW>
        all_y = [all_y; pts_snap(2,:)']; %#ok<AGROW>
        all_z = [all_z; pts_snap(3,:)']; %#ok<AGROW>
        for k = 1:M_pili
            X_p = results.X_pili_um{k}(:,:,ti);
            all_x = [all_x; X_p(:,1)]; %#ok<AGROW>
            all_y = [all_y; X_p(:,2)]; %#ok<AGROW>
            all_z = [all_z; X_p(:,3)]; %#ok<AGROW>
        end
    end
end
x_span = max(all_x) - min(all_x);
y_span = max(all_y) - min(all_y);
z_span = max(all_z) - min(all_z);
margin_x = max(0.10 * x_span, 0.15);
margin_y = max(0.10 * y_span, 0.15);
margin_z = max(0.10 * z_span, 0.15);
ax_lim_x = [min(all_x) - margin_x, max(all_x) + margin_x];
ax_lim_y = [min(all_y) - margin_y, max(all_y) + margin_y];
ax_lim_z = [min(all_z) - margin_z, max(all_z) + margin_z];

% Enforce minimum axis range
if strcmp(body_mesh_type, 'parametric_surf')
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

for si = 1:n_snapshots
    ti = snap_indices(si);
    if ti <= size(results.X_cm, 2)
        R_snap = quat_to_rotmat(results.q(:, ti));
        cm_um  = results.X_cm_um(:, ti);

        snap_alpha = 0.15 + 0.30 * (si / n_snapshots);

        % Body surface rendering (dispatch by mesh type)
        if strcmp(body_mesh_type, 'parametric_surf')
            pts_snap = R_snap * capsule_flat + cm_um;
            Xs = reshape(pts_snap(1,:), mesh_size);
            Ys = reshape(pts_snap(2,:), mesh_size);
            Zs = reshape(pts_snap(3,:), mesh_size);
            for ai = 1:numel(snap_axes)
                surf(snap_axes(ai), Xs, Ys, Zs, ...
                     'FaceColor', body_cmap(si,:), 'FaceAlpha', snap_alpha, ...
                     'EdgeColor', 'none');
            end
        else
            % Custom body: use trisurf with auto-generated triangulation.
            % NOTE: unlike surf(ax, ...), MATLAB's trisurf(tri, varargin)
            % binds its FIRST positional argument to the face matrix and
            % never checks it for an axes handle (axescheck only scans
            % varargin{2:end}), so trisurf(ax, F, x, y, z, ...) silently
            % misinterprets F as x/y/z data and crashes with
            % MATLAB:catenate:dimensionMismatch.  Target the axes via the
            % 'Parent' name-value pair instead, which trisurf handles
            % correctly (verified against R2026a).
            rotated_verts = (R_snap * custom_verts_um' + cm_um)';  % [V x 3]
            for ai = 1:numel(snap_axes)
                trisurf(params.body_faces, ...
                    rotated_verts(:,1), rotated_verts(:,2), rotated_verts(:,3), ...
                    'FaceColor', body_cmap(si,:), 'FaceAlpha', snap_alpha, ...
                    'EdgeColor', 'none', 'Parent', snap_axes(ai));
            end
        end

        % Pili polylines (per-pili gradient colours)
        for k = 1:M_pili
            X_p = results.X_pili_um{k}(:,:,ti);
            pk_c = pili_cmaps{k}(si,:);
            for ai = 1:numel(snap_axes)
                plot3(snap_axes(ai), X_p(:,1), X_p(:,2), X_p(:,3), ...
                      '-', 'Color', pk_c, 'LineWidth', 0.9);
            end
        end
    end
end

% Common formatting for all three snapshot views
for ai = 1:numel(snap_axes)
    axi = snap_axes(ai);
    xlim(axi, ax_lim_x);
    ylim(axi, ax_lim_y);
    zlim(axi, ax_lim_z);
    set(axi, 'DataAspectRatio', [1 1 1], 'TickDir', 'out');
    % Use uniform ambient lighting (no directional light source) so that
    % surfaces display their true FaceColor without Gouraud darkening.
    % A single headlight dims capsule faces whose normals are near-
    % perpendicular to the camera direction, making the body appear dark
    % regardless of the chosen colour gradient.  lighting('none') renders
    % pure FaceColor at the specified alpha, giving the intended bright
    % green/blue gradient appearance.
    lighting(axi, 'none');
    xlabel(axi, 'X (μm)', 'FontName', 'Arial', 'FontSize', 12);
    ylabel(axi, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 12);
    zlabel(axi, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 12);

end

view(ax_snap_3d, 45, 30);
view(ax_snap_xz, 0, 0);      % XZ plane: azimuth=0, elevation=0 -> X right, Z up
view(ax_snap_yz, 90, 0);     % YZ plane: azimuth=90, elevation=0 -> Y right, Z up
title(ax_snap_3d, '3D Body-Pili Snapshots', ...
      'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
title(ax_snap_xz, 'XZ View Snapshots', ...
      'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
title(ax_snap_yz, 'YZ View Snapshots', ...
      'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');

sgtitle(sprintf('Bacteria Elastohydrodynamic Simulation  [%s]', run_id), ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
drawnow;   % Flush sgtitle layout BEFORE measuring axis positions

% Gravity indicators: must be added AFTER sgtitle so that sgtitle's
% layout adjustment (which shifts all subplot Positions inward) has
% already been applied when we read ax.Position.
for ai = 1:numel(snap_axes)
    add_gravity_indicator_outside(fig_summary, snap_axes(ai), 12, 2.2);
end

% Save figure (raster + vector + HTML on R2026a+)
png_file  = fullfile(output_dir, 'summary_plot.png');
pdf_file  = fullfile(output_dir, 'summary_plot.pdf');
html_file = fullfile(output_dir, 'summary_plot.html');
export_figure_png_pdf_html(fig_summary, png_file, pdf_file, html_file);
close(fig_summary);
logf('[PLOT] Summary figure saved: %s  |  %s\n', png_file, pdf_file);
if isfile(html_file)
    logf('[PLOT] Interactive HTML: %s\n', html_file);
end

%% =====================================================================
%  3D ANIMATION (saved to MP4)
%  =====================================================================

if generate_animation
    anim_opts = struct('fps', animation_fps, 'skip', animation_skip, ...
                       'silent', silent, 'logf', logf, 'logf_file', logf_file);
    anim_info = render_motion_animation(results, params, output_dir, anim_opts);
    t_anim_elapsed = anim_info.elapsed_s;
    anim_filename  = anim_info.filename;
elseif timed_out && requested_generate_animation
    t_anim_elapsed = 0;
    logf('\n[ANIM] Skipped (integration timeout).\n\n');
end

%% =====================================================================
%  FLOW FIELD VISUALISATION (XZ Plane at Final State)
%  =====================================================================

if compute_flowfield
    logf('\n[FLOW] Computing XZ-plane flow field at final state...\n');
    t_flow_start = tic;
    flow_opts.grid_spacing_um = flowfield_grid_spacing;
    flow_opts.grid_range_um   = flowfield_grid_range;
    flow_opts.quiver_skip     = flowfield_quiver_skip;
    flow_opts.chunk_size      = flowfield_chunk_size;
    params.brown_force_dimless = results.flow_state.brown_force(:, end);
    compute_plot_flowfield(last_Z, params, output_dir, run_id, logf, flow_opts, T_achieved, silent);
    t_flow_elapsed = toc(t_flow_start);
    logf('[FLOW] Flow field visualisation completed in %s\n\n', format_duration(t_flow_elapsed));
elseif timed_out && requested_compute_flowfield
    t_flow_elapsed = 0;
    logf('\n[FLOW] Skipped (integration timeout).\n\n');
else
    t_flow_elapsed = 0;
    logf('\n[FLOW] Skipped (cfg.compute_flowfield = false).\n\n');
end

%% =====================================================================
%  FINAL SUMMARY
%  =====================================================================

total_sim_walltime = toc(tic_sim_start);

logf('\n================================================================\n');
if timed_out
    logf('  SIMULATION PARTIAL: INTEGRATION TIMEOUT\n');
else
    logf('  SIMULATION COMPLETE\n');
end
logf('================================================================\n');
logf('  Run ID          : %s\n',    run_id);
logf('  Timestamp       : %s\n',    run_timestamp);
logf('  Output directory : %s\n',    output_dir);
logf('  Physical time    : %.2f s\n', T_obs);
logf('  Sinking force    : %.4f pN\n', results.F_sink_pN);
logf('  Steady V_sink    : %.4e um/s\n', ...
    mean(V_sink_all(max(1, round(0.95*NT)):end)));
logf('  Basis rotations  : %d\n',    basis_changes);
logf('  Max |q| drift    : %.2e\n',  max(abs(sqrt(sum(results.q.^2,1))-1)));
logf('  ------------------------------------------------\n');
logf('  Timing breakdown:\n');
logf('    Initialisation  : %s\n',   format_duration(t_init_elapsed));
logf('    ODE integration : %s\n',   format_duration(elapsed_total));
if generate_animation
    logf('    Animation render: %s\n',   format_duration(t_anim_elapsed));
end
logf('    Flow field      : %s\n',   format_duration(t_flow_elapsed));
logf('    Total wall time : %s\n',   format_duration(total_sim_walltime));
logf('  ------------------------------------------------\n');
logf('  Files saved:\n');
logf('    Data   : simulation_results.mat\n');
logf('    CSV    : timeseries.csv\n');
logf('    Config : config.json\n');
if enable_brown
    logf('    Seed   : brownian_rng_seed.txt\n');
end
logf('    Plot   : summary_plot.png, summary_plot.pdf');
if isfile(fullfile(output_dir, 'summary_plot.html'))
    logf(', summary_plot.html\n');
else
    logf('\n');
end
if preview_initial
    logf('    Pose   : initial_pose.png, initial_pose.pdf');
    if isfile(fullfile(output_dir, 'initial_pose.html'))
        logf(', initial_pose.html\n');
    else
        logf('\n');
    end
end
if generate_animation
    logf('    Video  : %s\n', anim_filename);
end
if compute_flowfield
    logf('    Flow   : flowfield_panels.png, flowfield_panels.pdf');
    if isfile(fullfile(output_dir, 'flowfield_panels.html'))
        logf(', flowfield_panels.html\n');
    else
        logf('\n');
    end
else
    logf('    Flow   : (skipped)\n');
end
logf('    Log    : simulation.log\n');
logf('  ------------------------------------------------\n');
logf('  results fields added:\n');
logf('    tilt_angle  [rad] : body long-axis elevation angle vs time\n');
logf('    V_cm_um [um/s x3] : full translational velocity vs time\n');
logf('================================================================\n');
logf('  Output units: time [s], length [um], force [pN], angle [rad]\n');
logf('================================================================\n\n');

% Close log file
fclose(log_fid);
if ~silent
    fprintf('Log saved: %s\n', log_file);
end

catch ME_run
    % Any unhandled error above lands here: record full diagnostics in
    % simulation.log (which would otherwise just truncate silently),
    % release any figure/video handles opened before the error, then
    % rethrow so the caller (interactive script or batch driver) still
    % sees the original exception unchanged.
    err_report = getReport(ME_run, 'extended', 'hyperlinks', 'off');
    fatal_msg = sprintf(['\n================================================================\n' ...
        '[FATAL] Run terminated by an unhandled error\n' ...
        '================================================================\n' ...
        '[FATAL] Identifier: %s\n' ...
        '[FATAL] Message   : %s\n' ...
        '[FATAL] %s\n' ...
        '================================================================\n'], ...
        ME_run.identifier, ME_run.message, err_report);

    if ~isempty(fopen(log_fid))
        fprintf(log_fid, '%s', fatal_msg);
    end
    if ~silent
        fprintf('%s', fatal_msg);
    end

    % Best-effort release of figure/video handles opened before the
    % error; wrapped in its own try/catch so a cleanup failure never
    % masks the original error being rethrown below.
    try
        if exist('fig_preview', 'var') && isgraphics(fig_preview), close(fig_preview); end
        if exist('fig_summary', 'var') && isgraphics(fig_summary), close(fig_summary); end
    catch
        % Ignore cleanup failures; the original error is what matters.
    end

    % log_fid is also closed by cleanup_log_fid's onCleanup as a backstop,
    % but close it here too so it is released immediately (not just at
    % function-exit) in case downstream code inspects the file right away.
    if ~isempty(fopen(log_fid))
        fclose(log_fid);
    end

    rethrow(ME_run);
end

end  % run_simulation_HPC_Linux


%% ==================== Local Functions ====================

function safe_fclose(fid)
%SAFE_FCLOSE  Close a file identifier if it is still open; no-op/no-error
%  otherwise.  Used as an onCleanup target to guarantee simulation.log is
%  always released (normal return, error, or Ctrl-C), even though the
%  normal and catch paths above already close it explicitly.
    if isnumeric(fid) && isscalar(fid) && fid > 0 && ~isempty(fopen(fid))
        fclose(fid);
    end
end

function v = get_cfg_field(cfg, fname, default)
%GET_CFG_FIELD  Return cfg.(fname) if present and non-empty, else default
    if isstruct(cfg) && isfield(cfg, fname) && ~isempty(cfg.(fname))
        v = cfg.(fname);
    else
        v = default;
    end
end

function dual_log(fid, varargin)
%DUAL_LOG  Write formatted message to console and log file (with timestamp).
%
%  When the trailing logical flag echo_console is false (silent mode, used
%  by batch parallel runs), the console write is suppressed so that N
%  workers do not interleave their stdout.  The log file always receives
%  the message.
%
%  Trailing-arg convention: varargin = {msg_fmt, arg1, ..., echo_console}
%  where echo_console defaults to true when absent (preserves legacy calls).

    if nargin >= 2 && islogical(varargin{end})
        echo_console = varargin{end};
        varargin(end) = [];
    else
        echo_console = true;
    end
    msg = sprintf(varargin{:});
    ts  = char(datetime('now', 'Format', 'uuuu-MM-dd HH:mm:ss'));
    if echo_console
        fprintf('%s', msg);                % Console (no timestamp - cleaner)
    end
    fprintf(fid, '[%s] %s', ts, msg);     % File (with timestamp)
end

function log_file_only(fid, varargin)
%LOG_FILE_ONLY  Write formatted message to log file only (with timestamp)
    msg = sprintf(varargin{:});
    ts  = char(datetime('now', 'Format', 'uuuu-MM-dd HH:mm:ss'));
    fprintf(fid, '[%s] %s', ts, msg);
end

function s = fmt_field_value(S, fname)
%FMT_FIELD_VALUE  Compact string for a struct field used in [PARAMS] logging
    if ~isstruct(S) || ~isfield(S, fname)
        s = '(n/a)';
        return;
    end
    v = S.(fname);
    if isempty(v)
        s = '(n/a)';
    elseif isnumeric(v) && isscalar(v)
        s = sprintf('%.4g', v);
    elseif isnumeric(v)
        s = mat2str(v, 4);
    elseif ischar(v) || isstring(v)
        s = char(v);
    else
        s = '(complex)';
    end
end

function v = get_field_default(S, fname, default_val)
%GET_FIELD_DEFAULT  Return S.(fname) if present and non-empty, else default_val
    if isstruct(S) && isfield(S, fname) && ~isempty(S.(fname))
        v = S.(fname);
    else
        v = default_val;
    end
end

function str = format_duration(seconds)
%FORMAT_DURATION  Convert elapsed seconds to human-readable string
%  Examples:
%    format_duration(3.7)    -> '3.7 s'
%    format_duration(125.3)  -> '2 min 5.3 s'
%    format_duration(3725.1) -> '1 h 2 min 5.1 s'
    if seconds < 0
        str = '-- s';
        return;
    end
    h = floor(seconds / 3600);
    m = floor(mod(seconds, 3600) / 60);
    s = mod(seconds, 60);
    if h > 0
        str = sprintf('%d h %d min %.1f s', h, m, s);
    elseif m > 0
        str = sprintf('%d min %.1f s', m, s);
    else
        str = sprintf('%.1f s', s);
    end
end

function status = timeout_output(t, y, flag, tstart, tlim, silent)
%TIMEOUT_OUTPUT  Stop ode15s when the segment wall-clock limit is reached
    if toc(tstart) >= tlim
        status = 1;
    elseif silent
        status = 0;
    else
        status = odetpbar(t, y, flag);
    end
end


function [value, isterminal, direction] = singularity_event(~, Z, params, tstart, tlim)
%SINGULARITY_EVENT  ODE event: Euler angle singularity or timeout
    delta  = params.delta;
    N_segs = params.N_segments;
    M      = params.M_pili;

    % Check timeout
    elapsed     = toc(tstart);
    % Clamp at zero so a limit already exceeded before the first event
    % evaluation is still reported as an initial terminal event.  Returning
    % an always-negative value would provide no sign crossing to ode15s.
    timeout_val = max(tlim - elapsed, 0);

    % Check theta singularities
    min_margin = inf;
    [~, ~, pili_angles] = unpack_state(Z, N_segs);
    for kk = 1:M
        theta_kk  = pili_angles{kk}.theta;
        margin_kk = min([theta_kk; pi - theta_kk]);
        min_margin = min(min_margin, margin_kk);
    end

    singularity_val = min_margin - delta/2;

    % Return separate event channels so solver reports exact trigger reason:
    %   value(1) -> timeout
    %   value(2) -> Euler-angle singularity margin
    value      = [timeout_val; singularity_val];
    isterminal = [1; 1];
    direction  = [-1; -1];
end
