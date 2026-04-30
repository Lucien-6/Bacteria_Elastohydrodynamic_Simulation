function [results, params, output_dir] = run_simulation(cfg_or_path)
%RUN_SIMULATION  Execute coupled body-pili simulation with given configuration
%
%  Description:
%    Runs the full coupled body-pili elastohydrodynamic simulation:
%    initialisation, ODE integration, post-processing, visualisation,
%    and data export.  Encapsulates the logic formerly in the main script
%    so it can be called programmatically (e.g. for parameter sweeps).
%
%  Syntax:
%    [results, params, output_dir] = run_simulation(cfg)
%    [results, params, output_dir] = run_simulation(json_path)
%
%  Inputs:
%    cfg_or_path - Configuration struct (from load_sim_config or inline),
%                  OR a char/string path to a JSON config file
%
%  Outputs:
%    results    - Simulation results struct (physical units); [] if aborted
%    params     - Parameter struct (including scales); [] if aborted
%    output_dir - Path to output directory; '' if aborted
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17
%  Log banner version: 7.0.4 (2026-04-30) — keep in sync with README / CHANGELOG
%
%  Path: defensive addpath(project root) so run_simulation(cfg) works when
%  invoked without main_coupled_simulation.m. Unified policy = this line plus
%  addpath in main_coupled_simulation; post_processing/ and batch sweep/
%  retain separate addpath for standalone use.

addpath(fileparts(mfilename('fullpath')));

%% ---- Input parsing ----
if ischar(cfg_or_path) || isstring(cfg_or_path)
    cfg = load_sim_config(cfg_or_path);
else
    cfg = cfg_or_path;
end

results    = [];
params     = [];
output_dir = '';

%% ---- Unpack cfg to local variables ----
case_name              = cfg.case_name;
geometry_mode          = get_cfg_field(cfg, 'geometry_mode', 'parametric');

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
if isfield(cfg, 'fixman_midpoint')
    fixman_midpoint = logical(cfg.fixman_midpoint);
else
    fixman_midpoint = false;   % backward-compatible with cfg structs predating this field
end
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

% ---- Brownian RNG (twister): must run before save_sim_config and ODE integration ----
brown_rng_seed_used = NaN;
if enable_brown
    if isfield(cfg, 'brown_rng_seed') && isscalar(cfg.brown_rng_seed) ...
            && isnumeric(cfg.brown_rng_seed) && isfinite(cfg.brown_rng_seed) ...
            && cfg.brown_rng_seed > 0
        brown_rng_seed_used = round(double(cfg.brown_rng_seed));
    else
        brown_rng_seed_used = brown_rng_seed_from_clock();
    end
    rng(brown_rng_seed_used, 'twister');
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
    data_exts = {'.mat', '.csv', '.png', '.pdf', '.html', '.mp4', '.log', '.json', '.txt'};
    has_data  = any(cellfun(@(ext) any(endsWith({existing_files.name}, ext)), data_exts));

    if has_data
        if auto_confirm
            output_dir = [output_dir '-' run_timestamp];
            run_id     = [run_id '-' run_timestamp];
            mkdir(output_dir);
            fprintf('[OUTPUT] Auto-renamed to avoid overwrite: %s\n', output_dir);
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
                    fprintf('[OUTPUT] Cleared %d existing file(s) and proceeding.\n', ...
                        numel(existing_files));
                case 'Rename (append timestamp)'
                    output_dir = [output_dir '-' run_timestamp];
                    run_id     = [run_id '-' run_timestamp];
                    mkdir(output_dir);
                    fprintf('[OUTPUT] Renamed to: %s\n', output_dir);
                otherwise
                    fprintf('Simulation aborted: output directory conflict.\n');
                    return;
            end
        end
    end
else
    mkdir(output_dir);
end

% --- Initialise dual-channel logger (console + file) ---
log_file = fullfile(output_dir, 'simulation.log');
log_fid  = fopen(log_file, 'w');
if log_fid == -1
    error('Cannot create log file: %s', log_file);
end
logf = @(varargin) dual_log(log_fid, varargin{:});  %#ok<NFFT>  % persistent closure

% Log-file-only variant (no console echo) for high-frequency diagnostics
logf_file = @(varargin) log_file_only(log_fid, varargin{:});

logf('\n');
logf('================================================================\n');
logf('  Coupled Body-Pili Elastohydrodynamic Simulation\n');
logf('  Version 7.0.4\n');
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
    fid_seed = fopen(seed_txt, 'w');
    if fid_seed ~= -1
        fprintf(fid_seed, '%d\n', brown_rng_seed_used);
        fprintf(fid_seed, ...
            'Twister seed used for Brownian kicks. Re-run with the same seed via cfg.brown_rng_seed or JSON brownian.rng_seed.\n');
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
                    freq_str = fmt_field_value(am_kk, 'frequency_Hz');
                    wvl_str  = fmt_field_value(am_kk, 'wavelength');
                    logf(['[PARAMS]     active_moment: traveling_wave | ', ...
                          'amp = %s Nm, freq = %s Hz, wavelength = %s\n'], ...
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
logf('[PARAMS] Fluid: mu = %.2e Pa*s, rho_bac = %.1f kg/m^3, rho_fluid = %.3f kg/m^3\n', ...
    mu, rho_bac, rho_fluid);
logf('[PARAMS] Simulation: T_obs = %.1f s, dt = %.2e s, NT = %d, AbsTol = %.0e, RelTol = %.0e\n', ...
    T_obs, dt, NT, AbsTol, RelTol);
logf('[PARAMS] Numerics: delta = pi/%.0f rad, timeout = %d s\n', pi/delta, tlim);
logf('[PARAMS] Preview: %s\n', mat2str(preview_initial));
logf('[PARAMS] Animation: fps = %d, skip = %d\n', animation_fps, animation_skip);
logf('[PARAMS] Flowfield: grid = %.2f um, range = %.1f um, quiver_skip = %d\n', ...
    flowfield_grid_spacing, flowfield_grid_range, flowfield_quiver_skip);
logf('[PARAMS] Brownian: enable = %s, Temper = %.2f K, dt_B = %.4f s, fixman = %s\n', ...
    mat2str(enable_brown), Temper, dt_B, mat2str(fixman_midpoint));
if enable_brown
    logf('[PARAMS] Brownian RNG: seed = %d (twister), saved for reproducibility\n', ...
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
val_cfg = struct( ...
    'T_obs', T_obs, 'dt', dt, 'minor_axis', minor_axis, ...
    'major_axis', major_axis, 'N_body', N_body, 'seg_per_um', seg_per_um, ...
    'Pili_Matrix', Pili_Matrix, 'E', E, 'r_pili', r_pili, 'mu', mu, ...
    'rho_bac', rho_bac, 'rho_fluid', rho_fluid, ...
    'AbsTol', AbsTol, 'RelTol', RelTol, 'tlim', tlim, ...
    'enable_brown', enable_brown, 'dt_B', dt_B, ...
    'animation_fps', animation_fps, 'sigma_p', sigma_p, ...
    'fixman_midpoint', fixman_midpoint, ...
    'geometry_mode', geometry_mode);
if isfield(cfg, 'background_flow')
    val_cfg.background_flow = cfg.background_flow;
end
if isfield(cfg, 'pili_ext')
    val_cfg.pili_ext = cfg.pili_ext;
end
if strcmp(geometry_mode, 'custom') && isfield(cfg, 'custom_geometry')
    val_cfg.custom_geometry = cfg.custom_geometry;
end
validate_simulation_params(val_cfg, logf);
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
    'mu',           mu, ...
    'E',            E, ...
    'r_pili',       r_pili, ...
    'sigma',        sigma_p, ...
    'rho_bac',      rho_bac, ...
    'rho_fluid',    rho_fluid, ...
    'N_segments',   N_segments, ...
    'Pili_Matrix',  Pili_Matrix, ...
    'seg_per_um',   seg_per_um, ...
    'delta',        delta, ...
    'Temper',       Temper, ...
    'enable_brown', enable_brown, ...
    'fixman_midpoint', fixman_midpoint, ...
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
logf('[INIT]   epsilon = %.4e  (slenderness r_pili / L_0)\n',  params.epsilon);
if ~isempty(eps_body_cfg)
    logf('[INIT]   eps_body (RSM override, L0): %.6f  (h_mean = %.6f)\n', ...
        params.eps_body, params.body_h_mean);
end
logf('[INIT]   alpha_0 = %.4f rad (%.1f deg) — initial body tilt\n', ...
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
        logf('[INIT]          eps_pili = %.4f L0 (%.1f nm = r_pili), ds/eps = %.1f\n', ...
            params.epsilon, eps_nm, ds_nm / eps_nm);
    end
    dis_dim_nm = params.dim.dis * 1e9;
    logf('[INIT]   Gap : dis = 1.1*(eps_body + eps_pili) = %.1f nm\n', dis_dim_nm);
    logf('[INIT]   Note: eps_ij = eps_j (determined by force point j)\n');
    logf('[INIT]     Body force -> any point: eps = %.1f nm\n', ...
        params.eps_body * L_0_dim * 1e9);
    logf('[INIT]     Pili force -> any point: eps = %.1f nm (physical radius)\n', ...
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
        all_pts = [all_pts; X_pili_prev{kp}]; 
    end
    pv_min = min(all_pts);
    pv_max = max(all_pts);
    pv_span = pv_max - pv_min;
    pv_margin = max(0.10 * pv_span, 0.15);   % [1x3] per-axis margin

    % Create figure with 3 views
    fig_preview = figure('Name', 'Initial Pose Preview', ...
                         'NumberTitle', 'off', ...
                         'Position', [80, 80, 1500, 500], ...
                         'Color', 'w', 'Visible', 'on');
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
            legend_entries{end+1} = sprintf('Pili %d', kp); 
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
        logf('[PREVIEW] Auto-confirm enabled — skipping interactive dialog.\n');
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
        fprintf('Simulation aborted by user at initial pose preview.\n');
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

% --- Integration loop ---
clear odetpbar;                    % Reset persistent vars from prior runs
tic_total        = tic;
T_achieved       = 0;
last_Z           = Z0;
current_time_ind = 1;
basis_changes    = 0;
segment_count    = 0;
seg_walltimes    = [];          % Wall time per segment [s]
seg_dtstar       = [];          % Dimensionless time advanced per segment

% ---- Brownian dynamics node tracking ----
if params.brown_enabled
    dt_B_dimless   = dt_B / params.scales.T_0;
    T_brown_next   = dt_B_dimless;
    logf('[BROWN] Brownian dynamics enabled (full-system effective friction)\n');
    logf('[BROWN] dt_B = %.4f s = %.4e T*\n', dt_B, dt_B_dimless);
    if params.brown_fixman_midpoint
        logf('[BROWN] Fixman midpoint correction: ENABLED (2x mobility evaluation per step)\n');
    else
        logf('[BROWN] Fixman midpoint correction: DISABLED (single evaluation per step)\n');
    end
    brown_kick_count    = 0;
    brown_log_interval  = 100;
else
    dt_B_dimless = Inf;
    T_brown_next = Inf;
end

while T_achieved < T_obs_dimless
    Z = last_Z;
    segment_count = segment_count + 1;

    % Adaptive basis check for all pili
    [Z, params, rotated_flags] = adaptive_basis_check(Z, params);
    if any(rotated_flags)
        basis_changes = basis_changes + sum(rotated_flags);
        for kk = 1:length(rotated_flags)
            if rotated_flags(kk)
                logf('[ODE] Basis rotation triggered for pili %d at t* = %.4e\n', ...
                    kk, T_achieved);
            end
        end
    end

    % Determine integration target for this segment
    T_seg_target = min(T_brown_next, T_obs_dimless);

    % Setup ODE solver
    tstart_ode   = tic;
    eventFunc    = @(t_ev, Z_ev) singularity_event(t_ev, Z_ev, params, tstart_ode, tlim);
    progressFunc = @(t_p, y_p, fl_p, varargin) odetpbar(t_p, y_p, fl_p);

    ode_ops = odeset('OutputFcn',  progressFunc, ...
                     'Events',     eventFunc, ...
                     'Stats',      'off', ...
                     'AbsTol',     AbsTol, ...
                     'RelTol',     RelTol, ...
                     'Vectorized', 'off');

    dZ_func = @(t_ode, Z_ode) ode_coupled_rhs(t_ode, Z_ode, params);

    logf('[ODE] Segment %d: integrating t* = [%.4e, %.4e]\n', ...
        segment_count, T_achieved, T_seg_target);

    sol      = ode15s(dZ_func, [T_achieved, T_seg_target], Z, ode_ops);
    seg_wall = toc(tstart_ode);
    logf('[ODE] Segment %d completed (%.1f s, %d steps)\n', ...
        segment_count, seg_wall, length(sol.x));

    % Check if event was triggered
    if ~isempty(sol.xe)
        event_ids = unique(sol.ie(:)');
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
    end

    % Update state
    T_prev     = T_achieved;
    T_achieved = sol.x(end);
    Z_end      = sol.y(:, end);

    % Record segment timing for ETA estimation
    seg_walltimes(end+1) = seg_wall;            
    seg_dtstar(end+1)    = T_achieved - T_prev; 

    % Renormalise quaternion
    q_end         = Z_end(4:7);
    q_norm_before = norm(q_end);
    Z_end(4:7)    = q_end / q_norm_before;

    % Enforce geometric clamped BC (correct constraint drift)
    Z_end = enforce_clamped_geometry(Z_end, params);

    % ---- Apply Brownian kick at Brownian node (fixed dt_B) ----
    %  Full-system effective friction: dynamically computes the 6x6
    %  effective body mobility (including pili drag) at each Brownian
    %  step via probing the coupled linear system.  The 6D noise
    %  naturally includes translation-rotation coupling.
    brown_kick_applied = false;
    if params.brown_enabled
        brown_tol = dt_B_dimless * 1e-8;
        if T_achieved >= T_brown_next - brown_tol

            if params.brown_fixman_midpoint
                % --- Fixman midpoint correction (2x mobility evaluation) ---
                L_0_chol = compute_effective_body_brownian(Z_end, params);
                eta_6d   = randn(6, 1);
                dU_trial = L_0_chol * eta_6d * sqrt(dt_B);

                % Half-step predictor state
                Z_mid = Z_end;
                Z_mid(1:3) = Z_mid(1:3) + 0.5 * dU_trial(1:3) / params.scales.L_0;
                dTheta_half = 0.5 * dU_trial(4:6);
                ang_h = norm(dTheta_half);
                if ang_h > 1e-15
                    dq_h = [cos(ang_h/2); (sin(ang_h/2)/ang_h) * dTheta_half];
                else
                    dq_h = [1; 0; 0; 0];
                end
                Z_mid(4:7) = quat_multiply(dq_h, Z_mid(4:7));
                Z_mid(4:7) = Z_mid(4:7) / norm(Z_mid(4:7));
                Z_mid = enforce_clamped_geometry(Z_mid, params);

                % Corrector: recompute L at midpoint, SAME eta
                L_mid  = compute_effective_body_brownian(Z_mid, params);
                dU_dim = L_mid * eta_6d * sqrt(dt_B);
            else
                % --- Simple per-step evaluation (no drift correction) ---
                L_eff_step = compute_effective_body_brownian(Z_end, params);
                eta_6d     = randn(6, 1);
                dU_dim     = L_eff_step * eta_6d * sqrt(dt_B);
            end

            % -- Apply translational displacement (lab frame) --
            dX_lab     = dU_dim(1:3);
            Z_end(1:3) = Z_end(1:3) + dX_lab / params.scales.L_0;

            % -- Apply rotational displacement (lab-frame angle -> quaternion) --
            dTheta_lab = dU_dim(4:6);
            ang        = norm(dTheta_lab);
            if ang > 1e-15
                Delta_q = [cos(ang / 2); (sin(ang / 2) / ang) * dTheta_lab];
            else
                Delta_q = [1; 0; 0; 0];
            end
            Z_end(4:7) = quat_multiply(Delta_q, Z_end(4:7));
            Z_end(4:7) = Z_end(4:7) / norm(Z_end(4:7));

            % Re-enforce clamped BC (body quaternion changed)
            Z_end = enforce_clamped_geometry(Z_end, params);

            brown_kick_applied = true;
            brown_kick_count   = brown_kick_count + 1;
            T_brown_next       = (brown_kick_count + 1) * dt_B_dimless;

            if brown_kick_count <= 5 || mod(brown_kick_count, brown_log_interval) == 0
                logf('[BROWN] Kick #%d at t* = %.4e | |dX| = %.2e m | |dTheta| = %.2e rad\n', ...
                    brown_kick_count, T_achieved, norm(dX_lab), ang);
            end
        end
    end
    % ---- End Brownian kick ----

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

            % Override with post-kick state at Brownian node
            if brown_kick_applied && abs(valid_ts(i) - T_achieved) < 1e-14
                Z_i = Z_end;
            end

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
            end
        end
    end

    current_time_ind = find(ts_dimless <= T_achieved + 1e-14, 1, 'last');

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
end

elapsed_total = toc(tic_total);

logf('\n================================================================\n');
logf('  INTEGRATION COMPLETED\n');
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

% --- Sinking velocity [um/s] (Z-component via numerical gradient) ---
V_sink_all      = gradient(results.X_cm_um(3,:), results.t);
results.V_sink  = V_sink_all;

% --- Full translational velocity [um/s]: all 3 components ---
%  Row-by-row gradient avoids any dimension-convention ambiguity of the
%  multi-argument gradient(F, h1, h2) call on a 2-D matrix.
results.V_cm_um = [ gradient(results.X_cm_um(1,:), results.t); ...
                    gradient(results.X_cm_um(2,:), results.t); ...
                    gradient(results.X_cm_um(3,:), results.t) ];  % [3 x NT]

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
body_color    = [0.20 0.65 0.33];    % Forest green
body_alpha    = 0.35;                 % Body transparency
traj_color    = [0.85 0.25 0.25 0.45];  % Deeper red, 45 % opacity
                                         % without occluding pili under OpenGL depth testing

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

fig_summary = figure('Position', [100, 100, 1600, 1000], 'Color', 'w', ...
                     'Visible', 'on');
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

% Time-gradient colormaps: light (early) → dark (late)
body_cmap = [linspace(0.75, 0.10, n_snapshots)', ...   % R
             linspace(0.92, 0.55, n_snapshots)', ...   % G
             linspace(0.78, 0.25, n_snapshots)'];      % B  (light→dark green)

% Per-pili time-gradient colormaps: light (early) → dark (late)
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
        all_x = [all_x; pts_snap(1,:)']; 
        all_y = [all_y; pts_snap(2,:)']; 
        all_z = [all_z; pts_snap(3,:)']; 
        for k = 1:M_pili
            X_p = results.X_pili_um{k}(:,:,ti);
            all_x = [all_x; X_p(:,1)]; 
            all_y = [all_y; X_p(:,2)]; 
            all_z = [all_z; X_p(:,3)]; 
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
            % Custom body: use trisurf with auto-generated triangulation
            rotated_verts = (R_snap * custom_verts_um' + cm_um)';  % [V x 3]
            for ai = 1:numel(snap_axes)
                trisurf(snap_axes(ai), params.body_faces, ...
                    rotated_verts(:,1), rotated_verts(:,2), rotated_verts(:,3), ...
                    'FaceColor', body_cmap(si,:), 'FaceAlpha', snap_alpha, ...
                    'EdgeColor', 'none');
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

sgtitle(sprintf('Coupled Body-Pili Simulation  [%s]', run_id), ...
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
    logf('\n[ANIM] Generating 3D motion animation...\n');
    t_anim_start = tic;

    % Determine frame indices
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

    % Compute global axis limits from ALL frames (body + pili + CM)
    all_x = cm_traj_um(1,:)';
    all_y = cm_traj_um(2,:)';
    all_z = cm_traj_um(3,:)';
    for idx = 1:n_frames
        ti   = frame_indices(idx);
        R_fi = quat_to_rotmat(results.q(:, ti));
        if strcmp(body_mesh_type, 'parametric_surf')
            pts  = R_fi * capsule_flat + results.X_cm_um(:, ti);  % [3 x N_mesh]
        else
            pts  = R_fi * custom_verts_um' + results.X_cm_um(:, ti); % [3 x N_verts]
        end
        all_x = [all_x; pts(1,:)']; 
        all_y = [all_y; pts(2,:)']; 
        all_z = [all_z; pts(3,:)']; 
        for k = 1:M_pili
            Xp = results.X_pili_um{k}(:,:,ti);
            all_x = [all_x; Xp(:,1)]; 
            all_y = [all_y; Xp(:,2)]; 
            all_z = [all_z; Xp(:,3)]; 
        end
    end

    % Build per-axis fixed limits with independent margins
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

    logf('[ANIM] Global axis limits (per-axis): X=[%.2f, %.2f], Y=[%.2f, %.2f], Z=[%.2f, %.2f] um\n', ...
        ax_lim_x(1), ax_lim_x(2), ax_lim_y(1), ax_lim_y(2), ax_lim_z(1), ax_lim_z(2));

    % --- Compute figure dimensions to minimise left/right whitespace ----------
    % With DataAspectRatio=[1 1 1], MATLAB centres the equal-aspect data box
    % inside the axes Position rectangle.  When the projected data content is
    % taller than wide (e.g. a sinking bacterium with large Z range), a fixed
    % landscape figure leaves large blank bands on the sides.
    %
    % Fix: project the 3D bounding box onto the screen at view(45,25), derive
    % the natural content aspect ratio, and size the figure accordingly so the
    % data box fills it without padding.
    %
    % The axes inner area occupies approximately 70 % (width) x 78 % (height)
    % of the figure in a single-axes 3-D plot with title + axis labels.
    % Inverting these fractions converts projected data dimensions to figure px.
    az_v = 45;  el_v = 25;
    proj_w_data = diff(ax_lim_x) * abs(cosd(az_v)) + diff(ax_lim_y) * abs(sind(az_v));
    proj_h_data = (diff(ax_lim_x) * abs(sind(az_v)) + diff(ax_lim_y) * abs(cosd(az_v))) ...
                  * abs(sind(el_v)) + diff(ax_lim_z) * abs(cosd(el_v));
    axes_fill_w = 0.70;   % approx fraction of figure width used by axes
    axes_fill_h = 0.78;   % approx fraction of figure height used by axes
    % Ideal figure size from data projection (single-panel reference)
    anim_fig_h_ref = 1440;
    fig_w_single   = round(anim_fig_h_ref * (proj_w_data / proj_h_data) ...
                           * (axes_fill_h / axes_fill_w));
    % 1×3 panel layout: triple the single-panel width; enforce minimum and H.264 even
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
        logf('[ANIM] Ideal size %d x %d exceeds codec limit [%d x %d] — rescaled to %d x %d.\n', ...
            anim_fig_w, anim_fig_h, anim_max_w, anim_max_h, anim_fig_w_new, anim_fig_h_new);
        anim_fig_w = anim_fig_w_new;
        anim_fig_h = anim_fig_h_new;
    end
    logf('[ANIM] Figure size (logical): %d x %d px  (1x3 layout: 3D + XZ + YZ)\n', ...
        anim_fig_w, anim_fig_h);

    % Create (invisible) animation figure — 3-panel layout
    fig_anim = figure('Units', 'pixels', 'Position', [50, 50, anim_fig_w, anim_fig_h], ...
                      'Color', 'w', 'Visible', 'off');
    set(fig_anim, 'DefaultAxesFontName',  'Arial', ...
                  'DefaultAxesFontSize',  11, ...
                  'DefaultAxesLineWidth', 0.8);

    % Three axes in a single row: 3D (left), XZ (centre), YZ (right)
    ax_3d = axes(fig_anim, 'Position', [0.03, 0.08, 0.30, 0.84]);
    ax_xz = axes(fig_anim, 'Position', [0.36, 0.08, 0.30, 0.84]);
    ax_yz = axes(fig_anim, 'Position', [0.69, 0.08, 0.30, 0.84]);

    % Pre-configure 3D axes so Position is settled before gravity indicator
    xlim(ax_3d, ax_lim_x);  ylim(ax_3d, ax_lim_y);  zlim(ax_3d, ax_lim_z);
    set(ax_3d, 'DataAspectRatio', [1 1 1], 'TickDir', 'out');
    grid(ax_3d, 'on');  box(ax_3d, 'on');
    view(ax_3d, 45, 25);
    xlabel(ax_3d, 'X (μm)', 'FontName', 'Arial', 'FontSize', 12);
    ylabel(ax_3d, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 12);
    zlabel(ax_3d, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 12);
    title(ax_3d, sprintf('Time = %.4f s  (frame %d/%d)', 0.0, 1, n_frames), ...
        'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');

    % Pre-configure XZ projection axes (view along Y axis)
    xlim(ax_xz, ax_lim_x);  ylim(ax_xz, ax_lim_y);  zlim(ax_xz, ax_lim_z);
    set(ax_xz, 'DataAspectRatio', [1 1 1], 'TickDir', 'out', ...
        'Projection', 'orthographic');
    grid(ax_xz, 'on');  box(ax_xz, 'off');
    view(ax_xz, 0, 0);
    xlabel(ax_xz, 'X (μm)', 'FontName', 'Arial', 'FontSize', 11);
    zlabel(ax_xz, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
    title(ax_xz, 'XZ Projection', ...
        'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

    % Pre-configure YZ projection axes (view along X axis)
    xlim(ax_yz, ax_lim_x);  ylim(ax_yz, ax_lim_y);  zlim(ax_yz, ax_lim_z);
    set(ax_yz, 'DataAspectRatio', [1 1 1], 'TickDir', 'out', ...
        'Projection', 'orthographic');
    grid(ax_yz, 'on');  box(ax_yz, 'off');
    view(ax_yz, 90, 0);
    ylabel(ax_yz, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 11);
    zlabel(ax_yz, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
    title(ax_yz, 'YZ Projection', ...
        'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

    drawnow;
    add_gravity_indicator_outside(fig_anim, ax_3d, 13, 2.5);

    % --- Probe actual physical pixel size and enforce codec limits --------
    % On HiDPI / display scaling (e.g. 125%), getframe returns physical pixels
    % larger than the logical Position (e.g. 2250×1366 vs 1800×1200).  The codec
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
        logf('[ANIM] Physical size exceeds codec limit [%d x %d] — shrinking figure to %d x %d logical.\n', ...
            anim_max_w, anim_max_h, new_log_w, new_log_h);
        set(fig_anim, 'Position', [50, 50, new_log_w, new_log_h]);
        anim_fig_w = new_log_w;
        anim_fig_h = new_log_h;
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
        logf(['[ANIM] Captured frame still exceeds MPEG-4 constraints; ', ...
              'encoding will use %d x %d px and resize frames as needed.\n'], ...
              vid_w, vid_h);
    end

    % Now open VideoWriter — frame size is known
    mp4_file = fullfile(output_dir, 'bacteria_3d_motion.mp4');
    vid = VideoWriter(mp4_file, 'MPEG-4');
    vid.FrameRate = animation_fps;
    vid.Quality   = 100;
    open(vid);

    ax_all = [ax_3d, ax_xz, ax_yz];

    for fi = 1:n_frames
        ti = frame_indices(fi);

        % --- Compute geometry once per frame ---
        R_fi  = quat_to_rotmat(results.q(:, ti));
        cm_um = results.X_cm_um(:, ti);

        % --- Draw bacterium on all three axes ---
        for ai = 1:3
            ax = ax_all(ai);
            cla(ax);
            hold(ax, 'on');

            % Body surface rendering (dispatch by mesh type)
            if strcmp(body_mesh_type, 'parametric_surf')
                pts   = R_fi * capsule_flat + cm_um;       % [3 x N_mesh]
                Xsurf_f = reshape(pts(1,:), mesh_size);
                Ysurf_f = reshape(pts(2,:), mesh_size);
                Zsurf_f = reshape(pts(3,:), mesh_size);
                surf(ax, Xsurf_f, Ysurf_f, Zsurf_f, ...
                     'FaceColor', body_color, 'FaceAlpha', body_alpha, ...
                     'EdgeColor', 'none', 'FaceLighting', 'gouraud');
            else
                rotated_verts = (R_fi * custom_verts_um' + cm_um)';  % [V x 3]
                trisurf(ax, params.body_faces, ...
                    rotated_verts(:,1), rotated_verts(:,2), rotated_verts(:,3), ...
                    'FaceColor', body_color, 'FaceAlpha', body_alpha, ...
                    'EdgeColor', 'none', 'FaceLighting', 'gouraud');
            end

            for k = 1:M_pili
                Xp   = results.X_pili_um{k}(:,:,ti);
                pk_c = pili_colors(k, :);
                plot3(ax, Xp(:,1), Xp(:,2), Xp(:,3), ...
                      '-', 'Color', pk_c, 'LineWidth', 1.25);
            end

            plot3(ax, cm_um(1), cm_um(2), cm_um(3), ...
                  'r+', 'MarkerSize', 5, 'LineWidth', 1);

            if fi > 1
                plot3(ax, cm_traj_full_um(1,1:ti), ...
                          cm_traj_full_um(2,1:ti), ...
                          cm_traj_full_um(3,1:ti), ...
                          '-', 'Color', traj_color, 'LineWidth', 1.0);
            end

            hold(ax, 'off');
        end

        % --- Format 3D view ---
        camlight(ax_3d, 'headlight');
        camlight(ax_3d, 'left');
        lighting(ax_3d, 'gouraud');
        material(ax_3d, 'dull');
        xlim(ax_3d, ax_lim_x);  ylim(ax_3d, ax_lim_y);  zlim(ax_3d, ax_lim_z);
        set(ax_3d, 'DataAspectRatio', [1 1 1]);
        grid(ax_3d, 'on');  box(ax_3d, 'on');
        set(ax_3d, 'TickDir', 'out');
        view(ax_3d, 45, 25);
        xlabel(ax_3d, 'X (μm)', 'FontName', 'Arial', 'FontSize', 12);
        ylabel(ax_3d, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 12);
        zlabel(ax_3d, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 12);
        title(ax_3d, sprintf('Time = %.4f s  (frame %d/%d)', ...
            results.t(ti), fi, n_frames), ...
            'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');

        % --- Format XZ projection ---
        camlight(ax_xz, 'headlight');
        lighting(ax_xz, 'gouraud');
        material(ax_xz, 'dull');
        xlim(ax_xz, ax_lim_x);  ylim(ax_xz, ax_lim_y);  zlim(ax_xz, ax_lim_z);
        set(ax_xz, 'DataAspectRatio', [1 1 1], 'Projection', 'orthographic');
        grid(ax_xz, 'on');  box(ax_xz, 'off');
        set(ax_xz, 'TickDir', 'out');
        view(ax_xz, 0, 0);
        hold(ax_xz, 'on');
        plot3(ax_xz, [ax_lim_x(2) ax_lim_x(2)], ax_lim_y([1 1]), ax_lim_z, ...
            'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        plot3(ax_xz, ax_lim_x, ax_lim_y([1 1]), [ax_lim_z(2) ax_lim_z(2)], ...
            'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        hold(ax_xz, 'off');
        xlabel(ax_xz, 'X (μm)', 'FontName', 'Arial', 'FontSize', 11);
        zlabel(ax_xz, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
        title(ax_xz, 'XZ Projection', ...
            'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

        % --- Format YZ projection ---
        camlight(ax_yz, 'headlight');
        lighting(ax_yz, 'gouraud');
        material(ax_yz, 'dull');
        xlim(ax_yz, ax_lim_x);  ylim(ax_yz, ax_lim_y);  zlim(ax_yz, ax_lim_z);
        set(ax_yz, 'DataAspectRatio', [1 1 1], 'Projection', 'orthographic');
        grid(ax_yz, 'on');  box(ax_yz, 'off');
        set(ax_yz, 'TickDir', 'out');
        view(ax_yz, 90, 0);
        hold(ax_yz, 'on');
        plot3(ax_yz, ax_lim_x([1 1]), [ax_lim_y(2) ax_lim_y(2)], ax_lim_z, ...
            'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        plot3(ax_yz, ax_lim_x([1 1]), ax_lim_y, [ax_lim_z(2) ax_lim_z(2)], ...
            'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        hold(ax_yz, 'off');
        ylabel(ax_yz, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 11);
        zlabel(ax_yz, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
        title(ax_yz, 'YZ Projection', ...
            'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

        drawnow;

        % Capture frame (no rect — let MATLAB return natural physical pixels)
        frame = getframe(fig_anim);
        [fh, fw, ~] = size(frame.cdata);
        if fh ~= vid_h || fw ~= vid_w
            frame.cdata = normalize_anim_frame_cdata(frame.cdata, vid_h, vid_w);
        end
        writeVideo(vid, frame);

        % In-place progress bar (console only)
        render_progress_bar('[ANIM]   ', fi, n_frames, t_anim_start);

        % Log animation progress to file every 10% (file-only, no console echo)
        if mod(fi, max(1, floor(n_frames / 10))) == 0 || fi == n_frames
            logf_file('[ANIM] Frame %d / %d rendered\n', fi, n_frames);
        end
    end

    close(vid);
    close(fig_anim);

    t_anim_elapsed = toc(t_anim_start);
    vid_info = dir(mp4_file);
    logf('[ANIM] Animation saved: %s (%.2f MB, %d frames, %.1f s render)\n', ...
        mp4_file, vid_info.bytes/1e6, n_frames, t_anim_elapsed);
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
    compute_plot_flowfield(last_Z, params, output_dir, run_id, logf, flow_opts);
    t_flow_elapsed = toc(t_flow_start);
    logf('[FLOW] Flow field visualisation completed in %s\n\n', format_duration(t_flow_elapsed));
else
    t_flow_elapsed = 0;
    logf('\n[FLOW] Skipped (cfg.compute_flowfield = false).\n\n');
end

%% =====================================================================
%  FINAL SUMMARY
%  =====================================================================

total_sim_walltime = toc(tic_sim_start);

logf('\n================================================================\n');
logf('  SIMULATION COMPLETE\n');
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
    logf('    Video  : bacteria_3d_motion.mp4\n');
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
fprintf('Log saved: %s\n', log_file);

end  % run_simulation


%% ==================== Local Functions ====================

function v = get_cfg_field(cfg, fname, default)
%GET_CFG_FIELD  Return cfg.(fname) if present and non-empty, else default
    if isstruct(cfg) && isfield(cfg, fname) && ~isempty(cfg.(fname))
        v = cfg.(fname);
    else
        v = default;
    end
end

function dual_log(fid, varargin)
%DUAL_LOG  Write formatted message to both console and log file (with timestamp)
    msg = sprintf(varargin{:});
    ts  = char(datetime('now', 'Format', 'HH:mm:ss'));
    fprintf('%s', msg);                    % Console (no timestamp — cleaner)
    fprintf(fid, '[%s] %s', ts, msg);     % File (with timestamp)
end

function log_file_only(fid, varargin)
%LOG_FILE_ONLY  Write formatted message to log file only (with timestamp)
    msg = sprintf(varargin{:});
    ts  = char(datetime('now', 'Format', 'HH:mm:ss'));
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

function cdata = normalize_anim_frame_cdata(cdata, nh, nw)
%NORMALIZE_ANIM_FRAME_CDATA  Resize RGB(A) frame to exact [nh x nw] for VideoWriter.
%  On Windows Hi-DPI, getframe often returns physical pixels larger than figure Position;
%  without this, writeVideo fails ("frames must be WxH").
    if isempty(cdata)
        error('run_simulation:EmptyAnimFrame', 'getframe returned empty cdata');
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
        error('run_simulation:InvalidAnimFrameSize', ...
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

function add_gravity_indicator_outside(fig_handle, ax_handle, font_size, line_width)
%ADD_GRAVITY_INDICATOR_OUTSIDE  Draw gravity arrow flush outside axis right edge.
%  Uses figure-normalised annotation coordinates so marker location is
%  stable across views and axis limits.
%
%  Root-cause fix (v2.4):
%    With DataAspectRatio=[1 1 1] and a strongly elongated data range
%    (e.g. Z >> X,Y for a sinking bacterium), MATLAB renders a very narrow
%    axis box CENTRED inside the Position rectangle.  Using the Position
%    right edge (the old approach) places the arrow hundreds of pixels to
%    the right of the actual box edge.
%
%    This function instead computes the TRUE rendered right edge by:
%      1. Obtaining the figure size in pixels (aspect ratio).
%      2. Projecting the data bounding box (xlim x ylim x zlim) onto the
%         screen-horizontal direction via the orthographic approximation:
%           right   = [cos(az), sin(az), 0]
%           data_w  = dx|cos az| + dy|sin az|
%           data_h  = (dx|sin az|+dy|cos az|)|sin el| + dz|cos el|
%      3. Applying the equal-aspect scale:
%           px_per_unit = min(ax_w_px/data_w, ax_h_px/data_h)
%      4. Computing the centred right edge:
%           actual_right = ax_pos(1) + ax_pos(3)/2
%                        + data_w*px_per_unit / (2*fig_w_px)
%
%  IMPORTANT:
%    1. Must be called AFTER any sgtitle() on the parent figure.
%    2. Must be called AFTER the axes is fully configured (xlim/ylim/zlim,
%       DataAspectRatio, view, labels, title) so ax.Position is settled.

    drawnow;   % Flush pending layout so ax.Position is up-to-date

    ax_pos = get(ax_handle, 'Position');   % [left bottom width height] normalised

    % ---- Figure size in pixels (needed for aspect-ratio computation) ----
    old_units = get(fig_handle, 'Units');
    set(fig_handle, 'Units', 'pixels');
    fig_pos_px = get(fig_handle, 'Position');   % [x y w h] in screen pixels
    set(fig_handle, 'Units', old_units);
    fig_w_px = fig_pos_px(3);
    fig_h_px = fig_pos_px(4);

    % ---- Axis Position in pixels ----------------------------------------
    ax_w_px = ax_pos(3) * fig_w_px;
    ax_h_px = ax_pos(4) * fig_h_px;

    % ---- Data bounding-box ranges (equal physical units on all axes) ----
    dx = abs(diff(xlim(ax_handle)));
    dy = abs(diff(ylim(ax_handle)));
    dz = abs(diff(zlim(ax_handle)));

    % ---- Orthographic projection of the data box onto screen axes -------
    %  Screen right-vector:  r = [cos(az), sin(az), 0]
    %  Screen up-vector:     u = [-sin(az)*sin(el), cos(az)*sin(el), cos(el)]
    %
    %  Projected extents of the axis-aligned bounding box:
    %    data_w = dx|cos az| + dy|sin az|
    %    data_h = (dx|sin az|+dy|cos az|)*|sin el| + dz*|cos el|
    [az_deg, el_deg] = view(ax_handle);
    ca = abs(cos(az_deg * pi/180));   sa = abs(sin(az_deg * pi/180));
    ce = abs(cos(el_deg * pi/180));   se = abs(sin(el_deg * pi/180));

    data_w = dx * ca + dy * sa;
    data_h = (dx * sa + dy * ca) * se + dz * ce;
    data_w = max(data_w, 1e-12);   % guard against degenerate limits
    data_h = max(data_h, 1e-12);

    % ---- Equal-aspect scale factor: whichever axis fills first ----------
    px_per_unit = min(ax_w_px / data_w,  ax_h_px / data_h);

    % ---- Actual rendered half-width in normalised figure coords ---------
    % MATLAB centres the equal-aspect box inside the Position rectangle.
    actual_half_w = (data_w * px_per_unit) / (2 * fig_w_px);
    actual_right  = ax_pos(1) + ax_pos(3) / 2 + actual_half_w;

    % Place arrow 1 % of figure width to the right of the rendered box edge
    x_arrow = actual_right + 0.010;
    y_mid   = ax_pos(2) + 0.50 * ax_pos(4);
    y_top   = y_mid + 0.055;
    y_bot   = y_mid - 0.055;

    % Keep annotation within figure bounds
    x_arrow = min(max(x_arrow, 0.02), 0.98);
    y_top   = min(max(y_top,   0.05), 0.98);
    y_bot   = min(max(y_bot,   0.02), 0.95);
    if y_top <= y_bot
        y_top = min(0.98, y_bot + 0.06);
    end

    annotation(fig_handle, 'arrow', [x_arrow, x_arrow], [y_top, y_bot], ...
        'Color', [0.25 0.25 0.25], 'LineWidth', line_width, 'HeadLength', 12, ...
        'HeadWidth', 10);
    % "G" label centred directly below the arrow shaft
    lbl_w = 0.030;
    lbl_h = 0.030;
    annotation(fig_handle, 'textbox', ...
        [x_arrow - lbl_w/2, y_bot - lbl_h - 0.004, lbl_w, lbl_h], ...
        'String', '{\itG}', 'Interpreter', 'tex', 'FontName', 'Arial', ...
        'FontSize', font_size, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');
end

function [value, isterminal, direction] = singularity_event(~, Z, params, tstart, tlim)
%SINGULARITY_EVENT  ODE event: Euler angle singularity or timeout
    delta  = params.delta;
    N_segs = params.N_segments;
    M      = params.M_pili;

    % Check timeout
    elapsed     = toc(tstart);
    timeout_val = tlim - elapsed;

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
