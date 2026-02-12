%% main_coupled_simulation.m
%  Coupled Body-Pili Elastohydrodynamic Simulation
%
%  Description:
%    Simulates a capsule-shaped bacterial body with elastic Kirchhoff rod
%    pili sinking under gravity in Stokes flow.  The body is rigid; the
%    pili are inextensible elastic rods parameterised by Euler angles.
%
%    Hydrodynamic interactions (body-body, body-pili, pili-pili) are fully
%    coupled through the Regularised Stokeslet Method (RSM).
%
%    Time integration uses MATLAB's ode15s (implicit stiff solver) with
%    adaptive Euler angle basis rotation to avoid singularities.
%
%  Outputs (saved to output/run_YYYYMMDD_HHMMSS/):
%    simulation_results.mat  results + params structures
%    summary_plot.png        static overview figure
%    bacteria_3d_motion.mp4  3D animation
%    simulation.log          full run log
%
%  Output unit convention:
%    Time [s] | Length [um] | Force [pN] | Angle [rad]
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] Cortez et al. (2005), Physics of Fluids, 17, 031504.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Logging, data export, 3D animation;
%                         removed fault-tolerance retry; code formatting;
%                         added progress timer with ETA estimation;
%                         added in-segment ODE progress bar (ode_progress)

close all force;
clear; clc;

tic_sim_start = tic;        % Global wall clock for the entire run

%% =====================================================================
%  CASE IDENTIFICATION (must precede output directory setup)
%  =====================================================================
case_name = '';                    % User-defined case name (string)
                                   %   Non-empty : output saved to output/<case_name>/
                                   %   Empty ''  : output saved to output/run_YYYYMMDD_HHMMSS/

%% =====================================================================
%  OUTPUT DIRECTORY AND LOGGING SETUP
%  =====================================================================

run_timestamp = char(datetime('now', 'Format', 'uuuuMMdd''_''HHmmss'));
if isempty(strtrim(case_name))
    run_id     = ['run_' run_timestamp];
    output_dir = fullfile(pwd, 'output', run_id);
else
    run_id     = strtrim(case_name);
    output_dir = fullfile(pwd, 'output', run_id);
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% --- Initialise dual-channel logger (console + file) ---
log_file = fullfile(output_dir, 'simulation.log');
log_fid  = fopen(log_file, 'w');
if log_fid == -1
    error('Cannot create log file: %s', log_file);
end
logf = @(varargin) dual_log(log_fid, varargin{:});

logf('\n');
logf('================================================================\n');
logf('  Coupled Body-Pili Elastohydrodynamic Simulation\n');
logf('  Version 1.2\n');
logf('================================================================\n');
logf('  Run ID    : %s\n', run_id);
logf('  Timestamp : %s\n', run_timestamp);
logf('  Output dir: %s\n', output_dir);
logf('  Platform  : MATLAB %s on %s\n', version, computer);
logf('================================================================\n\n');

%% =====================================================================
%  PARAMETER SETTINGS
%  =====================================================================

logf('[PARAMS] Loading simulation parameters...\n');

% -------- Bacterial Body Geometry --------
minor_axis = 0.4e-6;              % [m]   Semi-minor axis (body radius)
major_axis = 1.25e-6;             % [m]   Semi-major axis (half body length)
N_body     = 1000;                % [-]   Number of body surface points
alpha_0    = pi/2;                % [rad] Initial tilt angle in XZ plane
                                  %       Angle between body long axis and +X
                                  %       pi/2 = vertical (+Z), 0 = horizontal (+X)
                                  %       -pi/2 = inverted (-Z)
                                  %       Range: [-pi/2, pi/2]

% -------- Pili Configuration --------
%  Each column defines one pilus:
%    Row 1 : Type (1=linear, 2=parabolic, 3=circular, 4=sinusoidal, 5=spiral)
%    Row 2 : Length [um]
%    Row 3 : Deflection angle from Z-axis [rad]
%    Row 4 : Deflection angle from X-axis [rad]
%    Row 5-7: Shape parameters (NaN if unused)
%    Row 8 : Attachment site (1-14)
Pili_Matrix = [1;                  % Type: linear (straight)
               3;                  % Length: 3 um
               0;                  % Deflection from Z: 0
               0;                  % Deflection from X: 0
               NaN; NaN; NaN;     % Shape parameters (unused)
               3];                 % Attachment: Pole 1

seg_per_um       = 150;            % [seg/um] Pili discretisation density
pili_lengths_um  = Pili_Matrix(2, :);
N_segments       = round(pili_lengths_um * seg_per_um);
N_segments       = max(N_segments, 10);

% -------- Material Properties --------
E       = 3.67e8;                 % [Pa]     Young's modulus of pili
r_pili  = 3e-9;                   % [m]      Pili radius
sigma_p = 0.2;                    % [-]      Poisson's ratio

% -------- Fluid Properties --------
mu = 1e-3;                        % [Pa*s]   Dynamic viscosity (water)

% -------- Density --------
rho_bac   = 1180.0;               % [kg/m^3] Bacteria density
rho_fluid = 998.232;              % [kg/m^3] Fluid density

% -------- Simulation Time --------
T_obs = 180.0;                     % [s]      Total observation time

% -------- Numerical Settings --------
AbsTol = 1e-4;                    % ODE absolute tolerance
RelTol = 1e-4;                    % ODE relative tolerance
delta  = pi / 20;                 % [rad]    Euler angle singularity threshold
tlim   = 120;                     % [s]      ODE timeout limit
dt     = 1e-1;                    % [s]      Output time step
NT     = round(T_obs / dt) + 1;  % [-]      Number of output time points (derived)

% -------- Animation Settings --------
generate_animation = true;
animation_fps      = 30;           % Frames per second for MP4
animation_skip     = 1;            % Frame skip interval (every Nth output point)

% --- Log all parameters ---
logf('[PARAMS] Body: a = %.2f um, b = %.2f um, N_body = %d, alpha_0 = %.4f rad (%.1f deg)\n', ...
    major_axis*1e6, minor_axis*1e6, N_body, alpha_0, rad2deg(alpha_0));
logf('[PARAMS] Pili: M = %d, seg_per_um = %d\n', ...
    size(Pili_Matrix, 2), seg_per_um);
for kk = 1:size(Pili_Matrix, 2)
    logf('[PARAMS]   Pili %d: type=%d, L=%.1f um, site=%d, N_seg=%d\n', ...
        kk, Pili_Matrix(1,kk), Pili_Matrix(2,kk), ...
        Pili_Matrix(8,kk), N_segments(kk));
end
logf('[PARAMS] Material: E = %.2e Pa, r_pili = %.1f nm, sigma = %.1f\n', ...
    E, r_pili*1e9, sigma_p);
logf('[PARAMS] Fluid: mu = %.2e Pa*s, rho_bac = %.1f kg/m^3, rho_fluid = %.3f kg/m^3\n', ...
    mu, rho_bac, rho_fluid);
logf('[PARAMS] Simulation: T_obs = %.1f s, dt = %.2e s, NT = %d, AbsTol = %.0e, RelTol = %.0e\n', ...
    T_obs, dt, NT, AbsTol, RelTol);
logf('[PARAMS] Numerics: delta = pi/%.0f rad, timeout = %d s\n', pi/delta, tlim);
logf('[PARAMS] Animation: fps = %d, skip = %d\n', animation_fps, animation_skip);
if isempty(strtrim(case_name))
    logf('[PARAMS] Case name: (auto — timestamp)\n');
else
    logf('[PARAMS] Case name: %s\n', strtrim(case_name));
end

%% =====================================================================
%  INITIALISATION
%  =====================================================================

logf('\n[INIT] Initialising simulation parameters...\n');
t_init_start = tic;

params = init_coupled_params( ...
    'N_body',      N_body, ...
    'minor_axis',  minor_axis, ...
    'major_axis',  major_axis, ...
    'alpha_0',     alpha_0, ...
    'mu',          mu, ...
    'E',           E, ...
    'r_pili',      r_pili, ...
    'sigma',       sigma_p, ...
    'rho_bac',     rho_bac, ...
    'rho_fluid',   rho_fluid, ...
    'N_segments',  N_segments, ...
    'Pili_Matrix', Pili_Matrix, ...
    'seg_per_um',  seg_per_um, ...
    'delta',       delta);

logf('[INIT] Building initial state vector...\n');
Z0 = init_coupled_state(params);

t_init_elapsed = toc(t_init_start);

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
for kk = 1:params.M_pili
    logf('[INIT]   Pili %d: ell = %.4f L_0, N_seg = %d, ds = %.4f L_0\n', ...
        kk, params.ell(kk), params.N_segments(kk), ...
        params.ell(kk) / params.N_segments(kk));
end
logf('[INIT]   alpha_0 = %.4f rad (%.1f deg) — initial body tilt\n', ...
    params.alpha_0, rad2deg(params.alpha_0));
logf('[INIT] State vector : %d DOFs\n', length(Z0));
logf('[INIT] Total points : %d (body + pili)\n', ...
    params.N_body + sum(params.N_segments + 1));
logf('[INIT] Initialisation completed in %.2f s\n\n', t_init_elapsed);

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
for k = 1:M_pili
    results.X_pili{k} = zeros(N_segs(k)+1, 3, NT);
    geom0 = reconstruct_geometry(Z0, params);
    results.X_pili{k}(:, :, 1) = geom0.X_pili{k};
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
        segment_count, T_achieved, T_obs_dimless);

    sol      = ode15s(dZ_func, [T_achieved, T_obs_dimless], Z, ode_ops);
    seg_wall = toc(tstart_ode);
    logf('[ODE] Segment %d completed (%.1f s, %d steps)\n', ...
        segment_count, seg_wall, length(sol.x));

    % Check if event was triggered
    if ~isempty(sol.xe)
        logf('[ODE]   Event triggered at t* = %.4e (singularity or timeout)\n', ...
            sol.xe(end));
    end

    % Update state
    T_prev     = T_achieved;
    T_achieved = sol.x(end);
    Z_end      = sol.y(:, end);

    % Record segment timing for ETA estimation
    seg_walltimes(end+1) = seg_wall;            %#ok<SAGROW>
    seg_dtstar(end+1)    = T_achieved - T_prev; %#ok<SAGROW>

    % Renormalise quaternion
    q_end         = Z_end(4:7);
    q_norm_before = norm(q_end);
    Z_end(4:7)    = q_end / q_norm_before;

    % Enforce geometric clamped BC (correct constraint drift)
    Z_end  = enforce_clamped_geometry(Z_end, params);
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
            end
        end
    end

    current_time_ind = find(ts_dimless <= T_achieved + 1e-14, 1, 'last');

    % Log progress with ETA
    progress_pct = min(100, T_achieved / T_obs_dimless * 100);
    elapsed_ode  = toc(tic_total);

    if progress_pct > 0 && progress_pct < 100
        % ETA from linear extrapolation: elapsed * (remaining / done)
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

% --- Sinking velocity [um/s] ---
V_sink_all      = gradient(results.X_cm_um(3,:), results.t);
results.V_sink  = V_sink_all;

logf('[POST] Physical unit conversion complete.\n');
logf('[POST] Sinking force     : %.4f pN\n',   results.F_sink_pN);
logf('[POST] Final Z position  : %.4e um\n',    results.X_cm_um(3, end));
logf('[POST] Mean V_sink       : %.4e um/s\n',  mean(V_sink_all));
logf('[POST] Steady V_sink     : %.4e um/s (last 10%%)\n', ...
    mean(V_sink_all(round(0.9*NT):end)));
logf('\n');

%% =====================================================================
%  SAVE DATA (.mat)
%  =====================================================================

logf('[SAVE] Saving simulation data...\n');

mat_file = fullfile(output_dir, 'simulation_results.mat');
save(mat_file, 'results', 'params', '-v7.3');
mat_info = dir(mat_file);
logf('[SAVE] Data saved: %s (%.2f MB)\n', mat_file, mat_info.bytes/1e6);

%% =====================================================================
%  PRE-COMPUTE CAPSULE SURFACE MESH (for visualisation)
%  =====================================================================
%  Generate a smooth parametric capsule (spherocylinder) mesh in the body-
%  local frame.  Units: micrometres [um].
%  The mesh comprises: left hemisphere + cylinder + right hemisphere,
%  concatenated along the meridional direction so that MATLAB's surf()
%  renders a seamless surface.
%  =====================================================================

logf('[PLOT] Building capsule surface mesh for visualisation...\n');

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

logf('[PLOT] Capsule mesh: %d x %d = %d vertices\n', ...
    mesh_size(1), mesh_size(2), numel(capsule_X0));

% --- Colour palette (journal style) ---
body_color    = [0.20 0.65 0.33];    % Forest green
body_alpha    = 0.35;                 % Body transparency
pili_color    = [0.20 0.40 0.80];    % Medium blue (uniform for all pili)
traj_color    = [1.00 0.55 0.55];    % Light red for CM trajectory

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
ylabel('Position ({\mu}m)',        'FontName', 'Arial', 'FontSize', 12);
title('Body Centre of Mass',      'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
lg = legend('X', 'Y', 'Z', 'Location', 'best');
set(lg, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off');
grid on; box on;

% --- Subplot 2: Sinking velocity ---
subplot(2, 3, 2);
plot(t_s, results.V_sink, 'k-', 'LineWidth', 1.8);
xlabel('Time (s)',                        'FontName', 'Arial', 'FontSize', 12);
ylabel('{\itV}_{sink} ({\mu}m/s)',        'FontName', 'Arial', 'FontSize', 12);
title('Sinking Velocity',                 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
grid on; box on;

% --- Subplot 3: Quaternion norm check ---
subplot(2, 3, 3);
q_norms = sqrt(sum(results.q.^2, 1));
plot(t_s, q_norms - 1, 'k-', 'LineWidth', 1.8);
xlabel('Time (s)',                         'FontName', 'Arial', 'FontSize', 12);
ylabel('|{\bfq}| \minus 1',               'FontName', 'Arial', 'FontSize', 12);
title('Quaternion Norm Deviation',         'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
grid on; box on;

% --- Subplot 4-6: 3D trajectory snapshots ---
ax_snap = subplot(2, 3, [4, 5, 6]);
hold on; grid on; box on;
view(45, 30);

n_snapshots  = min(10, NT);
snap_indices = round(linspace(1, NT, n_snapshots));

% Time-gradient colormaps: light (early) → dark (late)
body_cmap = [linspace(0.75, 0.10, n_snapshots)', ...   % R
             linspace(0.92, 0.55, n_snapshots)', ...   % G
             linspace(0.78, 0.25, n_snapshots)'];      % B  (light→dark green)
pili_cmap = [linspace(0.65, 0.10, n_snapshots)', ...   % R
             linspace(0.80, 0.30, n_snapshots)', ...   % G
             linspace(0.95, 0.70, n_snapshots)'];      % B  (light→dark blue)

% Draw CM trajectory in light red
plot3(ax_snap, X_um(1,:), X_um(2,:), X_um(3,:), ...
      '-', 'Color', traj_color, 'LineWidth', 1.5);

for si = 1:n_snapshots
    ti = snap_indices(si);
    if ti <= size(results.X_cm, 2)
        R_snap = quat_to_rotmat(results.q(:, ti));
        cm_um  = results.X_cm_um(:, ti);

        % Smooth capsule surface (gradient green)
        pts_snap = R_snap * capsule_flat + cm_um;
        Xs = reshape(pts_snap(1,:), mesh_size);
        Ys = reshape(pts_snap(2,:), mesh_size);
        Zs = reshape(pts_snap(3,:), mesh_size);
        snap_alpha = 0.15 + 0.30 * (si / n_snapshots);
        surf(ax_snap, Xs, Ys, Zs, ...
             'FaceColor', body_cmap(si,:), 'FaceAlpha', snap_alpha, ...
             'EdgeColor', 'none');

        % Pili polylines (gradient blue) with root markers
        for k = 1:M_pili
            X_p = results.X_pili_um{k}(:,:,ti);
            plot3(X_p(:,1), X_p(:,2), X_p(:,3), ...
                  '-', 'Color', pili_cmap(si,:), 'LineWidth', 1.8);
            plot3(X_p(1,1), X_p(1,2), X_p(1,3), ...
                  'o', 'Color', pili_cmap(si,:), 'MarkerSize', 4, ...
                  'MarkerFaceColor', pili_cmap(si,:));
        end
    end
end

% Lighting for Gouraud shading
camlight(ax_snap, 'headlight');
lighting(ax_snap, 'gouraud');
material(ax_snap, 'dull');

xlabel(ax_snap, 'X ({\mu}m)', 'FontName', 'Arial', 'FontSize', 12);
ylabel(ax_snap, 'Y ({\mu}m)', 'FontName', 'Arial', 'FontSize', 12);
zlabel(ax_snap, 'Z ({\mu}m)', 'FontName', 'Arial', 'FontSize', 12);
title(ax_snap, '3D Trajectory Snapshots', ...
      'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
axis equal;

% --- Gravity direction indicator ---
zL = ax_snap.ZLim;  xL = ax_snap.XLim;  yL = ax_snap.YLim;
g_arrow_len = 0.15 * (zL(2) - zL(1));
g_x0 = xL(1) + 0.08 * (xL(2) - xL(1));
g_y0 = yL(2) - 0.08 * (yL(2) - yL(1));
g_z0 = zL(2) - 0.06 * (zL(2) - zL(1));
quiver3(ax_snap, g_x0, g_y0, g_z0, 0, 0, -g_arrow_len, 0, ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
text(ax_snap, g_x0, g_y0, g_z0 - g_arrow_len - 0.04*(zL(2)-zL(1)), ...
     '{\itg}', 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.25 0.25 0.25]);

sgtitle(sprintf('Coupled Body-Pili Simulation  [%s]', run_id), ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
png_file = fullfile(output_dir, 'summary_plot.png');
exportgraphics(fig_summary, png_file, 'Resolution', 300);
logf('[PLOT] Summary figure saved: %s\n', png_file);

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

    % Compute global axis limits from ALL frames (capsule + pili + CM)
    all_x = cm_traj_um(1,:)';
    all_y = cm_traj_um(2,:)';
    all_z = cm_traj_um(3,:)';
    for idx = 1:n_frames
        ti   = frame_indices(idx);
        R_fi = quat_to_rotmat(results.q(:, ti));
        pts  = R_fi * capsule_flat + results.X_cm_um(:, ti);  % [3 x N_mesh]
        all_x = [all_x; pts(1,:)']; %#ok<AGROW>
        all_y = [all_y; pts(2,:)']; %#ok<AGROW>
        all_z = [all_z; pts(3,:)']; %#ok<AGROW>
        for k = 1:M_pili
            Xp = results.X_pili_um{k}(:,:,ti);
            all_x = [all_x; Xp(:,1)]; %#ok<AGROW>
            all_y = [all_y; Xp(:,2)]; %#ok<AGROW>
            all_z = [all_z; Xp(:,3)]; %#ok<AGROW>
        end
    end

    % Build equal-aspect-ratio fixed axis limits (cubic bounding box)
    x_span = max(all_x) - min(all_x);
    y_span = max(all_y) - min(all_y);
    z_span = max(all_z) - min(all_z);
    max_span = max([x_span, y_span, z_span]);
    margin   = 0.08 * max_span;
    if margin < 1e-6;  margin = 1;  end
    half_box = max_span / 2 + margin;
    x_mid = (max(all_x) + min(all_x)) / 2;
    y_mid = (max(all_y) + min(all_y)) / 2;
    z_mid = (max(all_z) + min(all_z)) / 2;
    ax_lim_x = [x_mid - half_box, x_mid + half_box];
    ax_lim_y = [y_mid - half_box, y_mid + half_box];
    ax_lim_z = [z_mid - half_box, z_mid + half_box];

    logf('[ANIM] Global axis limits (equal-aspect): X=[%.2f, %.2f], Y=[%.2f, %.2f], Z=[%.2f, %.2f] um\n', ...
        ax_lim_x(1), ax_lim_x(2), ax_lim_y(1), ax_lim_y(2), ax_lim_z(1), ax_lim_z(2));

    % Setup video writer
    mp4_file = fullfile(output_dir, 'bacteria_3d_motion.mp4');
    vid = VideoWriter(mp4_file, 'MPEG-4');
    vid.FrameRate = animation_fps;
    vid.Quality   = 95;
    open(vid);

    % Create (invisible) animation figure with OpenGL renderer for transparency
    fig_anim = figure('Position', [100, 100, 1080, 810], 'Color', 'w', ...
                      'Visible', 'off');
    set(fig_anim, 'DefaultAxesFontName',  'Arial', ...
                  'DefaultAxesFontSize',  11, ...
                  'DefaultAxesLineWidth', 0.8);
    ax_anim = axes(fig_anim);

    for fi = 1:n_frames
        ti = frame_indices(fi);
        cla(ax_anim);
        hold(ax_anim, 'on');

        % --- Draw CM trajectory: light red line from start to current ---
        if fi > 1
            plot3(ax_anim, cm_traj_um(1,1:fi), cm_traj_um(2,1:fi), ...
                  cm_traj_um(3,1:fi), '-', 'Color', traj_color, ...
                  'LineWidth', 1.8);
        end

        % --- Draw body: smooth capsule surface with transparency ---
        R_fi  = quat_to_rotmat(results.q(:, ti));
        cm_um = results.X_cm_um(:, ti);
        pts   = R_fi * capsule_flat + cm_um;       % [3 x N_mesh]
        Xsurf = reshape(pts(1,:), mesh_size);
        Ysurf = reshape(pts(2,:), mesh_size);
        Zsurf = reshape(pts(3,:), mesh_size);

        surf(ax_anim, Xsurf, Ysurf, Zsurf, ...
             'FaceColor', body_color, 'FaceAlpha', body_alpha, ...
             'EdgeColor', 'none', 'FaceLighting', 'gouraud');

        % --- Draw pili: blue polylines with root markers ---
        for k = 1:M_pili
            Xp = results.X_pili_um{k}(:,:,ti);
            plot3(ax_anim, Xp(:,1), Xp(:,2), Xp(:,3), ...
                  '-', 'Color', pili_color, 'LineWidth', 2.5);
            % Highlight root node (attachment point on body)
            plot3(ax_anim, Xp(1,1), Xp(1,2), Xp(1,3), ...
                  'o', 'Color', pili_color, 'MarkerSize', 5, ...
                  'MarkerFaceColor', pili_color);
        end

        % --- Mark body centre of mass ---
        plot3(ax_anim, cm_um(1), cm_um(2), cm_um(3), ...
              'r+', 'MarkerSize', 10, 'LineWidth', 2);

        % --- Gravity direction indicator (fixed corner position) ---
        g_arrow_len = 0.13 * (ax_lim_z(2) - ax_lim_z(1));
        g_ax = ax_lim_x(1) + 0.08 * (ax_lim_x(2) - ax_lim_x(1));
        g_ay = ax_lim_y(2) - 0.08 * (ax_lim_y(2) - ax_lim_y(1));
        g_az = ax_lim_z(2) - 0.06 * (ax_lim_z(2) - ax_lim_z(1));
        quiver3(ax_anim, g_ax, g_ay, g_az, 0, 0, -g_arrow_len, 0, ...
                'Color', [0.25 0.25 0.25], 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
        text(ax_anim, g_ax, g_ay, g_az - g_arrow_len - 0.03*(ax_lim_z(2)-ax_lim_z(1)), ...
             '{\itg}', 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'Color', [0.25 0.25 0.25]);

        % --- Lighting ---
        camlight(ax_anim, 'headlight');
        camlight(ax_anim, 'left');
        lighting(ax_anim, 'gouraud');
        material(ax_anim, 'dull');

        % --- Axis formatting (fixed equal-aspect limits) ---
        xlim(ax_anim, ax_lim_x);
        ylim(ax_anim, ax_lim_y);
        zlim(ax_anim, ax_lim_z);
        set(ax_anim, 'DataAspectRatio', [1 1 1]);
        grid(ax_anim, 'on');
        box(ax_anim,  'on');
        set(ax_anim, 'TickDir', 'out');
        view(ax_anim, 45, 25);
        xlabel(ax_anim, 'X ({\mu}m)', 'FontName', 'Arial', 'FontSize', 12);
        ylabel(ax_anim, 'Y ({\mu}m)', 'FontName', 'Arial', 'FontSize', 12);
        zlabel(ax_anim, 'Z ({\mu}m)', 'FontName', 'Arial', 'FontSize', 12);
        title(ax_anim, sprintf('Bacteria Sinking  |  {\\itt} = %.4f s  (frame %d/%d)', ...
            results.t(ti), fi, n_frames), ...
            'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');

        hold(ax_anim, 'off');
        drawnow;

        % Capture frame
        frame = getframe(fig_anim);
        writeVideo(vid, frame);

        % Progress log with ETA (every 10 %)
        if mod(fi, max(1, floor(n_frames/10))) == 0
            elapsed_anim = toc(t_anim_start);
            eta_anim_s   = elapsed_anim * (n_frames - fi) / fi;
            logf('[ANIM]   Frame %d/%d (%.0f%%) | Elapsed: %s | ETA: %s\n', ...
                fi, n_frames, fi/n_frames*100, ...
                format_duration(elapsed_anim), format_duration(eta_anim_s));
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
    mean(V_sink_all(round(0.9*NT):end)));
logf('  Basis rotations  : %d\n',    basis_changes);
logf('  Max |q| drift    : %.2e\n',  max(abs(sqrt(sum(results.q.^2,1))-1)));
logf('  ------------------------------------------------\n');
logf('  Timing breakdown:\n');
logf('    Initialisation  : %s\n',   format_duration(t_init_elapsed));
logf('    ODE integration : %s\n',   format_duration(elapsed_total));
if generate_animation
    logf('    Animation render: %s\n',   format_duration(t_anim_elapsed));
end
logf('    Total wall time : %s\n',   format_duration(total_sim_walltime));
logf('  ------------------------------------------------\n');
logf('  Files saved:\n');
logf('    Data  : simulation_results.mat\n');
logf('    Plot  : summary_plot.png\n');
if generate_animation
    logf('    Video : bacteria_3d_motion.mp4\n');
end
logf('    Log   : simulation.log\n');
logf('================================================================\n');
logf('  Output units: time [s], length [um], force [pN], angle [rad]\n');
logf('================================================================\n\n');

% Close log file
fclose(log_fid);
fprintf('Log saved: %s\n', log_file);


%% ==================== Local Functions ====================

function dual_log(fid, varargin)
%DUAL_LOG  Write formatted message to both console and log file
    msg = sprintf(varargin{:});
    fprintf('%s', msg);           % Console
    fprintf(fid, '%s', msg);      % File
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
    value      = min(timeout_val, singularity_val);
    isterminal = 1;
    direction  = -1;
end
