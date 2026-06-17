%% main_coupled_simulation.m
%  Bacteria Elastohydrodynamic Simulation - Entry Point
%
%  Description:
%    Lightweight entry script that defines the simulation configuration
%    (either inline or from a JSON file) and delegates execution to
%    run_simulation().
%
%    To run with inline parameters:  edit the defaults below, then Run.
%    To run from a JSON config file: set config_file to the path.
%    To automate parameter sweeps:   see batch sweep/ folder (batch_sweep_alpha.m etc.).
%
%    CURRENT INLINE DEMO (intentional, non-default physics):
%      case_name = 'OnePili-Rotate-Active-Test'.  A single free-torsion
%      flagellum driven by a constant motor torque.  Buoyancy is OFF on
%      purpose: rho_bac = rho_fluid = 998.232 kg/m^3 isolates motor-driven
%      rotation and propulsion from gravitational sinking.  This is NOT
%      the documented project default - the standard sinking workflow uses
%      rho_bac = 1180 kg/m^3.  Likewise T_obs = 20 s and dt = 0.01 s are
%      tuned for the high-frequency motor diagnostic (the JSON loader
%      fallback defaults in load_sim_config.m are T_obs = 30 s, dt = 0.1 s).
%      Edit the rho_bac / dt / T_obs / cfg
%      lines below to reproduce the standard sinking-bacterium scenario.
%
%  Output unit convention:
%    Time [s] | Length [um] | Velocity [um/s] | Force [pN] | Angle [rad] | Angular velocity [rad/s]
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] Cortez et al. (2005), Physics of Fluids, 17, 031504.
%    [3] Cortez (2018), J. Comput. Phys. 375, 783-796.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Documentation release: 7.8.0 (2026-06-11) - see README.md / CHANGELOG.md

close all force;
clear; clc;

% Path policy (project root): main entry here + defensive addpath inside
% run_simulation(). Post-processing and batch sweep scripts keep their own
% addpath calls when run standalone - that layout is intentional.
proj_root = fileparts(mfilename('fullpath'));
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));

%% =====================================================================
%  CONFIGURATION SOURCE
%  =====================================================================
%  Priority: config_file (if non-empty) > inline defaults below.
%  To use external config: set config_file to a valid JSON path.
%  To use inline parameters: leave config_file = '' (default).
%
%  --- Geometry Modes ---
%  'parametric' (default): Use capsule body + analytic pili shapes.
%    Configure body (minor_axis, major_axis, N_body) and pili (Pili_Matrix)
%    via the body/pili sections below.
%
%  'custom': Import arbitrary 3D point clouds for both body and pili.
%    Custom body and custom pili MUST be used together (no mixing).
%    Use configs/config_custom_example.json (which references geometry/*.csv) as
%    a template, or set cfg.geometry_mode = 'custom' and fill the inline
%    custom-mode block below.
config_file = '';    % e.g. 'configs/config_custom_example.json' for the custom demo

if ~isempty(config_file) && isfile(config_file)
    %% --- External JSON configuration ---
    fprintf('Loading configuration from: %s\n', config_file);
    cfg = load_sim_config(config_file);
else
    %% --- Inline parameter definitions ---

    % ===== Geometry Mode =====
    % Set to 'parametric' (default) for capsule body + analytic pili shapes.
    % Set to 'custom' to import arbitrary body/pili point clouds.
    cfg.geometry_mode = 'parametric';

    if strcmp(cfg.geometry_mode, 'parametric')
        % ===== Parametric Mode: Capsule Body + Analytic Pili =====

        cfg.case_name = 'OnePili-Rotate-Active-Fast-Test';

        % -------- Bacterial Body Geometry --------
        cfg.minor_axis = 0.4e-6;              % [m]   Semi-minor axis (body radius)
        cfg.major_axis = 1.25e-6;             % [m]   Semi-major axis (half body length)
        cfg.N_body     = 300;                % [-]   Number of body surface points
        cfg.alpha_0    = 0;                   % [rad] Initial tilt angle in XZ plane

        % -------- Pili Configuration (struct array) --------
        %  See parse_pili_config.m for full field documentation.
        %  For a bare body with NO pili, use: pili_config = struct([]);
        pili_config(1).type        = 'Spiral';
        pili_config(1).length_um   = 4;
        pili_config(1).shape_params = [0.1, 1, 2*pi];
        pili_config(1).deflection_z = 0;
        pili_config(1).deflection_x = 0;
        pili_config(1).site        = 1;
        % pili_config(2).type        = 'linear';
        % pili_config(2).length_um   = 4;
        % pili_config(2).deflection_z = pi/2;
        % pili_config(2).deflection_x = pi/2;
        % pili_config(2).site        = 5;

        % -------- Per-pili: Intrinsic Curvature --------
        %  Options: 'frozen' (default, keep initial geometry), 'zero',
        %           struct('type','constant','kappa',[k1 k2 k3]),
        %           struct('type','helix','kappa_d2',val,'kappa_d3',val)
        pili_config(1).intrinsic_curvature = 'frozen';
        % pili_config(2).intrinsic_curvature = 'frozen';
        
        % -------- Per-pili: Active Moment --------
        %  Options: struct('type','off')  (default)
        % pili_config(1).active_moment = struct('type','traveling_wave', ...
        %               'amplitude_Nm', 3e-19, 'omega_Hz', 50, ...
        %               'wavenumber', 2*pi, 'phase_rad', 0, ... 
        %               'd1_amp', 1, 'd2_amp', 0, 'd3_amp', 0);
        % pili_config(2).active_moment = struct('type','traveling_wave', ...
        %               'amplitude_Nm', 2e-19, 'omega_Hz', 50, ...
        %               'wavenumber', 2*pi, 'phase_rad', 0, ... 
        %               'd1_amp', 1, 'd2_amp', 0, 'd3_amp', 0);

        %  Custom:     struct('type','custom', 'fun', @my_moment_func)
        % pili_config(1).active_moment = struct('type', 'off');

        % -------- Per-pili: Root Coupling & Motor Torque --------
        %  root_coupling: 'clamped' (default) | 'free_torsion'
        %    clamped     - Position + all orientations locked to body (standard)
        %    free_torsion - Position + bending locked, twist (psi_1) free
        %  motor_torque: only meaningful for free_torsion pili
        %    struct('type','off')  (default)
        %    struct('type','constant','torque_pN_nm',2000,'max_omega_Hz',300)
        %      - Constant motor torque with linear speed limiting
        %      - torque_pN_nm: output torque [pN*nm] (positive=CCW, negative=CW)
        %      - max_omega_Hz: stall speed [Hz]
        % pili_config(1).root_coupling = 'clamped';
        pili_config(1).root_coupling = 'free_torsion';
        pili_config(1).motor_torque  = struct('type','constant', ...
            'torque_pN_nm', -1, 'max_omega_Hz', 500);

        [cfg.Pili_Matrix, cfg.pili_ext] = parse_pili_config(pili_config);

        % Legacy raw-matrix alternative (uncomment to use instead):
        %  cfg.Pili_Matrix = [1; 3; 0; 0; NaN; NaN; NaN; 1];
        %    Row 1: type, Row 2: length [um], Row 3: deflection Z [rad],
        %    Row 4: deflection X [rad], Rows 5-7: shape params, Row 8: site

        % -------- Background Flow --------
        %  Options: 'off' (default), 'uniform', 'shear_z'
        cfg.background_flow.type           = 'off';
        % cfg.background_flow.type         = 'shear_z';
        % cfg.background_flow.U_max        = 1e-6;       % [m/s]
        % cfg.background_flow.direction    = [1, 0, 0];
        % cfg.background_flow.shear_axis   = [0, 0, 1];
        % cfg.background_flow.shear_height_m = 3e-6;     % [m]

    else
        % ===== Custom Mode: Import arbitrary body/pili point clouds =====
        %  Custom body and custom pili MUST be used together (no mixing).
        %  The block below mirrors config_custom_example.json and uses the
        %  pre-shipped geometry/example_*.csv files so the example runs
        %  out of the box.  Edit paths / parameters for your own geometry.
        %  IMPORTANT: user-supplied pili centerlines MUST be pre-discretised
        %  with EQUAL-LENGTH segments; this is enforced at validation.

        cfg.case_name = 'Custom-Geometry-Run';

        % Placeholder parametric fields (required by run_simulation.m unpacking)
        cfg.Pili_Matrix = zeros(8, 0);
        cfg.pili_ext    = struct('intrinsic_curvature', {}, 'active_moment', {}, ...
                                     'root_coupling', {}, 'motor_torque', {});

        % -------- Custom body --------
        cfg.custom_geometry.body_points_file = 'geometry/example_body_capsule.csv';
        cfg.custom_geometry.L_0_override     = [];    % empty = auto (volume-based)
        cfg.custom_geometry.L_0_suggest_mode = 'volume';  % 'volume' | 'extent'
        cfg.eps_body                          = [];   % RSM regularisation override

        % -------- Initial orientation (quaternion / axis_angle / euler_zyx) --------
        %   Omit or set to [] to fall back on cfg.alpha_0 (legacy Y-axis tilt).
        cfg.custom_geometry.initial_orientation = struct( ...
            'euler_zyx', [0.0, 0.0, 0.0] );           % [yaw, pitch, roll] rad
        cfg.initial_orientation = cfg.custom_geometry.initial_orientation;  % lift to top level for run_simulation
        cfg.alpha_0 = 0;                               % ignored when initial_orientation set

        % -------- Attachment sites --------
        cfg.custom_geometry.attachment_sites(1).position  = [-1.25e-6; 0; 0];
        cfg.custom_geometry.attachment_sites(1).direction = [-1; 0; 0];
        cfg.custom_geometry.attachment_sites(1).label     = 'left_pole';
        cfg.custom_geometry.attachment_sites(2).position  = [ 1.25e-6; 0; 0];
        cfg.custom_geometry.attachment_sites(2).direction = [ 1; 0; 0];
        cfg.custom_geometry.attachment_sites(2).label     = 'right_pole';

        % -------- Custom pili (each centerline MUST have uniform segments) --------
        cfg.custom_geometry.pili(1).points_file        = 'geometry/example_pilus_linear.csv';
        cfg.custom_geometry.pili(1).site_index         = 1;
        cfg.custom_geometry.pili(1).deflection_z       = 0;
        cfg.custom_geometry.pili(1).deflection_x       = 0;
        cfg.custom_geometry.pili(1).intrinsic_curvature = 'frozen';
        cfg.custom_geometry.pili(1).active_moment       = struct('type', 'off');
        cfg.custom_geometry.pili(1).root_coupling       = 'clamped';
        cfg.custom_geometry.pili(1).motor_torque        = struct('type', 'off');

        cfg.custom_geometry.pili(2).points_file        = 'geometry/example_pilus_helix.csv';
        cfg.custom_geometry.pili(2).site_index         = 2;
        cfg.custom_geometry.pili(2).deflection_z       = 0.1;
        cfg.custom_geometry.pili(2).deflection_x       = 0;
        cfg.custom_geometry.pili(2).intrinsic_curvature = 'frozen';
        cfg.custom_geometry.pili(2).active_moment       = struct('type', 'off');
        cfg.custom_geometry.pili(2).root_coupling       = 'clamped';
        cfg.custom_geometry.pili(2).motor_torque        = struct('type', 'off');

        % Populate top-level cfg.pili_ext from per-pilus custom config so
        % the extended fields propagate identically to the parametric path.
        for kp = 1:numel(cfg.custom_geometry.pili)
            cfg.pili_ext(kp).intrinsic_curvature = cfg.custom_geometry.pili(kp).intrinsic_curvature;
            cfg.pili_ext(kp).active_moment       = cfg.custom_geometry.pili(kp).active_moment;
            cfg.pili_ext(kp).root_coupling       = cfg.custom_geometry.pili(kp).root_coupling;
            cfg.pili_ext(kp).motor_torque        = cfg.custom_geometry.pili(kp).motor_torque;
        end

        cfg.background_flow.type           = 'off';
    end

    % -------- Pili Discretisation (RSS, Cortez 2018) --------
    cfg.seg_per_um = 10;                  % [seg/um] RSS segment density

    % -------- Material Properties --------
    cfg.E       = 3.279e8;                % [Pa]   Young's modulus of pili
    cfg.r_pili  = 2e-8;                   % [m]    Pili radius
    cfg.sigma_p = 0.2;                    % [-]    Poisson's ratio

    % -------- Fluid Properties --------
    cfg.mu        = 1e-3;                 % [Pa*s]   Dynamic viscosity (water)
    % NOTE: rho_bac = rho_fluid disables buoyancy; intentional for the
    % motor-driven flagella demo above.  Use cfg.rho_bac = 1180 (kg/m^3,
    % load_sim_config.m default) to recover the standard sinking-bacterium physics.
    cfg.rho_bac   = 998.232;              % [kg/m^3] Bacteria density (= fluid: zero gravity)
    cfg.rho_fluid = 998.232;              % [kg/m^3] Fluid density

    % -------- Simulation Time & Numerics --------
    % NOTE: T_obs / dt below are tuned for the high-frequency motor demo;
    % the JSON loader fallback defaults are T_obs = 30 s, dt = 0.1 s for sinking studies.
    cfg.T_obs  = 5;                    % [s]   Total observation time
    cfg.AbsTol = 1e-4;                    % ODE absolute tolerance
    cfg.RelTol = 1e-4;                    % ODE relative tolerance
    cfg.delta  = pi / 20;                 % [rad] Euler angle singularity threshold
    cfg.tlim   = 3600;                    % [s]   ODE timeout limit
    cfg.dt     = 1e-2;                    % [s]   Output time step

    % -------- Brownian Motion --------
    cfg.enable_brown     = false;          % true = enable Brownian motion
    cfg.Temper           = 293.15;        % [K]   Fluid temperature (20 deg C)
    cfg.dt_B             = cfg.dt / 2;    % [s]   Brownian time step
    % Brownian thermal drift is handled by the strict Ito finite-difference
    % correction in run_simulation.m.
    % Optional: cfg.brown_rng_seed = positive integer for fixed twister seed; omit = MMddHHmmss at run start

    % -------- Preview & Animation --------
    cfg.preview_initial    = true;
    cfg.generate_animation = true;
    cfg.animation_fps      = 30;
    cfg.animation_skip     = 1;

    % -------- Automation --------
    cfg.auto_confirm = false;             % true = skip all interactive dialogs

    % -------- Flow Field Visualisation --------
    cfg.flowfield_grid_spacing = 0.05;    % [um]
    cfg.flowfield_grid_range   = 10.0;    % [um]
    cfg.flowfield_quiver_skip  = 5;
    cfg.flowfield_chunk_size   = 2000;
end

%% =====================================================================
%  RUN SIMULATION
%  =====================================================================

[results, params, output_dir] = run_simulation(cfg);

if isempty(results)
    fprintf('Simulation was aborted.\n');
end
