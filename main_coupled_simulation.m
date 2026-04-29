%% main_coupled_simulation.m
%  Coupled Body-Pili Elastohydrodynamic Simulation — Entry Point
%
%  Description:
%    Lightweight entry script that defines the simulation configuration
%    (either inline or from a JSON file) and delegates execution to
%    run_simulation().
%
%    To run with inline parameters:  edit the defaults below, then Run.
%    To run from a JSON config file: set config_file to the path.
%    To automate parameter sweeps:   see batch sweep/ folder (batch_sweep_alpha.m etc.).
%    Default inline demo: TwoPili-Active-HighFre-Test-Symmetry (edit the cfg block below).
%
%  Output unit convention:
%    Time [s] | Length [um] | Velocity [um/s] | Force [pN] | Angle [rad]
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] Cortez et al. (2005), Physics of Fluids, 17, 031504.
%    [3] Cortez (2018), J. Comput. Phys. 375, 783-796.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Documentation release: 7.0.0 (2026-04-29) — see README.md / CHANGELOG.md

close all force;
clear; clc;

% Path policy (project root): main entry here + defensive addpath inside
% run_simulation(). Post-processing and batch sweep scripts keep their own
% addpath calls when run standalone — that layout is intentional.
proj_root = fileparts(mfilename('fullpath'));
addpath(proj_root);

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
%    Use config_custom_example.json as a template, or set
%    cfg.geometry_mode = 'custom' and cfg.custom_geometry below.
config_file = '';    % e.g. 'configs/my_case.json' or 'config_custom_example.json'

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

        cfg.case_name = 'OnePili-Active-Test';

        % -------- Bacterial Body Geometry --------
        cfg.minor_axis = 0.4e-6;              % [m]   Semi-minor axis (body radius)
        cfg.major_axis = 1.25e-6;             % [m]   Semi-major axis (half body length)
        cfg.N_body     = 300;                % [-]   Number of body surface points
        cfg.alpha_0    = 0;                   % [rad] Initial tilt angle in XZ plane

        % -------- Pili Configuration (struct array) --------
        %  See parse_pili_config.m for full field documentation.
        %  For a bare body with NO pili, use: pili_config = struct([]);
        pili_config(1).type        = 'linear';
        pili_config(1).length_um   = 4;
        pili_config(1).deflection_z = 0;
        pili_config(1).deflection_x = 0;
        pili_config(1).site        = 5;
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
        pili_config(1).active_moment = struct('type','traveling_wave', ...
                      'amplitude_Nm', 3e-19, 'omega_Hz', 50, ...
                      'wavenumber', 2*pi, 'phase_rad', 0, ... 
                      'd1_amp', 1, 'd2_amp', 0, 'd3_amp', 0);
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
        pili_config(1).root_coupling = 'clamped';
        % pili_config(1).root_coupling = 'free_torsion';
        % pili_config(1).motor_torque  = struct('type','constant', ...
        %     'torque_pN_nm', 2000, 'max_omega_Hz', 300);

        [cfg.Pili_Matrix, cfg.pili_ext] = parse_pili_config(pili_config);

        % Legacy raw-matrix alternative (uncomment to use instead):
        %  cfg.Pili_Matrix = [1; 3; 0; 0; NaN; NaN; NaN; 1];

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
        %  Uncomment and edit the block below to configure custom geometry:

        cfg.case_name = 'Custom-Geometry-Run';
        cfg.alpha_0   = 0;                    % [rad] Initial tilt angle in XZ plane
        % cfg.eps_body = [];                  % [L0] RSM regularisation override (empty = auto)

        % Placeholder parametric fields (required by run_simulation.m unpacking)
        cfg.Pili_Matrix = zeros(8, 0);
        cfg.pili_ext    = struct('intrinsic_curvature', {}, 'active_moment', {}, ...
                                     'root_coupling', {}, 'motor_torque', {});

        % cfg.custom_geometry.body_points_file = 'geometry/my_body_points.csv';
        % cfg.custom_geometry.L_0_override = 1.0e-6;
        % cfg.custom_geometry.attachment_sites(1).position = [-1.0e-6, 0, 0];
        % cfg.custom_geometry.attachment_sites(1).direction = [-1, 0, 0];
        % cfg.custom_geometry.attachment_sites(1).label = 'left_pole';
        % cfg.custom_geometry.pili(1).points_file = 'geometry/my_pilus_centerline.csv';
        % cfg.custom_geometry.pili(1).site_index = 1;
        % cfg.custom_geometry.pili(1).deflection_z = 0;
        % cfg.custom_geometry.pili(1).deflection_x = 0;

        cfg.background_flow.type           = 'off';
    end

    % -------- Pili Discretisation (RSS, Cortez 2018) --------
    cfg.seg_per_um = 10;                  % [seg/um] RSS segment density

    % -------- Material Properties --------
    cfg.E       = 3.279e8;                % [Pa]   Young's modulus of pili
    cfg.r_pili  = 3e-9;                   % [m]    Pili radius
    cfg.sigma_p = 0.2;                    % [-]    Poisson's ratio

    % -------- Fluid Properties --------
    cfg.mu        = 1e-3;                 % [Pa*s]   Dynamic viscosity (water)
    cfg.rho_bac   = 998.232;              % [kg/m^3] Bacteria density
    cfg.rho_fluid = 998.232;             % [kg/m^3] Fluid density

    % -------- Simulation Time & Numerics --------
    cfg.T_obs  = 60.0;                    % [s]   Total observation time
    cfg.AbsTol = 1e-4;                    % ODE absolute tolerance
    cfg.RelTol = 1e-4;                    % ODE relative tolerance
    cfg.delta  = pi / 20;                 % [rad] Euler angle singularity threshold
    cfg.tlim   = 3600;                    % [s]   ODE timeout limit
    cfg.dt     = 1e-2;                    % [s]   Output time step

    % -------- Brownian Motion --------
    cfg.enable_brown     = false;          % true = enable Brownian motion
    cfg.Temper           = 293.15;        % [K]   Fluid temperature (20 deg C)
    cfg.dt_B             = cfg.dt / 2;    % [s]   Brownian time step
    cfg.fixman_midpoint  = false;         % true = Fixman midpoint drift correction (2x cost)
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
