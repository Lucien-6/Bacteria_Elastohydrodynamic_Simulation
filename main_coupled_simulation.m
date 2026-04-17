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

close all force;
clear; clc;

%% =====================================================================
%  CONFIGURATION SOURCE
%  =====================================================================
%  Priority: config_file (if non-empty) > inline defaults below.
%  To use external config: set config_file to a valid JSON path.
%  To use inline parameters: leave config_file = '' (default).
config_file = '';    % e.g. 'configs/my_case.json'

if ~isempty(config_file) && isfile(config_file)
    %% --- External JSON configuration ---
    fprintf('Loading configuration from: %s\n', config_file);
    cfg = load_sim_config(config_file);
else
    %% --- Inline parameter definitions ---
    cfg.case_name = 'TwoPili-Active-HighFre-Test-new1';

    % -------- Bacterial Body Geometry --------
    cfg.minor_axis = 0.4e-6;              % [m]   Semi-minor axis (body radius)
    cfg.major_axis = 1.25e-6;             % [m]   Semi-major axis (half body length)
    cfg.N_body     = 200;                % [-]   Number of body surface points
    cfg.alpha_0    = 0;                   % [rad] Initial tilt angle in XZ plane

    % -------- Pili Configuration (struct array) --------
    %  See parse_pili_config.m for full field documentation.
    %  For a bare body with NO pili, use: pili_config = struct([]);
    pili_config(1).type        = 'linear';
    pili_config(1).length_um   = 3;
    pili_config(1).deflection_z = pi/2;
    pili_config(1).deflection_x = pi/2;
    pili_config(1).site        = 3;
    pili_config(2).type        = 'linear';
    pili_config(2).length_um   = 3;
    pili_config(2).deflection_z = pi/2;
    pili_config(2).deflection_x = pi/2;
    pili_config(2).site        = 5;

    % -------- Per-pili: Intrinsic Curvature --------
    %  Options: 'frozen' (default, keep initial geometry), 'zero',
    %           struct('type','constant','kappa',[k1 k2 k3]),
    %           struct('type','helix','kappa_d2',val,'kappa_d3',val)
    pili_config(1).intrinsic_curvature = 'frozen';
    pili_config(2).intrinsic_curvature = 'frozen';
    
    % -------- Per-pili: Active Moment --------
    %  Options: struct('type','off')  (default)
    pili_config(1).active_moment = struct('type','traveling_wave', ...
                  'amplitude_Nm', 2e-19, 'omega_Hz', 50, ...
                  'wavenumber', 2*pi, 'phase_rad', 0, ... 
                  'd1_amp', 1, 'd2_amp', 0, 'd3_amp', 0);
    pili_config(2).active_moment = struct('type','traveling_wave', ...
                  'amplitude_Nm', 2e-19, 'omega_Hz', 50, ...
                  'wavenumber', 2*pi, 'phase_rad', pi, ... 
                  'd1_amp', 1, 'd2_amp', 0, 'd3_amp', 0);

    %  Custom:     struct('type','custom', 'fun', @my_moment_func)
    % pili_config(1).active_moment = struct('type', 'off');

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
    cfg.T_obs  = 30.0;                    % [s]   Total observation time
    cfg.AbsTol = 1e-4;                    % ODE absolute tolerance
    cfg.RelTol = 1e-4;                    % ODE relative tolerance
    cfg.delta  = pi / 20;                 % [rad] Euler angle singularity threshold
    cfg.tlim   = 360;                     % [s]   ODE timeout limit
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
