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
    cfg.case_name = 'TwoPili-S28-L3-Brown-New';

    % -------- Bacterial Body Geometry --------
    cfg.minor_axis = 0.4e-6;              % [m]   Semi-minor axis (body radius)
    cfg.major_axis = 1.25e-6;             % [m]   Semi-major axis (half body length)
    cfg.N_body     = 500;                % [-]   Number of body surface points
    cfg.alpha_0    = 0;                   % [rad] Initial tilt angle in XZ plane

    % -------- Pili Configuration (struct array) --------
    %  See parse_pili_config.m for full field documentation.
    %  For a bare body with NO pili, use: pili_config = struct([]);
    pili_config(1).type        = 'linear';
    pili_config(1).length_um   = 3;
    pili_config(1).deflection_z = 0;
    pili_config(1).deflection_x = 0;
    pili_config(1).site        = 2;
    pili_config(2).type        = 'linear';
    pili_config(2).length_um   = 3;
    pili_config(2).deflection_z = 0;
    pili_config(2).deflection_x = 0;
    pili_config(2).site        = 8;
    cfg.Pili_Matrix = parse_pili_config(pili_config);

    % Legacy raw-matrix alternative (uncomment to use instead):
    %  cfg.Pili_Matrix = [1; 3; 0; 0; NaN; NaN; NaN; 1];

    % -------- Pili Discretisation (RSS, Cortez 2018) --------
    cfg.seg_per_um = 15;                  % [seg/um] RSS segment density

    % -------- Material Properties --------
    cfg.E       = 3.279e8;                % [Pa]   Young's modulus of pili
    cfg.r_pili  = 3e-9;                   % [m]    Pili radius
    cfg.sigma_p = 0.2;                    % [-]    Poisson's ratio

    % -------- Fluid Properties --------
    cfg.mu        = 1e-3;                 % [Pa*s]   Dynamic viscosity (water)
    cfg.rho_bac   = 1180.0;              % [kg/m^3] Bacteria density
    cfg.rho_fluid = 998.232;             % [kg/m^3] Fluid density

    % -------- Simulation Time & Numerics --------
    cfg.T_obs  = 30.0;                    % [s]   Total observation time
    cfg.AbsTol = 1e-4;                    % ODE absolute tolerance
    cfg.RelTol = 1e-4;                    % ODE relative tolerance
    cfg.delta  = pi / 20;                 % [rad] Euler angle singularity threshold
    cfg.tlim   = 120;                     % [s]   ODE timeout limit
    cfg.dt     = 1e-1;                    % [s]   Output time step

    % -------- Brownian Motion --------
    cfg.enable_brown     = true;          % true = enable Brownian motion
    cfg.Temper           = 293.15;        % [K]   Fluid temperature (20 deg C)
    cfg.dt_B             = cfg.dt / 2;    % [s]   Brownian time step
    cfg.fixman_midpoint  = false;         % true = Fixman midpoint drift correction (2x cost)

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
