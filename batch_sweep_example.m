%% batch_sweep_example.m
%  Example: Automated Parameter Sweep using run_simulation()
%
%  Description:
%    Demonstrates how to programmatically sweep a parameter (here: initial
%    tilt angle alpha_0) using run_simulation().  Each run is assigned a
%    unique case name, and results are collected for cross-comparison.
%
%    Key points:
%      - auto_confirm = true  to suppress all interactive dialogs
%      - preview_initial = false  to skip the pose preview figure
%      - generate_animation = false  to reduce per-run wall time
%      - Each run produces its own output directory with full logs and data
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17

close all force;
clear; clc;

%% ---- Define base configuration ----
cfg_base.case_name  = '';           % Will be auto-set per run below
cfg_base.minor_axis = 0.4e-6;
cfg_base.major_axis = 1.25e-6;
cfg_base.N_body     = 1000;
cfg_base.alpha_0    = 0;            % Overridden in the loop

pili_cfg(1).type      = 'linear';
pili_cfg(1).length_um = 3;
pili_cfg(1).site      = 1;
cfg_base.Pili_Matrix  = parse_pili_config(pili_cfg);
cfg_base.seg_per_um   = 15;

cfg_base.E       = 1.768e7;
cfg_base.r_pili  = 3e-9;
cfg_base.sigma_p = 0.2;

cfg_base.mu        = 1e-3;
cfg_base.rho_bac   = 1180.0;
cfg_base.rho_fluid = 998.232;

cfg_base.T_obs  = 30.0;
cfg_base.AbsTol = 1e-4;
cfg_base.RelTol = 1e-4;
cfg_base.delta  = pi / 20;
cfg_base.tlim   = 120;
cfg_base.dt     = 0.1;

cfg_base.enable_brown = true;
cfg_base.Temper       = 293.15;
cfg_base.dt_B         = cfg_base.dt / 2;

cfg_base.preview_initial    = false;    % No interactive preview
cfg_base.generate_animation = false;    % Skip animation for speed
cfg_base.animation_fps      = 30;
cfg_base.animation_skip     = 1;

cfg_base.auto_confirm = true;           % Headless mode

cfg_base.flowfield_grid_spacing = 0.05;
cfg_base.flowfield_grid_range   = 10.0;
cfg_base.flowfield_quiver_skip  = 5;
cfg_base.flowfield_chunk_size   = 2000;

%% ---- Define sweep values ----
alpha_values = deg2rad([0, 15, 30, 45, 60, 75, 90]);
n_runs       = numel(alpha_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);

fprintf('========================================\n');
fprintf('  Parameter Sweep: alpha_0\n');
fprintf('  %d runs scheduled\n', n_runs);
fprintf('========================================\n\n');

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    cfg_i = cfg_base;
    cfg_i.alpha_0   = alpha_values(ii);
    cfg_i.case_name = sprintf('sweep_alpha_%02d_deg', round(rad2deg(alpha_values(ii))));

    fprintf('\n--- Run %d/%d: alpha_0 = %.1f deg ---\n', ...
        ii, n_runs, rad2deg(alpha_values(ii)));

    [res_i, ~, dir_i] = run_simulation(cfg_i);
    sweep_results{ii} = res_i;
    sweep_dirs{ii}    = dir_i;
end

t_sweep_total = toc(t_sweep_start);

%% ---- Summary ----
fprintf('\n========================================\n');
fprintf('  Sweep Complete  (%d runs, %.1f s total)\n', n_runs, t_sweep_total);
fprintf('========================================\n');

for ii = 1:n_runs
    if ~isempty(sweep_results{ii})
        Vs = sweep_results{ii}.V_sink;
        fprintf('  alpha_0 = %5.1f deg  |  Steady V_sink = %.4e um/s\n', ...
            rad2deg(alpha_values(ii)), mean(Vs(round(0.9*numel(Vs)):end)));
    else
        fprintf('  alpha_0 = %5.1f deg  |  ABORTED\n', rad2deg(alpha_values(ii)));
    end
end
fprintf('========================================\n');
