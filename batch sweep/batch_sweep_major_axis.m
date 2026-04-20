%% batch_sweep_major_axis.m
%  Automated Parameter Sweep: Bacterial Body Major Axis Length
%
%  Description:
%    Sweeps the semi-major axis of the bacterial ellipsoidal body while
%    keeping all other parameters fixed.  The aspect ratio changes with
%    each run because minor_axis is held constant.  Model is body-only
%    (no pili): empty Pili_Matrix.  Results compared via V_sink steady mean.
%
%    Key points:
%      - No pili: Pili_Matrix = zeros(8, 0)
%      - auto_confirm = true  to suppress all interactive dialogs
%      - preview_initial = false  to skip the pose preview figure
%      - generate_animation = true; animation_skip = 5
%      - Each run produces its own output directory with full logs and data
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-18

close all force;
clear; clc;

%% ---- Define base configuration ----
cfg_base.case_name  = '';           % Will be auto-set per run below
cfg_base.minor_axis = 0.4e-6;
cfg_base.major_axis = 1.25e-6;     % Overridden in the loop
cfg_base.N_body     = 200;
cfg_base.alpha_0    = 0;

cfg_base.Pili_Matrix  = zeros(8, 0);   % No pili (body only)
cfg_base.seg_per_um   = 10;

cfg_base.E       = 3.279e8;
cfg_base.r_pili  = 3e-9;
cfg_base.sigma_p = 0.2;

cfg_base.mu        = 1e-3;
cfg_base.rho_bac   = 1180.0;
cfg_base.rho_fluid = 998.232;

cfg_base.T_obs  = 200.0;
cfg_base.AbsTol = 1e-4;
cfg_base.RelTol = 1e-4;
cfg_base.delta  = pi / 20;
cfg_base.tlim   = 3600;
cfg_base.dt     = 0.1;

cfg_base.enable_brown = false;
cfg_base.Temper       = 293.15;
cfg_base.dt_B         = cfg_base.dt / 2;
cfg_base.fixman_midpoint  = false;

cfg_base.preview_initial    = false;    % No interactive preview
cfg_base.generate_animation = true;
cfg_base.animation_fps      = 30;
cfg_base.animation_skip     = 5;

cfg_base.auto_confirm = true;           % Headless mode

cfg_base.flowfield_grid_spacing = 0.05;
cfg_base.flowfield_grid_range   = 10.0;
cfg_base.flowfield_quiver_skip  = 5;
cfg_base.flowfield_chunk_size   = 2000;

%% ---- Define sweep values ----
% Semi-major axis [m]; minor_axis fixed at 0.4 um; grid 0.4:0.2:1.6 um
major_values = (0.4:0.2:1.6) * 1e-6;
n_runs       = numel(major_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);

fprintf('========================================\n');
fprintf('  Parameter Sweep: major_axis (body only, no pili)\n');
fprintf('  %d runs scheduled\n', n_runs);
fprintf('========================================\n\n');

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    cfg_i = cfg_base;
    cfg_i.major_axis = major_values(ii);
    cfg_i.case_name  = sprintf('sweep-major-%04.0f-nm', major_values(ii) * 1e9);

    fprintf('\n--- Run %d/%d: major_axis = %.2f um ---\n', ...
        ii, n_runs, major_values(ii) * 1e6);

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
        fprintf('  major_axis = %5.2f um  |  Steady V_sink = %.4e um/s\n', ...
            major_values(ii) * 1e6, mean(Vs(max(1, round(0.95 * numel(Vs))):end)));
    else
        fprintf('  major_axis = %5.2f um  |  ABORTED\n', major_values(ii) * 1e6);
    end
end
fprintf('========================================\n');
