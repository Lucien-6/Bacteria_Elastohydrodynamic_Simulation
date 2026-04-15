%% batch_sweep_seg_per_um.m
%  Automated Parameter Sweep: Pili Segments per Micron (seg_per_um)
%
%  Description:
%    Sweeps seg_per_um (discrete segments per micron along each pilus),
%    which sets the elastic rod resolution, while keeping body and fluid
%    parameters fixed.  Values < 5 trigger a validator warning (RSS accuracy).
%    Results compared via V_sink steady mean.
%
%    Key points:
%      - Edit seg_values below for the sweep grid (typically integers >= 5)
%      - auto_confirm = true; preview_initial = false
%      - generate_animation = true; animation_skip = 5
%      - enable_brown = false
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-30

close all force;
clear; clc;

%% ---- Define base configuration ----
cfg_base.case_name  = '';
cfg_base.minor_axis = 0.4e-6;
cfg_base.major_axis = 1.25e-6;
cfg_base.N_body     = 200;
cfg_base.alpha_0    = 0;

pili_cfg(1).type      = 'linear';
pili_cfg(1).length_um = 3;
pili_cfg(1).site      = 1;
cfg_base.Pili_Matrix  = parse_pili_config(pili_cfg);
cfg_base.seg_per_um   = 10;          % Overridden each run

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
cfg_base.tlim   = 360;
cfg_base.dt     = 0.1;

cfg_base.enable_brown = false;
cfg_base.Temper       = 293.15;
cfg_base.dt_B         = cfg_base.dt / 2;
cfg_base.fixman_midpoint  = false;

cfg_base.preview_initial    = false;
cfg_base.generate_animation = true;
cfg_base.animation_fps      = 30;
cfg_base.animation_skip     = 5;

cfg_base.auto_confirm = true;

cfg_base.flowfield_grid_spacing = 0.05;
cfg_base.flowfield_grid_range   = 10.0;
cfg_base.flowfield_quiver_skip  = 5;
cfg_base.flowfield_chunk_size   = 2000;

%% ---- Define sweep values ----
%  Segments per micron; >= 5 recommended (validator warning if lower)
seg_values = 10:5:35;
n_runs     = numel(seg_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);

fprintf('========================================\n');
fprintf('  Parameter Sweep: seg_per_um\n');
fprintf('  %d runs scheduled\n', n_runs);
fprintf('========================================\n\n');

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    sg = seg_values(ii);
    cfg_i             = cfg_base;
    cfg_i.seg_per_um  = sg;
    cfg_i.case_name   = sprintf('sweep-segperum-%d', sg);

    fprintf('\n--- Run %d/%d: seg_per_um = %g ---\n', ii, n_runs, sg);

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
        fprintf('  seg_per_um = %5g  |  Steady V_sink = %.4e um/s\n', ...
            seg_values(ii), mean(Vs(max(1, round(0.95 * numel(Vs))):end)));
    else
        fprintf('  seg_per_um = %5g  |  ABORTED\n', seg_values(ii));
    end
end
fprintf('========================================\n');
