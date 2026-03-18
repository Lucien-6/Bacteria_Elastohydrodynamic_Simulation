%% batch_sweep_pili_count.m
%  Automated Parameter Sweep: Number of Pili
%
%  Description:
%    Sweeps the total number of pili attached to the bacterial body (1-6).
%    Each run adds one additional linear pilus of identical geometry
%    (length 3 um), assigning consecutive attachment sites 1, 2, ..., N.
%    All body, fluid and mechanical parameters remain fixed.  Results are
%    compared via the steady-state sinking speed V_sink.
%
%    Key points:
%      - auto_confirm = true  to suppress all interactive dialogs
%      - preview_initial = false  to skip the pose preview figure
%      - generate_animation = false  to reduce per-run wall time
%      - Each run produces its own output directory with full logs and data
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-18

close all force;
clear; clc;

%% ---- Define base configuration ----
cfg_base.case_name  = '';
cfg_base.minor_axis = 0.4e-6;
cfg_base.major_axis = 1.25e-6;
cfg_base.N_body     = 1000;
cfg_base.alpha_0    = 0;

cfg_base.seg_per_um = 15;

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

cfg_base.preview_initial    = false;
cfg_base.generate_animation = false;
cfg_base.animation_fps      = 30;
cfg_base.animation_skip     = 1;

cfg_base.auto_confirm = true;

cfg_base.flowfield_grid_spacing = 0.05;
cfg_base.flowfield_grid_range   = 10.0;
cfg_base.flowfield_quiver_skip  = 5;
cfg_base.flowfield_chunk_size   = 2000;

%% ---- Define sweep values ----
count_values = 1:6;     % Number of pili per run (sites 1:N are used)
n_runs       = numel(count_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);

fprintf('========================================\n');
fprintf('  Parameter Sweep: pili count\n');
fprintf('  %d runs scheduled\n', n_runs);
fprintf('========================================\n\n');

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    N = count_values(ii);

    % Build pili struct array: N identical pili at consecutive sites
    pili_cfg_i = struct('type', 'linear', 'length_um', 3, 'site', num2cell(1:N));

    cfg_i             = cfg_base;
    cfg_i.Pili_Matrix = parse_pili_config(pili_cfg_i);
    cfg_i.case_name   = sprintf('sweep_pcount_%02d', N);

    fprintf('\n--- Run %d/%d: pili count = %d ---\n', ii, n_runs, N);

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
        fprintf('  pili count = %d  |  Steady V_sink = %.4e um/s\n', ...
            count_values(ii), mean(Vs(round(0.9 * numel(Vs)):end)));
    else
        fprintf('  pili count = %d  |  ABORTED\n', count_values(ii));
    end
end
fprintf('========================================\n');
