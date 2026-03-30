%% batch_sweep_pili_site.m
%  Automated Parameter Sweep: Pilus Attachment Site
%
%  Description:
%    Sweeps the attachment site index (1-14) of a single linear pilus
%    while keeping its length and all body/fluid parameters fixed.
%    Different sites correspond to distinct spatial positions on the
%    bacterial body surface, allowing assessment of how pilus insertion
%    location affects the sinking dynamics.  Results are compared via
%    the steady-state sinking speed V_sink.
%
%    Key points:
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
cfg_base.case_name  = '';
cfg_base.minor_axis = 0.4e-6;
cfg_base.major_axis = 1.25e-6;
cfg_base.N_body     = 1000;
cfg_base.alpha_0    = 0;

cfg_base.seg_per_um = 15;

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
cfg_base.tlim   = 120;
cfg_base.dt     = 0.1;

cfg_base.enable_brown = false;
cfg_base.Temper       = 293.15;
cfg_base.dt_B         = cfg_base.dt / 2;

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
site_values = 1:14;     % All valid attachment site indices
n_runs      = numel(site_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);

fprintf('========================================\n');
fprintf('  Parameter Sweep: pili attachment site\n');
fprintf('  %d runs scheduled\n', n_runs);
fprintf('========================================\n\n');

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    pili_cfg_i(1).type      = 'linear';
    pili_cfg_i(1).length_um = 3;
    pili_cfg_i(1).site      = site_values(ii);

    cfg_i             = cfg_base;
    cfg_i.Pili_Matrix = parse_pili_config(pili_cfg_i);
    cfg_i.case_name   = sprintf('sweep-site-%02d', site_values(ii));

    fprintf('\n--- Run %d/%d: pili site = %d ---\n', ...
        ii, n_runs, site_values(ii));

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
        fprintf('  pili site = %2d  |  Steady V_sink = %.4e um/s\n', ...
            site_values(ii), mean(Vs(round(0.9 * numel(Vs)):end)));
    else
        fprintf('  pili site = %2d  |  ABORTED\n', site_values(ii));
    end
end
fprintf('========================================\n');
