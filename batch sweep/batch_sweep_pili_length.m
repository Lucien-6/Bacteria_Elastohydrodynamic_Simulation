%% batch_sweep_pili_length.m
%  Automated Parameter Sweep: Pilus Length (single or multiple pili)
%
%  Description:
%    Sweeps the contour length [um] shared by all linear pili while keeping
%    other parameters fixed.  When N_pili > 1, each run uses N_pili identical
%    filaments at consecutive attachment sites 1..N_pili (all same length).
%    The number of discrete segments per filament scales with length via
%    seg_per_um.  Results are compared via the steady-state sinking speed V_sink.
%
%    Key points:
%      - Set N_pili (1-14) below; default 1 recovers the original single-pilus sweep
%      - auto_confirm = true  to suppress all interactive dialogs
%      - preview_initial = false  to skip the pose preview figure
%      - generate_animation = false  to reduce per-run wall time
%      - Each run produces its own output directory with full logs and data
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-18
%  Updated: 2026-03-30 — multi-pili support (equal length per run, sites 1:N_pili)

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

%% ---- Multi-pili: equal length per run (sites 1 .. N_pili) ----
N_pili = 1;   % integer 1-14; all pili share length_um from the sweep each run

if N_pili < 1 || N_pili > 14 || floor(N_pili) ~= N_pili
    error('batch_sweep_pili_length:InvalidNPili', ...
        'N_pili must be an integer in 1..14 (attachment sites 1..N_pili).');
end

%% ---- Define sweep values ----
length_values = [1, 2, 3, 4, 5, 6];    % Pilus contour length [um], all pili equal
n_runs        = numel(length_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);

fprintf('========================================\n');
fprintf('  Parameter Sweep: pili length (N_pili = %d, equal length per run)\n', N_pili);
fprintf('  %d runs scheduled\n', n_runs);
fprintf('========================================\n\n');

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    L = length_values(ii);
    % N_pili linear pili, same length, consecutive sites (same pattern as batch_sweep_pili_count.m)
    pili_cfg_i = struct('type', 'linear', 'length_um', L, 'site', num2cell(1:N_pili));

    cfg_i             = cfg_base;
    cfg_i.Pili_Matrix = parse_pili_config(pili_cfg_i);
    cfg_i.case_name   = sprintf('sweep_plen_%02d_um_np%02d', L, N_pili);

    fprintf('\n--- Run %d/%d: length = %d um, N_pili = %d ---\n', ...
        ii, n_runs, L, N_pili);

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
        fprintf('  length = %2d um  (N_pili=%d)  |  Steady V_sink = %.4e um/s\n', ...
            length_values(ii), N_pili, mean(Vs(round(0.9 * numel(Vs)):end)));
    else
        fprintf('  length = %2d um  (N_pili=%d)  |  ABORTED\n', ...
            length_values(ii), N_pili);
    end
end
fprintf('========================================\n');
