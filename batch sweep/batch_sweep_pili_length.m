%% batch_sweep_pili_length.m
%  Automated Parameter Sweep: Pilus Length (single or multiple pili)
%
%  Description:
%    Sweeps the contour length [um] shared by all linear pili while keeping
%    other parameters fixed.  Each pilus uses an attachment site index from
%    site_list (1-14); sites may be chosen freely (e.g. [1 6 11] or 1:3).
%    All pili in a run share the same length_um; segment count scales with
%    length via seg_per_um.  Results are compared via steady-state V_sink.
%
%    Key points:
%      - Set site_list below (one entry per pilus); e.g. 1 = original single
%        pilus at site 1; [1 2 3 6] = four pili at those sites
%      - Duplicate sites are allowed but trigger init_coupled_params warning
%      - auto_confirm = true  to suppress all interactive dialogs
%      - preview_initial = false  to skip the pose preview figure
%      - generate_animation = true; animation_skip = 5
%      - Each run produces its own output directory with full logs and data
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-18
%  Updated: 2026-03-30 — site_list for arbitrary per-pilus attachment sites
%           2026-04-07 — cfg_base.fixman_midpoint

close all force;
clear; clc;

%% ---- Define base configuration ----
cfg_base.case_name  = '';
cfg_base.minor_axis = 0.4e-6;
cfg_base.major_axis = 1.25e-6;
cfg_base.N_body     = 200;
cfg_base.alpha_0    = 0;

cfg_base.seg_per_um = 10;

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

%% ---- Per-pilus attachment sites (indices 1-14, arbitrary order) ----
%  One integer per pilus; length(site_list) = number of pili.
%  Examples:  1          -> single pilus at site 1
%             1:6        -> six pili at consecutive sites 1..6
%             [1 6 11]   -> three pili at left pole, right y+, middle z-
site_list = 1;

site_list = site_list(:).';   % row vector
n_pili    = numel(site_list);
if n_pili < 1
    error('batch_sweep_pili_length:EmptySiteList', ...
        'site_list must contain at least one attachment site index.');
end
if any(site_list < 1 | site_list > 14 | floor(site_list) ~= site_list)
    error('batch_sweep_pili_length:InvalidSite', ...
        'Each site must be an integer in 1..14 (see parse_pili_config).');
end

%% ---- Define sweep values ----
length_values = [1, 2, 3, 4, 5, 6];    % Pilus contour length [um], all pili equal
n_runs        = numel(length_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);

fprintf('========================================\n');
fprintf('  Parameter Sweep: pili length (n_pili = %d, sites = [%s])\n', ...
    n_pili, num2str(site_list));
fprintf('  %d runs scheduled\n', n_runs);
fprintf('========================================\n\n');

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    L = length_values(ii);
    % Same length for each pilus; site per pilus from site_list
    pili_cfg_i = struct('type', 'linear', 'length_um', L, 'site', num2cell(site_list));

    cfg_i             = cfg_base;
    cfg_i.Pili_Matrix = parse_pili_config(pili_cfg_i);
    site_tag = sprintf('%d-', site_list);
    site_tag = site_tag(1:end-1);   % e.g. "1-6-11"
    cfg_i.case_name   = sprintf('sweep-plen-%02d-um-sites-%s', L, site_tag);

    fprintf('\n--- Run %d/%d: length = %d um, sites [%s] ---\n', ...
        ii, n_runs, L, num2str(site_list));

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
        fprintf('  length = %2d um  (sites [%s])  |  Steady V_sink = %.4e um/s\n', ...
            length_values(ii), num2str(site_list), ...
            mean(Vs(max(1, round(0.95 * numel(Vs))):end)));
    else
        fprintf('  length = %2d um  (sites [%s])  |  ABORTED\n', ...
            length_values(ii), num2str(site_list));
    end
end
fprintf('========================================\n');
