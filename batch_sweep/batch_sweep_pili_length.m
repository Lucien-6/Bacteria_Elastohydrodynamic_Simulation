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
%    Output locations:
%      - Per-run products: <project_root>/output/<case_name>/
%      - Timestamped sweep log: beside this driver under batch_sweep/
%      - cfg_base.compute_flowfield is omitted, so run_simulation() uses its
%        backward-compatible default true and generates per-run flow fields.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-18
%  Updated: 2026-05-06 — Lucien — added geometry_mode='parametric' (required by
%                        save_sim_config since v6/v7 custom-geometry feature),
%                        project-root cd-block, and pili_ext capture from
%                        parse_pili_config for forward compatibility.
%  Updated: 2026-07-13

close all force;
clear; clc;

%% ---- Project root (robust to current directory) ----
this_file = mfilename('fullpath');
if isempty(this_file)
    proj_root = pwd;
else
    proj_root = fileparts(fileparts(this_file));
end
if ~isfolder(proj_root)
    error('batch_sweep_pili_length:BadRoot', 'Cannot resolve project root.');
end
cd(proj_root);
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));

%% ---- Define base configuration ----
cfg_base.case_name     = '';
cfg_base.geometry_mode = 'parametric'; % required by save_sim_config (v6+)
cfg_base.minor_axis    = 0.4e-6;
cfg_base.major_axis    = 1.25e-6;
cfg_base.N_body        = 200;
cfg_base.alpha_0       = 0;

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
cfg_base.tlim   = 3600;
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

%% ---- Summary log setup (persistent, timestamped) ----
[sweep_fid, sweep_log_path] = open_sweep_summary_log(mfilename('fullpath'));
run_status   = repmat("PENDING", n_runs, 1);
run_walltime = nan(n_runs, 1);
run_errmsg   = repmat("", n_runs, 1);
if sweep_fid ~= -1
    fprintf(sweep_fid, 'script    : batch_sweep_pili_length\n');
    fprintf(sweep_fid, 'sweep_var : pilus length [um] (shared across %d pili)\n', n_pili);
    fprintf(sweep_fid, 'sites     : [%s]\n', num2str(site_list));
    fprintf(sweep_fid, 'n_runs    : %d\n', n_runs);
    fprintf(sweep_fid, 'values_um : %s\n', mat2str(length_values));
    fprintf(sweep_fid, 'started   : %s\n\n', char(datetime('now', 'Format', 'uuuu-MM-dd HH:mm:ss')));
end

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    L = length_values(ii);
    % Same length for each pilus; site per pilus from site_list
    pili_cfg_i = struct('type', 'linear', 'length_um', L, 'site', num2cell(site_list));

    cfg_i             = cfg_base;
    [cfg_i.Pili_Matrix, cfg_i.pili_ext] = parse_pili_config(pili_cfg_i);
    site_tag = sprintf('%d-', site_list);
    site_tag = site_tag(1:end-1);   % e.g. "1-6-11"
    cfg_i.case_name   = sprintf('sweep-plen-%02d-um-sites-%s', L, site_tag);

    fprintf('\n--- Run %d/%d: length = %d um, sites [%s] ---\n', ...
        ii, n_runs, L, num2str(site_list));

    t_run_start = tic;
    try
        [res_i, ~, dir_i] = run_simulation(cfg_i);
        run_walltime(ii)  = toc(t_run_start);
        sweep_results{ii} = res_i;
        sweep_dirs{ii}    = dir_i;
        run_status(ii) = simulation_result_status(res_i);
    catch ME
        run_walltime(ii)  = toc(t_run_start);
        run_status(ii)    = "FAILED";
        run_errmsg(ii)    = string(ME.message);
        sweep_results{ii} = [];
        sweep_dirs{ii}    = '';
        fprintf('\n  [FAIL] length = %d um: %s\n', L, ME.message);
        if ~isempty(ME.stack)
            fprintf('         Location: %s (line %d)\n', ...
                ME.stack(1).name, ME.stack(1).line);
        end
    end

    if sweep_fid ~= -1
        fprintf(sweep_fid, '[%02d/%02d] length = %2d um  status = %-8s  wall = %8.2f s\n', ...
            ii, n_runs, L, run_status(ii), run_walltime(ii));
        if run_status(ii) == "FAILED"
            fprintf(sweep_fid, '         Identifier : %s\n', ME.identifier);
            fprintf(sweep_fid, '         Message    : %s\n', run_errmsg(ii));
            if ~isempty(ME.stack)
                fprintf(sweep_fid, '         Location   : %s (line %d)\n', ...
                    ME.stack(1).name, ME.stack(1).line);
            end
        end
    end
end

t_sweep_total = toc(t_sweep_start);

%% ---- Summary ----
n_ok  = sum(run_status == "OK");
n_err = sum(run_status == "FAILED");
n_abt = sum(run_status == "ABORTED" | run_status == "TIMEOUT");

fprintf('\n========================================\n');
fprintf('  Sweep Complete  (%d runs, %.1f s total)\n', n_runs, t_sweep_total);
fprintf('  OK: %d  |  FAILED: %d  |  ABORTED: %d\n', n_ok, n_err, n_abt);
fprintf('========================================\n');

for ii = 1:n_runs
    if ~isempty(sweep_results{ii}) && run_status(ii) == "OK"
        Vs = sweep_results{ii}.V_sink;
        fprintf('  length = %2d um  (sites [%s])  |  Steady V_sink = %.4e um/s\n', ...
            length_values(ii), num2str(site_list), ...
            mean(Vs(max(1, round(0.95 * numel(Vs))):end)));
    else
        fprintf('  length = %2d um  (sites [%s])  |  %s\n', ...
            length_values(ii), num2str(site_list), run_status(ii));
    end
end
fprintf('========================================\n');

if sweep_fid ~= -1
    fprintf(sweep_fid, '\nelapsed_total_s: %.2f\n', t_sweep_total);
    fprintf(sweep_fid, 'success_count  : %d / %d  (FAILED: %d, ABORTED: %d)\n', ...
        n_ok, n_runs, n_err, n_abt);
    fclose(sweep_fid);
    fprintf('\nSweep summary log: %s\n', sweep_log_path);
end

%% ---- Local helpers ----
function [fid, log_path] = open_sweep_summary_log(script_fullpath)
    ts_str  = char(datetime('now', 'Format', 'uuuuMMdd''_''HHmmss'));
    [d, nm] = fileparts(script_fullpath);
    if isempty(d)
        d = fullfile(pwd, 'batch sweep');
    end
    if ~isfolder(d), mkdir(d); end
    log_path = fullfile(d, sprintf('%s_summary_%s.log', nm, ts_str));
    fid = fopen(log_path, 'w', 'n', 'UTF-8');
    if fid == -1
        warning('batch_sweep:LogOpenFailed', ...
            'Cannot open sweep summary log: %s', log_path);
        log_path = '';
    end
end
