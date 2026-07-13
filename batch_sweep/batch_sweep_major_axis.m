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
%    Output locations:
%      - Per-run products: <project_root>/output/<case_name>/
%      - Timestamped sweep log: beside this driver under batch_sweep/
%      - cfg_base.compute_flowfield is omitted, so run_simulation() uses its
%        backward-compatible default true and generates per-run flow fields.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-18
%  Updated: 2026-05-06 — Lucien — added geometry_mode='parametric' (required by
%                        save_sim_config since v6/v7 custom-geometry feature)
%                        and project-root cd-block for robust standalone runs.
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
    error('batch_sweep_major_axis:BadRoot', 'Cannot resolve project root.');
end
cd(proj_root);
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));

%% ---- Define base configuration ----
cfg_base.case_name     = '';           % Will be auto-set per run below
cfg_base.geometry_mode = 'parametric'; % required by save_sim_config (v6+)
cfg_base.minor_axis    = 0.4e-6;
cfg_base.major_axis    = 1.25e-6;      % Overridden in the loop
cfg_base.N_body        = 200;
cfg_base.alpha_0       = 0;

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

%% ---- Summary log setup (persistent, timestamped) ----
[sweep_fid, sweep_log_path] = open_sweep_summary_log(mfilename('fullpath'));
run_status   = repmat("PENDING", n_runs, 1);
run_walltime = nan(n_runs, 1);
run_errmsg   = repmat("", n_runs, 1);
if sweep_fid ~= -1
    fprintf(sweep_fid, 'script    : batch_sweep_major_axis\n');
    fprintf(sweep_fid, 'sweep_var : major_axis [um]\n');
    fprintf(sweep_fid, 'n_runs    : %d\n', n_runs);
    fprintf(sweep_fid, 'values_um : %s\n', mat2str(major_values * 1e6, 4));
    fprintf(sweep_fid, 'started   : %s\n\n', char(datetime('now', 'Format', 'uuuu-MM-dd HH:mm:ss')));
end

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    cfg_i = cfg_base;
    cfg_i.major_axis = major_values(ii);
    cfg_i.case_name  = sprintf('sweep-major-%04.0f-nm', major_values(ii) * 1e9);

    fprintf('\n--- Run %d/%d: major_axis = %.2f um ---\n', ...
        ii, n_runs, major_values(ii) * 1e6);

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
        fprintf('\n  [FAIL] major_axis = %.2f um: %s\n', ...
            major_values(ii) * 1e6, ME.message);
    end

    if sweep_fid ~= -1
        fprintf(sweep_fid, '[%02d/%02d] major_axis = %5.2f um  status = %-8s  wall = %8.2f s\n', ...
            ii, n_runs, major_values(ii) * 1e6, run_status(ii), run_walltime(ii));
        if run_status(ii) == "FAILED"
            fprintf(sweep_fid, '         error: %s\n', run_errmsg(ii));
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
        fprintf('  major_axis = %5.2f um  |  Steady V_sink = %.4e um/s\n', ...
            major_values(ii) * 1e6, mean(Vs(max(1, round(0.95 * numel(Vs))):end)));
    else
        fprintf('  major_axis = %5.2f um  |  %s\n', major_values(ii) * 1e6, run_status(ii));
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
