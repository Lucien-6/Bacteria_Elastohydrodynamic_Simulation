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
%    Output locations:
%      - Per-run products: <project_root>/output/<case_name>/
%      - Timestamped sweep log: beside this driver under batch_sweep/
%      - cfg_base.compute_flowfield is omitted, so run_simulation() uses its
%        backward-compatible default true and generates per-run flow fields.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-30
%  Updated: 2026-05-06 — Lucien — added geometry_mode='parametric' (required by
%                        save_sim_config since v6/v7 custom-geometry feature),
%                        project-root cd-block, and pili_ext capture from
%                        parse_pili_config for forward compatibility.
%  Updated: 2026-07-12

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
    error('batch_sweep_seg_per_um:BadRoot', 'Cannot resolve project root.');
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

pili_cfg(1).type      = 'linear';
pili_cfg(1).length_um = 3;
pili_cfg(1).site      = 1;
[cfg_base.Pili_Matrix, cfg_base.pili_ext] = parse_pili_config(pili_cfg);
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

%% ---- Summary log setup (persistent, timestamped) ----
[sweep_fid, sweep_log_path] = open_sweep_summary_log(mfilename('fullpath'));
run_status   = repmat("PENDING", n_runs, 1);
run_walltime = nan(n_runs, 1);
run_errmsg   = repmat("", n_runs, 1);
if sweep_fid ~= -1
    fprintf(sweep_fid, 'script    : batch_sweep_seg_per_um\n');
    fprintf(sweep_fid, 'sweep_var : seg_per_um [-]\n');
    fprintf(sweep_fid, 'n_runs    : %d\n', n_runs);
    fprintf(sweep_fid, 'values    : %s\n', mat2str(seg_values));
    fprintf(sweep_fid, 'started   : %s\n\n', char(datetime('now', 'Format', 'uuuu-MM-dd HH:mm:ss')));
end

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    sg = seg_values(ii);
    cfg_i             = cfg_base;
    cfg_i.seg_per_um  = sg;
    cfg_i.case_name   = sprintf('sweep-segperum-%d', sg);

    fprintf('\n--- Run %d/%d: seg_per_um = %g ---\n', ii, n_runs, sg);

    t_run_start = tic;
    try
        [res_i, ~, dir_i] = run_simulation(cfg_i);
        run_walltime(ii)  = toc(t_run_start);
        sweep_results{ii} = res_i;
        sweep_dirs{ii}    = dir_i;
        if isempty(res_i)
            run_status(ii) = "ABORTED";
        else
            run_status(ii) = "OK";
        end
    catch ME
        run_walltime(ii)  = toc(t_run_start);
        run_status(ii)    = "FAILED";
        run_errmsg(ii)    = string(ME.message);
        sweep_results{ii} = [];
        sweep_dirs{ii}    = '';
        fprintf('\n  [FAIL] seg_per_um = %g: %s\n', sg, ME.message);
    end

    if sweep_fid ~= -1
        fprintf(sweep_fid, '[%02d/%02d] seg_per_um = %5g  status = %-8s  wall = %8.2f s\n', ...
            ii, n_runs, sg, run_status(ii), run_walltime(ii));
        if run_status(ii) == "FAILED"
            fprintf(sweep_fid, '         error: %s\n', run_errmsg(ii));
        end
    end
end

t_sweep_total = toc(t_sweep_start);

%% ---- Summary ----
n_ok  = sum(run_status == "OK");
n_err = sum(run_status == "FAILED");
n_abt = sum(run_status == "ABORTED");

fprintf('\n========================================\n');
fprintf('  Sweep Complete  (%d runs, %.1f s total)\n', n_runs, t_sweep_total);
fprintf('  OK: %d  |  FAILED: %d  |  ABORTED: %d\n', n_ok, n_err, n_abt);
fprintf('========================================\n');

for ii = 1:n_runs
    if ~isempty(sweep_results{ii}) && run_status(ii) == "OK"
        Vs = sweep_results{ii}.V_sink;
        fprintf('  seg_per_um = %5g  |  Steady V_sink = %.4e um/s\n', ...
            seg_values(ii), mean(Vs(max(1, round(0.95 * numel(Vs))):end)));
    else
        fprintf('  seg_per_um = %5g  |  %s\n', seg_values(ii), run_status(ii));
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
