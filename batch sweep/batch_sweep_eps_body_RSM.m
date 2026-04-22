%% batch_sweep_eps_body_RSM.m
%  Parameter Sweep: RSM body regularisation eps_body (no pili, sphere)
%
%  Description:
%    Sweeps dimensionless RSM regularisation length eps_body on a spherical
%    body (major_axis = minor_axis = 0.4 um) with no pili.  Grid values are
%      eps_body = k * h_mean,   k = 0.1 : 0.1 : 0.5
%    i.e. from h_mean/10 to h_mean/2 in steps of h_mean/10 (5 runs).
%    Brownian motion is disabled.
%
%    Requires project root on the MATLAB path (this script adds it below).
%    Path policy: main entry + run_simulation() use addpath(project root);
%    batch/postprocess keep local addpath for standalone runs (see
%    main_coupled_simulation.m header).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-07

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
    error('batch_sweep_eps_body_RSM:BadRoot', 'Cannot resolve project root.');
end
cd(proj_root);
addpath(proj_root);

%% ---- Geometry and fixed numerics ----
R_sphere = 0.4e-6;   % Sphere radius [m] (semi-axes)

N_body   = 200;      % Surface discretisation (edit to match your study)

%% ---- Base configuration (no pili, no Brownian) ----
cfg_base = struct();
cfg_base.case_name  = '';
cfg_base.minor_axis = R_sphere;
cfg_base.major_axis = R_sphere;
cfg_base.N_body     = N_body;
cfg_base.alpha_0    = 0;

cfg_base.Pili_Matrix  = zeros(8, 0);
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
cfg_base.fixman_midpoint = false;

cfg_base.preview_initial    = false;
cfg_base.generate_animation = false;
cfg_base.animation_fps      = 30;
cfg_base.animation_skip     = 5;

cfg_base.auto_confirm = true;

cfg_base.flowfield_grid_spacing = 0.05;
cfg_base.flowfield_grid_range   = 10.0;
cfg_base.flowfield_quiver_skip  = 5;
cfg_base.flowfield_chunk_size   = 2000;

%% ---- h_mean from mesh (same N_body and sphere as sweep) ----
params_ref = init_coupled_params( ...
    'N_body',     N_body, ...
    'minor_axis', R_sphere, ...
    'major_axis', R_sphere, ...
    'Pili_Matrix', zeros(8, 0));

h_mean = params_ref.body_h_mean;

k_vals = 0.1 : 0.1 : 0.5;
eps_body_values = k_vals * h_mean;
n_runs = numel(eps_body_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);
steady_V      = nan(n_runs, 1);
run_status    = repmat("PENDING", n_runs, 1);
run_walltime  = nan(n_runs, 1);
run_errmsg    = repmat("", n_runs, 1);

fprintf('========================================\n');
fprintf('  RSM eps_body sweep (sphere R = %.2f um, no pili, no Brownian)\n', R_sphere * 1e6);
fprintf('  N_body = %d, h_mean = %.6f L0\n', N_body, h_mean);
fprintf('  eps_body = k * h_mean, k = [%s]\n', num2str(k_vals, '%.1f '));
fprintf('  %d runs\n', n_runs);
fprintf('========================================\n\n');

%% ---- Summary log setup (persistent, timestamped) ----
[sweep_fid, sweep_log_path] = open_sweep_summary_log(mfilename('fullpath'));
if sweep_fid ~= -1
    fprintf(sweep_fid, 'script       : batch_sweep_eps_body_RSM\n');
    fprintf(sweep_fid, 'sweep_var    : eps_body = k * h_mean\n');
    fprintf(sweep_fid, 'N_body       : %d\n', N_body);
    fprintf(sweep_fid, 'h_mean (L0)  : %.6g\n', h_mean);
    fprintf(sweep_fid, 'k_values     : %s\n', mat2str(k_vals, 4));
    fprintf(sweep_fid, 'eps_body vals: %s\n', mat2str(eps_body_values, 4));
    fprintf(sweep_fid, 'n_runs       : %d\n', n_runs);
    fprintf(sweep_fid, 'started      : %s\n\n', char(datetime('now', 'Format', 'uuuu-MM-dd HH:mm:ss')));
end

t_sweep_start = tic;

for ii = 1:n_runs
    cfg_i = cfg_base;
    cfg_i.eps_body = eps_body_values(ii);
    cfg_i.case_name = sprintf('sweep-epsRSM-k%02d', round(100 * k_vals(ii)));

    fprintf('\n--- Run %d/%d: eps_body/h_mean = %.2f (eps_body = %.6f L0) ---\n', ...
        ii, n_runs, k_vals(ii), eps_body_values(ii));

    t_run_start = tic;
    try
        [res_i, ~, dir_i] = run_simulation(cfg_i);
        run_walltime(ii)  = toc(t_run_start);
        sweep_results{ii} = res_i;
        sweep_dirs{ii}    = dir_i;

        if ~isempty(res_i) && isfield(res_i, 'V_sink')
            Vs = res_i.V_sink(:);
            nt = numel(Vs);
            steady_V(ii) = mean(Vs(max(1, round(0.95 * nt)):end));
            run_status(ii) = "OK";
        else
            run_status(ii) = "ABORTED";
        end
    catch ME
        run_walltime(ii)  = toc(t_run_start);
        run_status(ii)    = "FAILED";
        run_errmsg(ii)    = string(ME.message);
        sweep_results{ii} = [];
        sweep_dirs{ii}    = '';
        fprintf('\n  [FAIL] k = %.2f: %s\n', k_vals(ii), ME.message);
    end

    if sweep_fid ~= -1
        fprintf(sweep_fid, '[%02d/%02d] k = %.2f  eps = %.6g  status = %-8s  V_sink = %12.4g  wall = %8.2f s\n', ...
            ii, n_runs, k_vals(ii), eps_body_values(ii), run_status(ii), steady_V(ii), run_walltime(ii));
        if run_status(ii) == "FAILED"
            fprintf(sweep_fid, '         error: %s\n', run_errmsg(ii));
        end
    end
end

t_sweep_total = toc(t_sweep_start);

%% ---- Summary table ----
n_ok  = sum(run_status == "OK");
n_err = sum(run_status == "FAILED");
n_abt = sum(run_status == "ABORTED");

fprintf('\n========================================\n');
fprintf('  Sweep complete  (%d runs, %.1f s total)\n', n_runs, t_sweep_total);
fprintf('  OK: %d  |  FAILED: %d  |  ABORTED: %d\n', n_ok, n_err, n_abt);
fprintf('========================================\n');

sum_tbl = table( ...
    k_vals(:), ...
    eps_body_values(:), ...
    steady_V(:), ...
    'VariableNames', {'k_eps_over_h_mean', 'eps_body_dimless', 'steady_V_sink_um_s'});

disp(sum_tbl);

out_csv = fullfile(proj_root, 'batch sweep', 'batch_sweep_eps_body_RSM_summary.csv');
try
    writetable(sum_tbl, out_csv);
    fprintf('\nSummary written: %s\n', out_csv);
catch ME
    warning('batch_sweep_eps_body_RSM:WriteCSV', '%s', ME.message);
end

out_mat = fullfile(proj_root, 'batch sweep', 'batch_sweep_eps_body_RSM_results.mat');
try
    save(out_mat, 'sum_tbl', 'sweep_results', 'sweep_dirs', 'h_mean', 'N_body', ...
        'R_sphere', 'k_vals', 'eps_body_values', ...
        'run_status', 'run_walltime', 'run_errmsg', '-v7.3');
    fprintf('Full results saved: %s\n', out_mat);
catch ME
    warning('batch_sweep_eps_body_RSM:WriteMAT', '%s', ME.message);
end

if sweep_fid ~= -1
    fprintf(sweep_fid, '\nelapsed_total_s: %.2f\n', t_sweep_total);
    fprintf(sweep_fid, 'success_count  : %d / %d  (FAILED: %d, ABORTED: %d)\n', ...
        n_ok, n_runs, n_err, n_abt);
    fprintf(sweep_fid, 'summary_csv    : %s\n', out_csv);
    fprintf(sweep_fid, 'results_mat    : %s\n', out_mat);
    fclose(sweep_fid);
    fprintf('\nSweep summary log: %s\n', sweep_log_path);
end

fprintf('========================================\n');

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
