%% batch_sweep_N_body_eps_RSM.m
%  Nested sweep: N_body x (eps_body = k * h_mean) for RSM calibration
%
%  Description:
%    For each body surface discretisation count N_body, computes h_mean from
%    init_coupled_params (same geometry as batch_sweep_eps_body_RSM.m: sphere
%    R = 0.4 um, no pili), then sweeps dimensionless regularisation
%      eps_body = k * h_mean
%    over a grid of k values.  Steady vertical velocity V_sink [um/s] is
%    recorded for each run.
%
%    For each N_body, the script finds k_star that minimises the absolute
%    difference between steady V_sink and the Stokes settling speed for the
%    same sphere (same rho_bac, rho_fluid, mu as cfg).  This tests whether an
%    approximately constant optimal ratio k = eps_body / h_mean exists across
%    mesh refinements.
%
%    Outputs:
%      - batch_sweep_N_body_eps_RSM_grid.csv   — full (N_body, k) table
%      - batch_sweep_N_body_eps_RSM_summary.csv — one row per N_body: h_mean, k_star, error
%      - batch_sweep_N_body_eps_RSM_results.mat — matrices + cell results
%
%    Notes:
%      - Set cfg_base.compute_flowfield = false and generate_animation = false
%        to keep run time manageable (edit if needed).
%      - Reduce cfg_base.T_obs for a quick pilot; restore for production.
%
%  Requires project root on the MATLAB path (or run from project root).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-13

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
    error('batch_sweep_N_body_eps_RSM:BadRoot', 'Cannot resolve project root.');
end
cd(proj_root);
addpath(proj_root);

%% ---- Geometry (sphere) and physical constants [SI] ----
R_sphere = 0.4e-6;   % Sphere radius [m] (semi-axes)

mu        = 1e-3;
rho_bac   = 1180.0;
rho_fluid = 998.232;
g         = 9.81;
Delta_rho = rho_bac - rho_fluid;

% Match capsule volume formula for sphere (major_axis = minor_axis = R)
b = R_sphere;
a = R_sphere;
V_capsule = pi * b^2 * (4/3 * b + 2 * (a - b));
F_sink    = Delta_rho * V_capsule * g;

% Stokes drag: F = 6*pi*mu*R*|V|; balance F_sink = F_drag (settling: V_z < 0)
V_stokes_m_s   = -F_sink / (6 * pi * mu * R_sphere);
V_stokes_um_s  = V_stokes_m_s * 1e6;

%% ---- Sweep grids (edit) ----
% N_body: integers >= 10 (validator)
N_body_values = 200:200:1600;
% k = eps_body / h_mean (same convention as batch_sweep_eps_body_RSM.m)
k_vals = 0.1 : 0.1 : 0.5;

n_N = numel(N_body_values);
n_k = numel(k_vals);
n_total = n_N * n_k;

%% ---- Base configuration (no pili, no Brownian) ----
cfg_base = struct();
cfg_base.case_name  = '';
cfg_base.minor_axis = R_sphere;
cfg_base.major_axis = R_sphere;
cfg_base.alpha_0    = 0;

cfg_base.Pili_Matrix  = zeros(8, 0);
cfg_base.seg_per_um   = 10;

cfg_base.E       = 3.279e8;
cfg_base.r_pili  = 3e-9;
cfg_base.sigma_p = 0.2;

cfg_base.mu        = mu;
cfg_base.rho_bac   = rho_bac;
cfg_base.rho_fluid = rho_fluid;

cfg_base.T_obs  = 200.0;
cfg_base.AbsTol = 1e-4;
cfg_base.RelTol = 1e-4;
cfg_base.delta  = pi / 20;
cfg_base.tlim   = 360;
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
cfg_base.compute_flowfield      = false;

%% ---- Preallocate ----
steady_V_grid = nan(n_N, n_k);
h_mean_list   = nan(n_N, 1);
sweep_cell    = cell(n_N, n_k);

fprintf('========================================\n');
fprintf('  Nested sweep: N_body x k (eps_body = k * h_mean)\n');
fprintf('  Sphere R = %.2f um | no pili | no Brownian | no flowfield/anim\n', R_sphere * 1e6);
fprintf('  Stokes reference V_sink = %.6e um/s (same F_sink as simulation)\n', V_stokes_um_s);
fprintf('  N_body in [%d, %d] (%d values); k in [%.2f, %.2f] (%d values); total runs = %d\n', ...
    min(N_body_values), max(N_body_values), n_N, k_vals(1), k_vals(end), n_k, n_total);
fprintf('========================================\n\n');

t_sweep_start = tic;
run_id = 0;
rows_grid = nan(n_total, 6);

for iN = 1:n_N
    Nb = N_body_values(iN);

    params_ref = init_coupled_params( ...
        'N_body',      Nb, ...
        'minor_axis',  R_sphere, ...
        'major_axis',  R_sphere, ...
        'Pili_Matrix', zeros(8, 0));

    h_m = params_ref.body_h_mean;
    h_mean_list(iN) = h_m;

    fprintf('--- N_body = %d  |  h_mean = %.6f L0  |  default eps/h_mean would be %.4f ---\n', ...
        Nb, h_m, params_ref.eps_body / h_m);

    for ik = 1:n_k
        run_id = run_id + 1;
        k = k_vals(ik);
        eps_body_i = k * h_m;

        cfg_i = cfg_base;
        cfg_i.N_body     = Nb;
        cfg_i.eps_body   = eps_body_i;
        cfg_i.case_name  = sprintf('sweep-nbody-%d-k%03d', Nb, round(100 * k));

        fprintf('  [%d/%d] N_body=%d  k=%.2f  eps_body=%.6f L0\n', ...
            run_id, n_total, Nb, k, eps_body_i);

        [res_i, ~, ~] = run_simulation(cfg_i);
        sweep_cell{iN, ik} = res_i;

        if ~isempty(res_i) && isfield(res_i, 'V_sink')
            Vs = res_i.V_sink(:);
            nt = numel(Vs);
            steady_V_grid(iN, ik) = mean(Vs(round(0.9 * nt):end));
        end

        rows_grid(run_id, :) = [Nb, h_m, k, eps_body_i, steady_V_grid(iN, ik), ...
            abs(steady_V_grid(iN, ik) - V_stokes_um_s)];
    end
    fprintf('\n');
end

t_sweep_total = toc(t_sweep_start);

%% ---- Per-N_body: k_star minimising |V_steady - V_Stokes| ----
err_grid = abs(steady_V_grid - V_stokes_um_s);
k_star   = nan(n_N, 1);
min_err  = nan(n_N, 1);
V_at_star = nan(n_N, 1);

for iN = 1:n_N
    row = err_grid(iN, :);
    if all(~isfinite(row))
        continue;
    end
    [min_err(iN), idx_min] = min(row, [], 'omitnan');
    if isempty(idx_min) || isnan(min_err(iN))
        continue;
    end
    k_star(iN)    = k_vals(idx_min);
    V_at_star(iN) = steady_V_grid(iN, idx_min);
end

%% ---- Tables ----
grid_tbl = array2table(rows_grid, 'VariableNames', { ...
    'N_body', 'h_mean_dimless', 'k_eps_over_h_mean', 'eps_body_dimless', ...
    'steady_V_sink_um_s', 'abs_err_vs_Stokes_um_s'});

summary_tbl = table( ...
    N_body_values(:), ...
    h_mean_list(:), ...
    k_star(:), ...
    min_err(:), ...
    V_at_star(:), ...
    repmat(V_stokes_um_s, n_N, 1), ...
    'VariableNames', {'N_body', 'h_mean_dimless', 'k_star', ...
    'min_abs_err_vs_Stokes_um_s', 'steady_V_at_k_star_um_s', 'V_Stokes_ref_um_s'});

fprintf('========================================\n');
fprintf('  Sweep complete  (%.1f s wall time)\n', t_sweep_total);
fprintf('========================================\n');
fprintf('  k_star vs Stokes (min |V - V_Stokes| per N_body):\n');
disp(summary_tbl);

out_csv_grid = fullfile(proj_root, 'batch sweep', 'batch_sweep_N_body_eps_RSM_grid.csv');
out_csv_sum  = fullfile(proj_root, 'batch sweep', 'batch_sweep_N_body_eps_RSM_summary.csv');

try
    writetable(grid_tbl, out_csv_grid);
    fprintf('\nGrid table written: %s\n', out_csv_grid);
catch ME
    warning('batch_sweep_N_body_eps_RSM:WriteCSVGrid', '%s', ME.message);
end

try
    writetable(summary_tbl, out_csv_sum);
    fprintf('Summary table written: %s\n', out_csv_sum);
catch ME
    warning('batch_sweep_N_body_eps_RSM:WriteCSVSummary', '%s', ME.message);
end

out_mat = fullfile(proj_root, 'batch sweep', 'batch_sweep_N_body_eps_RSM_results.mat');
try
    save(out_mat, ...
        'grid_tbl', 'summary_tbl', 'steady_V_grid', 'err_grid', 'sweep_cell', ...
        'N_body_values', 'k_vals', 'h_mean_list', 'k_star', 'V_stokes_um_s', ...
        'R_sphere', 'F_sink', 'V_stokes_m_s', 't_sweep_total', '-v7.3');
    fprintf('Full results saved: %s\n', out_mat);
catch ME
    warning('batch_sweep_N_body_eps_RSM:WriteMAT', '%s', ME.message);
end

fprintf('========================================\n');
