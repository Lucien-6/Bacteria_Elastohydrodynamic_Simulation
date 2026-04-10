%% batch_sweep_N_body.m
%  Automated Parameter Sweep: Body Surface Discretisation Count N_body
%
%  Description:
%    Sweeps N_body on a spherical body of radius R = 0.4 um (semi-axes
%    major_axis = minor_axis = 0.4e-6 m) with no pili (empty Pili_Matrix).
%    N_body must be >= 10 (validator).  Results compared via V_sink steady mean.
%
%    Key points:
%      - Sphere geometry: major_axis == minor_axis (validator allows equality)
%      - No pili: Pili_Matrix = zeros(8, 0)
%      - Edit N_body_values below for the sweep grid
%      - auto_confirm = true; preview_initial = false
%      - generate_animation = true; animation_skip = 5
%      - enable_brown = false
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-30
%  Updated: 2026-03-30 — sphere R=0.4 um, no pili
%           2026-04-07 — cfg_base.fixman_midpoint (run_simulation also defaults if omitted)

close all force;
clear; clc;

%% ---- Define base configuration ----
R_sphere = 0.4e-6;   % Sphere radius [m]

cfg_base.case_name  = '';
cfg_base.minor_axis = R_sphere;
cfg_base.major_axis = R_sphere;
cfg_base.N_body     = 200;           % Overridden each run
cfg_base.alpha_0    = 0;

cfg_base.Pili_Matrix  = zeros(8, 0); % No pili
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
%  Integers >= 10
N_body_values = 200:200:1600;
n_runs        = numel(N_body_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);

fprintf('========================================\n');
fprintf('  Parameter Sweep: N_body (sphere R = 0.4 um, no pili)\n');
fprintf('  %d runs scheduled\n', n_runs);
fprintf('========================================\n\n');

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    Nb = N_body_values(ii);
    cfg_i          = cfg_base;
    cfg_i.N_body   = Nb;
    cfg_i.case_name = sprintf('sweep-nbody-%d', Nb);

    fprintf('\n--- Run %d/%d: N_body = %d ---\n', ii, n_runs, Nb);

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
        fprintf('  N_body = %6d  |  Steady V_sink = %.4e um/s\n', ...
            N_body_values(ii), mean(Vs(round(0.9 * numel(Vs)):end)));
    else
        fprintf('  N_body = %6d  |  ABORTED\n', N_body_values(ii));
    end
end
fprintf('========================================\n');
