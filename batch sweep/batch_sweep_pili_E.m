%% batch_sweep_pili_E.m
%  Automated Parameter Sweep: Pilus Young's Modulus
%
%  Description:
%    Sweeps the bending Young's modulus E of the pili from 1e5 Pa to 1e9 Pa
%    (tenfold per step: five values) while keeping all geometry and fluid parameters
%    fixed.  Higher E yields a stiffer, less deformable pilus.  Results
%    are compared via the steady-state sinking speed V_sink.
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
cfg_base.N_body     = 200;
cfg_base.alpha_0    = 0;

pili_cfg(1).type      = 'linear';
pili_cfg(1).length_um = 3;
pili_cfg(1).site      = 1;
cfg_base.Pili_Matrix  = parse_pili_config(pili_cfg);
cfg_base.seg_per_um   = 10;

cfg_base.E       = 3.279e8;     % Overridden in the loop
cfg_base.r_pili  = 3e-9;
cfg_base.sigma_p = 0.2;

cfg_base.mu        = 1e-3;
cfg_base.rho_bac   = 1180.0;
cfg_base.rho_fluid = 998.232;

cfg_base.T_obs  = 360.0;
cfg_base.AbsTol = 1e-4;
cfg_base.RelTol = 1e-4;
cfg_base.delta  = pi / 20;
cfg_base.tlim   = 3600;
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
% Young's modulus [Pa]: 1e5 to 1e9, multiplicative step 10
E_values = 10.^(5:9);
n_runs   = numel(E_values);

sweep_results = cell(n_runs, 1);
sweep_dirs    = cell(n_runs, 1);

fprintf('========================================\n');
fprintf('  Parameter Sweep: pili Young''s modulus E\n');
fprintf('  %d runs scheduled\n', n_runs);
fprintf('========================================\n\n');

%% ---- Execute sweep ----
t_sweep_start = tic;

for ii = 1:n_runs
    cfg_i           = cfg_base;
    cfg_i.E         = E_values(ii);
    % Compact case name: remove '+' and replace '.' with 'p' (e.g. 1.77e+07 -> 1p77e07)
    raw_str         = sprintf('%.2e', E_values(ii));
    raw_str         = strrep(raw_str, '+', '');
    raw_str         = strrep(raw_str, '.', 'p');
    cfg_i.case_name = sprintf('sweep-E-%s', raw_str);

    fprintf('\n--- Run %d/%d: E = %.3e Pa ---\n', ii, n_runs, E_values(ii));

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
        fprintf('  E = %.3e Pa  |  Steady V_sink = %.4e um/s\n', ...
            E_values(ii), mean(Vs(max(1, round(0.95 * numel(Vs))):end)));
    else
        fprintf('  E = %.3e Pa  |  ABORTED\n', E_values(ii));
    end
end
fprintf('========================================\n');
