%% batch_repeat_brownian.m
%  Independent repeated simulations: same configuration, Brownian motion ON
%
%  Description:
%    Runs the same physical setup multiple times with independent RNG streams
%    for Brownian kicks.  Each repeat is saved under:
%        output/<case_name_base>/001/
%        output/<case_name_base>/002/
%        ...
%    Subfolder names are zero-padded three digits (001, 002, ..., 999+).
%
%    Requirements:
%      - Brownian motion must be enabled (enable_brown = true).  The script
%        errors at startup if this is not satisfied.
%      - This script forces preview_initial = false (no interactive initial
%        pose preview) and auto_confirm = true (no overwrite dialogs).
%
%    Usage:
%      1. Set case_name_base, n_repeats, and either config_json OR the inline
%         cfg_base block below.
%      2. Run this script from the Coupled_Body_Pili project directory
%         (same as other batch scripts).
%      3. Optional: rng_seed_base controls reproducible seeds:
%            cfg.brown_rng_seed = rng_seed_base + k for repeat k (run_simulation sets rng).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-03

close all force;
clear; clc;

% Ensure working directory is project root (contains run_simulation.m)
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
if isfile(fullfile(projectRoot, 'run_simulation.m'))
    cd(projectRoot);
elseif ~isfile('run_simulation.m')
    error('batch_repeat_brownian:WrongDirectory', ...
        ['Could not find run_simulation.m. Place this script in batch sweep/ ', ...
         'under Coupled_Body_Pili, or cd to Coupled_Body_Pili before running.']);
end

%% ========================================================================
%  USER SETTINGS
%  ========================================================================

% --- Root folder name under output/ (required).  Repeats go to:
%     output/<case_name_base>/001, 002, ...
case_name_base = 'Sphere-Brown';

% --- Number of independent repeats (>= 1)
n_repeats = 50;

% --- RNG: seed for repeat k is rng_seed_base + k (reproducible batch runs)
rng_seed_base = 202604091;

% --- Configuration source (choose one):
%     A)  Set config_json to a path to load via load_sim_config (recommended)
%     B)  Leave config_json empty ('') to use the inline cfg_base block below.

config_json = '';   % e.g. fullfile('configs', 'OnePili-P1-L3.json')  (after cd to project root)

%% ========================================================================
%  INLINE CONFIG (used only when config_json is empty)
%  ========================================================================

cfg_base = struct();
cfg_base.case_name  = '';   % Overwritten per repeat; ignored here
cfg_base.minor_axis = 0.4e-6;
cfg_base.major_axis = 0.4e-6;
cfg_base.N_body     = 250;
cfg_base.alpha_0    = pi / 2;

pili_cfg = struct([]);
cfg_base.Pili_Matrix = parse_pili_config(pili_cfg);

cfg_base.seg_per_um = 10;
cfg_base.E          = 3.279e8;
cfg_base.r_pili     = 3e-9;
cfg_base.sigma_p    = 0.2;

cfg_base.mu        = 1e-3;
cfg_base.rho_bac   = 1180.0;
cfg_base.rho_fluid = 998.232;

cfg_base.T_obs  = 200.0;
cfg_base.AbsTol = 1e-4;
cfg_base.RelTol = 1e-4;
cfg_base.delta  = pi / 20;
cfg_base.tlim   = 120;
cfg_base.dt     = 0.1;

% Brownian MUST be on (also checked after load)
cfg_base.enable_brown = true;
cfg_base.Temper       = 293.15;
cfg_base.dt_B         = cfg_base.dt / 2;
cfg_base.fixman_midpoint = false;

cfg_base.preview_initial    = false;
cfg_base.generate_animation = true;
cfg_base.animation_fps      = 30;
cfg_base.animation_skip     = 5;

cfg_base.auto_confirm = true;

cfg_base.flowfield_grid_spacing = 0.05;
cfg_base.flowfield_grid_range   = 10.0;
cfg_base.flowfield_quiver_skip  = 5;
cfg_base.flowfield_chunk_size   = 2000;

%% ---- Load configuration ----
if isempty(strtrim(case_name_base))
    error('batch_repeat_brownian:EmptyCaseName', ...
        'case_name_base must be a non-empty string (folder name under output/).');
end

if n_repeats < 1 || floor(n_repeats) ~= n_repeats
    error('batch_repeat_brownian:InvalidNRepeats', ...
        'n_repeats must be a positive integer (>= 1).');
end

if ~isempty(strtrim(config_json))
    if ~isfile(config_json)
        error('batch_repeat_brownian:ConfigNotFound', ...
            'Config file not found:\n  %s', config_json);
    end
    cfg = load_sim_config(config_json);
else
    cfg = cfg_base;
end

% --- Require Brownian ON ---
if ~isfield(cfg, 'enable_brown') || ~logical(cfg.enable_brown)
    error('batch_repeat_brownian:BrownianRequired', ...
        ['This script is only for Brownian simulations.  Set enable_brown = true ', ...
         '(inline) or brownian.enable = true in JSON.']);
end

% --- Batch defaults: no preview, headless ---
cfg.preview_initial = false;
cfg.auto_confirm    = true;

%% ---- Summary log under output/<case_name_base>/ ----
out_root = fullfile(pwd, 'output', case_name_base);
if ~isfolder(out_root)
    mkdir(out_root);
end

log_ts   = char(datetime('now', 'Format', 'uuuuMMdd''_''HHmmss'));
sum_path = fullfile(out_root, sprintf('batch_repeat_summary_%s.txt', log_ts));
fid_sum  = fopen(sum_path, 'w', 'n', 'UTF-8');
if fid_sum == -1
    sum_path = '';
end

%% ---- Header ----
sep = repmat('=', 1, 60);
fprintf('\n%s\n', sep);
fprintf('  Batch repeat (Brownian): %d independent runs\n', n_repeats);
fprintf('  Output root: %s\n', out_root);
fprintf('  Brownian seeds: cfg.brown_rng_seed = rng_seed_base + k, k = 1..%d  (base = %d)\n', ...
    n_repeats, rng_seed_base);
fprintf('%s\n\n', sep);

if fid_sum ~= -1
    fprintf(fid_sum, 'batch_repeat_brownian\n');
    fprintf(fid_sum, 'case_name_base: %s\n', case_name_base);
    fprintf(fid_sum, 'n_repeats: %d\n', n_repeats);
    fprintf(fid_sum, 'rng_seed_base: %d\n', rng_seed_base);
    fprintf(fid_sum, 'started: %s\n\n', log_ts);
end

run_status  = repmat("PENDING", n_repeats, 1);
run_outdir  = repmat("", n_repeats, 1);
run_errmsg  = repmat("", n_repeats, 1);
run_wall    = nan(n_repeats, 1);

t_batch = tic;

for k = 1:n_repeats
    sub = sprintf('%03d', k);
    fprintf('%s\n  Repeat %d/%d  ->  output/%s/%s\n%s\n', ...
        repmat('-', 1, 56), k, n_repeats, case_name_base, sub, repmat('-', 1, 56));

    cfg_k = cfg;
    cfg_k.case_name = fullfile(case_name_base, sub);

    % Independent, reproducible twister seed for this repeat (see run_simulation.m)
    cfg_k.brown_rng_seed = rng_seed_base + k;

    t_run = tic;
    try
        [res_k, ~, dir_k] = run_simulation(cfg_k);
        run_wall(k) = toc(t_run);

        if isempty(res_k)
            run_status(k) = "FAILED";
            run_errmsg(k) = "run_simulation returned empty results";
            run_outdir(k) = string(dir_k);
            fprintf('  [WARN] Empty results.\n');
        else
            run_status(k) = "OK";
            run_outdir(k) = string(dir_k);
            fprintf('  [OK]   Wall time %.1f s  |  %s\n', run_wall(k), dir_k);
        end
    catch ME
        run_wall(k)     = toc(t_run);
        run_status(k)   = "FAILED";
        run_errmsg(k)   = string(ME.message);
        run_outdir(k)   = "";
        fprintf('  [ERR]  %s\n', ME.message);
    end

    if fid_sum ~= -1
        fprintf(fid_sum, '%s  seed=%d  status=%s  wall=%.2f s\n', ...
            sub, rng_seed_base + k, run_status(k), run_wall(k));
        if run_status(k) == "FAILED"
            fprintf(fid_sum, '  error: %s\n', run_errmsg(k));
        elseif strlength(run_outdir(k)) > 0
            fprintf(fid_sum, '  out: %s\n', run_outdir(k));
        end
    end
end

elapsed = toc(t_batch);

fprintf('\n%s\n', sep);
fprintf('  Batch complete: %.1f s total\n', elapsed);
n_ok = sum(run_status == "OK");
fprintf('  Success: %d / %d\n', n_ok, n_repeats);
fprintf('%s\n', sep);

if fid_sum ~= -1
    fprintf(fid_sum, '\nelapsed_total_s: %.2f\n', elapsed);
    fprintf(fid_sum, 'success_count: %d / %d\n', n_ok, n_repeats);
    fclose(fid_sum);
    fprintf('\n  Summary written: %s\n', sum_path);
end
