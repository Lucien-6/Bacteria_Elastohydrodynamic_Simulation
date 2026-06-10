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
%      3. Brownian seed mode (see brown_rng_seed_mode below):
%            'fixed'     — cfg.brown_rng_seed = rng_seed_base + k (reproducible batch).
%            'timestamp' — SERIAL branch: omit cfg.brown_rng_seed; run_simulation uses
%                          MMddHHmmss (brown_rng_seed_from_clock) at each run start.
%                          PARALLEL branch: seeds are pre-assigned as clock_base + k so
%                          workers cannot collide on the same MMddHHmmss second.
%            'default'   — alias for 'timestamp' (pre-rename value, retained for
%                          backward compatibility; normalised silently at runtime).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-03
%  Updated: 2026-05-06 — Lucien — optional parfor batch path, pool sizing,
%                        PCT fallback; cfg.geometry_mode='parametric' and
%                        pili_ext capture for save_sim_config compatibility.
clear; clc;

% Ensure working directory is project root (contains run_simulation.m)
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
if isfile(fullfile(projectRoot, 'run_simulation.m'))
    cd(projectRoot);
elseif isfile('run_simulation.m')
    projectRoot = pwd;
else
    error('batch_repeat_brownian:WrongDirectory', ...
        ['Could not find run_simulation.m. Place this script in batch_sweep/ ', ...
         'under Coupled_Body_Pili, or cd to Coupled_Body_Pili before running.']);
end
addpath(genpath(fullfile(projectRoot, 'src')));

%% ========================================================================
%  USER SETTINGS
%  ========================================================================

% --- Root folder name under output/ (required).  Repeats go to:
%     output/<case_name_base>/001, 002, ...
case_name_base = 'OnePili-S1-L3-Brown';

% --- Number of independent repeats (>= 1)
n_repeats = 100;

% --- RNG (Brownian twister; see run_simulation.m):
%     'fixed'     — per repeat: cfg.brown_rng_seed = rng_seed_base + k (reproducible).
%     'timestamp' — SERIAL: unset cfg.brown_rng_seed; each run uses MMddHHmmss at start.
%                   PARALLEL: pre-set cfg.brown_rng_seed = clock_base + k (unique per worker).
%     'default'   — alias for 'timestamp' (backward-compatible pre-rename value).
brown_rng_seed_mode = 'timestamp';   % 'timestamp' | 'fixed'   ('default' = alias for 'timestamp')

% --- Used only when brown_rng_seed_mode == 'fixed'
rng_seed_base = 202604091;

% --- Parallel execution (Lucien, 2026-05-06) ---
%     use_parallel        : true  -> parfor branch (auto-sized pool, default)
%                           false -> legacy serial for-loop branch (same physics as pre-parallel)
%     mem_per_worker_gb   : conservative per-worker memory budget [GB].
%                           Default 1.5 covers N_body in 200..300 with M_pili <= 2.
%                           Lower this for tighter packing on small-RAM machines.
%     max_resource_frac   : cap on BOTH physical-core count and available memory.
%                           0.75 means we never request > 75% of either resource.
%     worker_blas_threads : per-worker maxNumCompThreads.  1 prevents the
%                           N_workers x N_BLAS_threads "thread storm" that would
%                           otherwise oversubscribe the CPU.  Raise to 2 only if
%                           you have far more cores than concurrent workers.
use_parallel        = true;
mem_per_worker_gb   = 2.0;
max_resource_frac   = 0.75;
worker_blas_threads = 1;

% --- Configuration source (choose one):
%     A)  Set config_json to a path to load via load_sim_config (recommended)
%     B)  Leave config_json empty ('') to use the inline cfg_base block below.

config_json = '';   % e.g. fullfile('configs', 'OnePili-P1-L3.json')  (after cd to project root)

%% ========================================================================
%  INLINE CONFIG (used only when config_json is empty)
%  ========================================================================

cfg_base = struct();
cfg_base.case_name    = '';   % Overwritten per repeat; ignored here
cfg_base.geometry_mode = 'parametric';   % required by save_sim_config (v6+)
cfg_base.minor_axis   = 0.4e-6;
cfg_base.major_axis   = 1.25e-6;
cfg_base.N_body       = 300;
cfg_base.alpha_0      = 0;

pili_cfg(1).type        = 'linear';
pili_cfg(1).length_um   = 3;
pili_cfg(1).site        = 1;
[cfg_base.Pili_Matrix, cfg_base.pili_ext] = parse_pili_config(pili_cfg);

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
cfg_base.tlim   = 3600;
cfg_base.dt     = 0.1;

% Brownian MUST be on (also checked after load)
cfg_base.enable_brown = true;
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

%% ---- Load configuration ----
if isempty(strtrim(case_name_base))
    error('batch_repeat_brownian:EmptyCaseName', ...
        'case_name_base must be a non-empty string (folder name under output/).');
end

if n_repeats < 1 || floor(n_repeats) ~= n_repeats
    error('batch_repeat_brownian:InvalidNRepeats', ...
        'n_repeats must be a positive integer (>= 1).');
end

valid_seed_modes = {'fixed', 'timestamp', 'default'};
if ~any(strcmpi(brown_rng_seed_mode, valid_seed_modes))
    error('batch_repeat_brownian:InvalidSeedMode', ...
        'brown_rng_seed_mode must be ''timestamp'', ''fixed'', or ''default'' (alias for ''timestamp'').');
end
brown_rng_seed_mode = lower(brown_rng_seed_mode);
if strcmp(brown_rng_seed_mode, 'default')
    brown_rng_seed_mode = 'timestamp';   % normalise backward-compatible alias
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

%% ========================================================================
%  PARALLEL SCHEDULING DECISION
%  ========================================================================
%  Decide whether to actually run in parallel and, if so, how many workers
%  to spin up.  The double cap is enforced here (CPU and memory both
%  capped at max_resource_frac of the machine total).  All branches set
%  parallel_info so the header / summary log can report what happened.

parallel_info = struct( ...
    'use_parallel',        use_parallel, ...
    'mem_per_worker_gb',   mem_per_worker_gb, ...
    'max_resource_frac',   max_resource_frac, ...
    'worker_blas_threads', worker_blas_threads, ...
    'n_physical_cores',    NaN, ...
    'available_mem_gb',    NaN, ...
    'cpu_cap',             NaN, ...
    'mem_cap',             NaN, ...
    'binding',             "n/a", ...
    'n_pool',              1, ...
    'fallback_reason',     "");

if use_parallel
    if ~license('test', 'Distrib_Computing_Toolbox')
        warning('batch_repeat_brownian:NoPCT', ...
            'Parallel Computing Toolbox is not licensed on this MATLAB. Falling back to serial.');
        use_parallel = false;
        parallel_info.fallback_reason = "PCT license unavailable";
    end
end

if use_parallel
    [n_pool_decided, parallel_info] = decide_pool_size( ...
        n_repeats, mem_per_worker_gb, max_resource_frac, parallel_info, ...
        worker_blas_threads);
    parallel_info.n_pool = n_pool_decided;
    % n_pool may legitimately be 1 (resource cap or n_repeats==1); still use parfor
    % so the parallel code path stays testable (plan §11).
end

%% ========================================================================
%  PRE-COMPUTE BROWNIAN SEEDS (parallel branch only)
%  ========================================================================
%  In parallel mode we MUST hand each worker a unique seed before parfor
%  starts: brown_rng_seed_from_clock() returns MMddHHmmss and multiple
%  workers booting in the same second would collide.  Both modes therefore
%  produce a deterministic offset list  seeds(k) = base + k  in this branch.
%  The serial branch keeps the original behaviour (no pre-computation;
%  run_simulation generates the timestamp seed itself), so 'fixed' serial
%  runs remain bit-exact with the previous version.

seeds = [];
seed_origin_label = '';
if use_parallel
    seeds = zeros(n_repeats, 1);
    if strcmp(brown_rng_seed_mode, 'fixed')
        for k = 1:n_repeats
            seeds(k) = rng_seed_base + k;
        end
        seed_origin_label = sprintf('fixed (base=%d, +k)', rng_seed_base);
    else
        seed_base = brown_rng_seed_from_clock();
        for k = 1:n_repeats
            seeds(k) = seed_base + k;
        end
        seed_origin_label = sprintf('timestamp+offset (base=%d, +k)', seed_base);
    end
end

%% ========================================================================
%  GET-OR-CREATE PARPOOL  (and configure worker state)
%  ========================================================================

pool_handle = [];
if use_parallel
    pool_handle = ensure_parpool(parallel_info.n_pool);
    if isempty(pool_handle)
        warning('batch_repeat_brownian:PoolCreateFailed', ...
            'parpool creation failed. Falling back to serial.');
        use_parallel = false;
        parallel_info.fallback_reason = "parpool creation failed";
    else
        parallel_info.n_pool = pool_handle.NumWorkers;
        % Force figures invisible on every worker (run_simulation creates
        % a 'Visible','on' summary figure that would otherwise pop up).
        f_vis = parfevalOnAll(pool_handle, @() set(0,'DefaultFigureVisible','off'), 0);
        % Cap per-worker BLAS threads to prevent N_workers x N_blas_threads
        % oversubscription on the host CPU.
        f_thr = parfevalOnAll(pool_handle, @apply_max_num_comp_threads, 0, worker_blas_threads);
        wait(f_vis);
        wait(f_thr);
    end
end
% Worker defaults (figure visibility, BLAS threads) persist for the lifetime of this
% parpool — intentional for batch use.  To restore MATLAB defaults, run delete(gcp('nocreate')).

%% ========================================================================
%  OUTPUT ROOT + SUMMARY LOG
%  ========================================================================

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

%% ---- Header (console) ----
sep = repmat('=', 1, 60);
fprintf('\n%s\n', sep);
fprintf('  Batch repeat (Brownian): %d independent runs\n', n_repeats);
fprintf('  Output root: %s\n', out_root);
if use_parallel
    fprintf('  Mode: PARALLEL  (n_pool = %d, %s-bound)\n', ...
        parallel_info.n_pool, char(parallel_info.binding));
    fprintf(['  Physical cores: %d  (cap %.0f%% -> %d) | ', ...
             'Available memory: %.2f GB  (per-worker %g GB, cap %.0f%% -> %d)\n'], ...
        parallel_info.n_physical_cores, max_resource_frac*100, parallel_info.cpu_cap, ...
        parallel_info.available_mem_gb, mem_per_worker_gb, max_resource_frac*100, parallel_info.mem_cap);
    fprintf('  BLAS threads / worker: %d\n', worker_blas_threads);
else
    if strlength(parallel_info.fallback_reason) > 0
        fprintf('  Mode: SERIAL  (fallback: %s)\n', char(parallel_info.fallback_reason));
    else
        fprintf('  Mode: SERIAL  (use_parallel = false)\n');
    end
end
if strcmp(brown_rng_seed_mode, 'fixed')
    fprintf('  Brownian seeds: FIXED -- cfg.brown_rng_seed = rng_seed_base + k  (base = %d, k = 1..%d)\n', ...
        rng_seed_base, n_repeats);
else
    fprintf('  Brownian seeds: TIMESTAMP -- time-based MMddHHmmss per run\n');
    if use_parallel
        fprintf('                    parallel mode: pre-computed as base + k  (base = %d)\n', ...
            seeds(1) - 1);
    end
end
fprintf('%s\n\n', sep);

%% ---- Header (summary file) ----
if fid_sum ~= -1
    fprintf(fid_sum, 'batch_repeat_brownian\n');
    fprintf(fid_sum, 'case_name_base: %s\n', case_name_base);
    fprintf(fid_sum, 'n_repeats: %d\n', n_repeats);
    fprintf(fid_sum, 'brown_rng_seed_mode: %s\n', brown_rng_seed_mode);
    if strcmpi(brown_rng_seed_mode, 'fixed')
        fprintf(fid_sum, 'rng_seed_base: %d\n', rng_seed_base);
    end
    fprintf(fid_sum, 'execution_mode: %s\n', char(string(use_parallel)));
    if use_parallel
        fprintf(fid_sum, 'n_pool: %d\n',               parallel_info.n_pool);
        fprintf(fid_sum, 'n_physical_cores: %d\n',     parallel_info.n_physical_cores);
        if isfield(parallel_info, 'core_probe')
            fprintf(fid_sum, 'core_probe: %s\n',       char(parallel_info.core_probe));
        end
        fprintf(fid_sum, 'available_mem_gb: %.2f\n',   parallel_info.available_mem_gb);
        fprintf(fid_sum, 'mem_per_worker_gb: %g\n',    mem_per_worker_gb);
        fprintf(fid_sum, 'max_resource_frac: %g\n',    max_resource_frac);
        fprintf(fid_sum, 'cpu_cap: %d\n',              parallel_info.cpu_cap);
        fprintf(fid_sum, 'mem_cap: %d\n',              parallel_info.mem_cap);
        fprintf(fid_sum, 'binding: %s\n',              char(parallel_info.binding));
        fprintf(fid_sum, 'worker_blas_threads: %d\n',  worker_blas_threads);
        fprintf(fid_sum, 'seed_origin: %s\n',          seed_origin_label);
    elseif strlength(parallel_info.fallback_reason) > 0
        fprintf(fid_sum, 'parallel_fallback_reason: %s\n', char(parallel_info.fallback_reason));
    end
    fprintf(fid_sum, 'started: %s\n\n', log_ts);
end

%% ========================================================================
%  SLICE-PRE-ALLOCATED RESULT ARRAYS (parfor requires this layout)
%  ========================================================================
run_status = repmat("PENDING", n_repeats, 1);
run_outdir = repmat("",        n_repeats, 1);
run_errmsg = repmat("",        n_repeats, 1);
run_wall   = nan(n_repeats, 1);
run_seed   = zeros(n_repeats, 1);   % 0 = unknown (timestamp serial branch)

t_batch = tic;

%% ========================================================================
%  EXECUTION
%  ========================================================================

if use_parallel
    %% ---- PARALLEL BRANCH (parfor + DataQueue progress) ----
    progress_q = parallel.pool.DataQueue;
    print_progress(struct('reset', true), n_repeats);   % reset persistent counter
    afterEach(progress_q, @(msg) print_progress(msg, n_repeats));

    fprintf('Launching parfor with %d workers...\n\n', parallel_info.n_pool);

    parfor k = 1:n_repeats
        cfg_k                 = cfg;                                         %#ok<PFBNS>
        cfg_k.case_name       = fullfile(case_name_base, sprintf('%03d', k));
        cfg_k.brown_rng_seed  = seeds(k);
        cfg_k.preview_initial = false;
        cfg_k.auto_confirm    = true;
        run_seed(k) = seeds(k);

        t_run = tic;
        try
            [res_k, ~, dir_k] = run_simulation(cfg_k);
            wall_k = toc(t_run);
            run_wall(k) = wall_k;

            if isempty(res_k)
                run_status(k) = "FAILED";
                run_errmsg(k) = "run_simulation returned empty results";
                run_outdir(k) = string(dir_k);
                send(progress_q, struct('k', k, 'status', "FAILED", ...
                    'msg', 'empty results', 'wall', wall_k));
            else
                run_status(k) = "OK";
                run_outdir(k) = string(dir_k);
                send(progress_q, struct('k', k, 'status', "OK", ...
                    'msg', '', 'wall', wall_k));
            end
        catch ME
            wall_k = toc(t_run);
            run_wall(k)   = wall_k;
            run_status(k) = "FAILED";
            run_errmsg(k) = string(ME.message);
            run_outdir(k) = "";
            send(progress_q, struct('k', k, 'status', "FAILED", ...
                'msg', ME.message, 'wall', wall_k));
        end
    end

else
    %% ---- SERIAL BRANCH (legacy for-loop; summary lines written after the loop) ----
    for k = 1:n_repeats
        sub = sprintf('%03d', k);
        fprintf('%s\n  Repeat %d/%d  ->  output/%s/%s\n%s\n', ...
            repmat('-', 1, 56), k, n_repeats, case_name_base, sub, repmat('-', 1, 56));

        cfg_k = cfg;
        cfg_k.case_name = fullfile(case_name_base, sub);

        if strcmp(brown_rng_seed_mode, 'fixed')
            cfg_k.brown_rng_seed = rng_seed_base + k;
            run_seed(k) = rng_seed_base + k;
        else
            % Timestamp policy: let run_simulation use brown_rng_seed_from_clock (MMddHHmmss).
            if isfield(cfg_k, 'brown_rng_seed')
                cfg_k = rmfield(cfg_k, 'brown_rng_seed');
            end
        end

        res_k = [];
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

            if ~strcmp(brown_rng_seed_mode, 'fixed') && ...
                    ~isempty(res_k) && isfield(res_k, 'brown_rng_seed')
                run_seed(k) = res_k.brown_rng_seed;
            end
        catch ME
            run_wall(k)   = toc(t_run);
            run_status(k) = "FAILED";
            run_errmsg(k) = string(ME.message);
            run_outdir(k) = "";
            fprintf('  [ERR]  %s\n', ME.message);
        end
    end
end

elapsed = toc(t_batch);

%% ========================================================================
%  AGGREGATE SUMMARY WRITES (main thread, in k-order)
%  ========================================================================
%  Both branches converge here so the summary file format stays identical
%  regardless of execution mode.  In serial mode the per-iteration
%  fprintf to fid_sum has been moved here for symmetry; the printed
%  console output (legacy [OK]/[ERR] lines) is preserved inside the loop.

if fid_sum ~= -1
    for k = 1:n_repeats
        if run_seed(k) > 0
            seed_tag = sprintf('seed=%d', run_seed(k));
            if ~strcmp(brown_rng_seed_mode, 'fixed')
                seed_tag = [seed_tag, ' (timestamp)']; %#ok<AGROW>
            end
        else
            seed_tag = 'seed=n/a';
        end
        fprintf(fid_sum, '%03d  %s  status=%s  wall=%.2f s\n', ...
            k, seed_tag, run_status(k), run_wall(k));

        if run_status(k) == "FAILED" && strlength(run_errmsg(k)) > 0
            fprintf(fid_sum, '  error: %s\n', run_errmsg(k));
        elseif strlength(run_outdir(k)) > 0
            fprintf(fid_sum, '  out: %s\n', run_outdir(k));
        end
    end
end

%% ========================================================================
%  FINAL REPORT
%  ========================================================================

n_ok = sum(run_status == "OK");

fprintf('\n%s\n', sep);
fprintf('  Batch complete: %.1f s total\n', elapsed);
fprintf('  Success: %d / %d\n', n_ok, n_repeats);
if use_parallel
    sum_wall = sum(run_wall(~isnan(run_wall)));
    speedup  = sum_wall / max(elapsed, eps);
    fprintf('  Parallel scheduling report:\n');
    fprintf('    n_pool                  : %d\n',          parallel_info.n_pool);
    fprintf('    physical cores detected : %d\n',          parallel_info.n_physical_cores);
    fprintf('    available memory (start): %.2f GB\n',     parallel_info.available_mem_gb);
    fprintf('    cap                     : %.0f%% (cpu_cap=%d, mem_cap=%d, binding=%s)\n', ...
        max_resource_frac*100, parallel_info.cpu_cap, parallel_info.mem_cap, char(parallel_info.binding));
    fprintf('    sum of per-case wall    : %.1f s\n',      sum_wall);
    fprintf('    parallel speed-up est.  : %.2fx\n',       speedup);
end
fprintf('%s\n', sep);

if fid_sum ~= -1
    fprintf(fid_sum, '\nelapsed_total_s: %.2f\n', elapsed);
    fprintf(fid_sum, 'success_count: %d / %d\n', n_ok, n_repeats);
    if use_parallel
        sum_wall = sum(run_wall(~isnan(run_wall)));
        fprintf(fid_sum, 'sum_per_case_wall_s: %.2f\n',       sum_wall);
        fprintf(fid_sum, 'parallel_speedup_estimate: %.2f\n', sum_wall / max(elapsed, eps));
    end
    fclose(fid_sum);
    fprintf('\n  Summary written: %s\n', sum_path);
end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function apply_max_num_comp_threads(n)
%APPLY_MAX_NUM_COMP_THREADS  Worker-side wrapper for parfevalOnAll BLAS cap.
    maxNumCompThreads(n);
end


function [n_pool, info] = decide_pool_size(n_repeats, mem_per_worker_gb, max_frac, info, blas_threads)
%DECIDE_POOL_SIZE  Probe machine resources and compute capped worker count.
%
%  blas_threads is recorded in the console probe (plan: pass worker BLAS policy
%  into the sizing routine).  Optional 6th-arg logf was omitted — fprintf here
%  matches every other batch script in this project.
%
%  Returns the smaller of:
%    - floor(n_physical_cores * max_frac)
%    - floor(available_mem_gb * max_frac / mem_per_worker_gb)
%    - n_repeats
%  and labels which constraint binds (used for the diagnostic header).

    if nargin < 5 || isempty(blas_threads)
        blas_threads = 1;
    end

    %% --- Physical core count (avoid hyper-threading inflation) ---
    n_phys = NaN;
    probe_note = "";
    try
        if ispc
            [stat, out] = system('wmic cpu get NumberOfCores /value');
            if stat == 0
                tok = regexp(out, 'NumberOfCores=(\d+)', 'tokens');
                if ~isempty(tok)
                    cores_per_cpu = cellfun(@(t) str2double(t{1}), tok);
                    n_phys = sum(cores_per_cpu);
                    probe_note = "wmic";
                end
            end
        end
    catch
        n_phys = NaN;
    end

    if ~isfinite(n_phys) || n_phys < 1
        try
            if ispc
                ps_cmd = ['powershell -NoProfile -Command ' ...
                    '"(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum"'];
                [st_ps, out_ps] = system(ps_cmd);
                if st_ps == 0
                    n_cim = str2double(strtrim(out_ps));
                    if isfinite(n_cim) && n_cim >= 1
                        n_phys = n_cim;
                        probe_note = "CIM";
                    end
                end
            end
        catch
        end
    end

    if ~isfinite(n_phys) || n_phys < 1
        try
            n_env = str2double(getenv('NUMBER_OF_PROCESSORS'));
            if isfinite(n_env) && n_env >= 1
                % Windows: logical CPU count — halve as a conservative physical estimate.
                n_phys = max(1, floor(n_env / 2));
                probe_note = "NUMBER_OF_PROCESSORS/2";
            end
        catch
        end
    end

    if ~isfinite(n_phys) || n_phys < 1
        try
            n_logical = feature('numCores');
            n_phys = max(1, floor(n_logical / 2));
            probe_note = "feature(numCores)/2";
        catch
            try
                n_phys = max(1, maxNumCompThreads);
                probe_note = "maxNumCompThreads";
            catch
                n_phys = 1;
                probe_note = "fallback=1";
            end
        end
    end

    info.n_physical_cores = n_phys;
    info.core_probe = probe_note;

    %% --- Available physical memory (Windows: memory()) ---
    avail_gb = NaN;
    try
        [~, sys] = memory;
        avail_gb = sys.PhysicalMemory.Available / 1e9;
    catch
        avail_gb = NaN;
    end
    if ~isfinite(avail_gb)
        avail_gb = 4.0;
        fprintf('  [WARN] memory() unavailable; assuming 4 GB available.\n');
    end
    info.available_mem_gb = avail_gb;

    %% --- Two caps + n_repeats ---
    n_cpu = max(1, floor(n_phys   * max_frac));
    n_mem = max(1, floor(avail_gb * max_frac / mem_per_worker_gb));
    info.cpu_cap = n_cpu;
    info.mem_cap = n_mem;

    n_pool = min([n_cpu, n_mem, n_repeats]);

    if n_pool == n_repeats && n_repeats <= n_cpu && n_repeats <= n_mem
        info.binding = "n_repeats";
    elseif n_cpu <= n_mem
        info.binding = "cpu";
    else
        info.binding = "mem";
    end

    fprintf('\n  [POOL] Resource probe:\n');
    fprintf('         physical cores : %d  (%s)  (cap %.0f%% -> %d workers)\n', ...
        n_phys, char(probe_note), max_frac*100, n_cpu);
    fprintf('         available mem  : %.2f GB  (per-worker %g GB, cap %.0f%% -> %d workers)\n', ...
        avail_gb, mem_per_worker_gb, max_frac*100, n_mem);
    fprintf('         n_repeats      : %d\n', n_repeats);
    fprintf('         BLAS threads   : %d / worker (after parfevalOnAll init)\n', blas_threads);
    fprintf('         => n_pool = %d  (binding: %s)\n\n', n_pool, char(info.binding));
end


function pool = ensure_parpool(n_target)
%ENSURE_PARPOOL  Reuse-or-create a process-based parpool of the requested size.
%
%  Tries 'Processes' (R2022a+) first, falls back to 'local' on older MATLAB.
%  Returns [] if pool creation fails.

    pool = [];
    p = gcp('nocreate');
    if ~isempty(p) && p.NumWorkers == n_target
        pool = p;
        return;
    end
    if ~isempty(p)
        try %#ok<TRYNC>
            delete(p);
        end
    end
    try
        try
            pool = parpool('Processes', n_target);
        catch
            pool = parpool('local', n_target);
        end
    catch ME
        fprintf('  [WARN] parpool creation failed: %s\n', ME.message);
        pool = [];
    end
end


function print_progress(msg, n_total)
%PRINT_PROGRESS  Streaming progress callback for parfor's DataQueue.
%
%  Maintains a persistent (per-MATLAB-session) counter and timer so that
%  workers can fire-and-forget completion notices and the main thread
%  prints them in completion order with elapsed/ETA estimates.
%
%  Reset between batches with:  print_progress(struct('reset', true), n_total)

    persistent n_done t0;

    if isstruct(msg) && isfield(msg, 'reset') && msg.reset
        n_done = 0;
        t0     = tic;
        return;
    end
    if isempty(n_done) || isempty(t0)
        n_done = 0;
        t0     = tic;
    end
    if ~isstruct(msg) || ~isfield(msg, 'k')
        return;   % defensive: ignore unexpected messages
    end

    n_done  = n_done + 1;
    elapsed = toc(t0);
    eta_s   = elapsed * max(0, n_total - n_done) / max(n_done, 1);

    if msg.status == "OK"
        fprintf('  [%3d/%3d] OK     k=%03d  wall=%6.1f s  | elapsed=%6.1f s  ETA=%6.1f s\n', ...
            n_done, n_total, msg.k, msg.wall, elapsed, eta_s);
    else
        err_one_line = regexprep(char(msg.msg), '[\r\n]+', ' ');
        if numel(err_one_line) > 120
            err_one_line = [err_one_line(1:117), '...'];
        end
        fprintf('  [%3d/%3d] FAILED k=%03d  wall=%6.1f s  | err: %s\n', ...
            n_done, n_total, msg.k, msg.wall, err_one_line);
    end
end


