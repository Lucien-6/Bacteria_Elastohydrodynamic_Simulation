%% batch_repeat_brownian_HPC_Linux.m
%  HPC/Linux fork of batch_repeat_brownian.m -- independent repeated
%  simulations (same configuration, Brownian motion ON), hardened for
%  headless, single-node-exclusive Linux HPC cluster execution.
%
%  This file is a maintained fork of batch_repeat_brownian.m, NOT
%  auto-synced with it -- port any future physics/batch-logic fix from the
%  original here manually.  Differences from batch_repeat_brownian.m:
%    1. Calls run_simulation_HPC_Linux (not run_simulation) -- see that
%       file's header for why (VideoWriter 'MPEG-4' is Windows/Mac-only
%       and hard-errors on Linux; the HPC_Linux copy uses 'Motion JPEG
%       AVI' on Linux instead).
%    2. decide_pool_size gained a Linux-native physical-core probe
%       (lscpu) and a Linux-native available-memory probe (/proc/meminfo
%       MemAvailable), and the previous feature('numCores')/2 fallback
%       -- which incorrectly halved a value that is already a physical
%       (not logical) core count -- is fixed to no longer halve it.
%    3. The single max_resource_frac cap is split into cpu_resource_frac
%       and mem_resource_frac (memory is more conservative by default),
%       see USER SETTINGS below for why.
%    4. ensure_parpool relocates the local cluster's JobStorageLocation to
%       a node-local, job-unique temp folder on Linux (see that function's
%       header) instead of the MATLAB default under the user's home
%       directory, which on clusters with distributed/network home storage
%       is a documented cause of slow pool startup and, when concurrent
%       jobs share it, race conditions / pool-start failures.
%    5. generate_animation / compute_flowfield are forced OFF for every
%       repeat except a small sample (first, middle, last by default --
%       see heavy_output_sample_idx below), because these features hold
%       significant extra memory per worker (frame buffers, video
%       encoding, flow-field grids) and were implicated in an observed
%       failure mode on a real HPC run: workers silently killed (Linux
%       OOM-killer signature: process count in `top` suddenly drops,
%       remaining processes idle) partway through a long parfor batch,
%       after which the batch hangs indefinitely (parfor has no built-in
%       timeout / dead-worker recovery -- the client just waits forever
%       for a result that will never arrive).  Shrinking per-worker memory
%       footprint reduces the likelihood of hitting this; it does not
%       eliminate it (a full fix would replace parfor with a
%       parfeval+fetchNext task queue with explicit per-task timeouts --
%       deliberately NOT done in this fork per 2026-07-03 discussion;
%       revisit if OOM/hangs recur after this fix).
%    6. Live sub-case progress: run_simulation_HPC_Linux accepts optional
%       cfg.progress_queue / cfg.progress_tag and, when set, sends a
%       throttled (>= 15 s apart per repeat) 'heartbeat' progress message
%       through the same DataQueue used for completion messages -- this
%       is the only way to see a still-running repeat's progress while
%       cfg.silent = true suppresses all of its normal console echo.
%       dispatch_progress_msg (renamed from on_run_complete) routes
%       'heartbeat' vs 'done' messages; heartbeats are console-only
%       (never written to the persistent summary file).
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
%            'timestamp' — SERIAL branch: omit cfg.brown_rng_seed; run_simulation_HPC_Linux uses
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
%  Updated: 2026-07-03 — Lucien — summary reproducibility (config source,
%                        batch_id snapshot, brown_rng_policy, fallback probe);
%                        FAILED outdir/log_hint; footer finished/failed/pending;
%                        unified empty-results msg; catch partial outdir.
%  Forked:  2026-07-03 — Lucien — HPC_Linux variant per WX2 partition review
%                        (Linux, MATLAB R2023a, exclusive single node,
%                        distributed/network storage): see items 1-5 above.
clear; clc;

% Ensure working directory is project root (contains run_simulation_HPC_Linux.m)
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
if isfile(fullfile(projectRoot, 'run_simulation_HPC_Linux.m'))
    cd(projectRoot);
elseif isfile('run_simulation_HPC_Linux.m')
    projectRoot = pwd;
else
    error('batch_repeat_brownian_HPC_Linux:WrongDirectory', ...
        ['Could not find run_simulation_HPC_Linux.m. Place this script in batch_sweep/ ', ...
         'under Coupled_Body_Pili, or cd to Coupled_Body_Pili before running.']);
end
addpath(genpath(fullfile(projectRoot, 'src')));

%% ========================================================================
%  USER SETTINGS
%  ========================================================================

% --- Root folder name under output/ (required).  Repeats go to:
%     output/<case_name_base>/001, 002, ...
case_name_base = 'TwoPili-S29-L3-Brown';

% --- Number of independent repeats (>= 1)
n_repeats = 100;

% --- RNG (Brownian twister; see run_simulation_HPC_Linux.m):
%     'fixed'     — per repeat: cfg.brown_rng_seed = rng_seed_base + k (reproducible).
%     'timestamp' — SERIAL: unset cfg.brown_rng_seed; each run uses MMddHHmmss at start.
%                   PARALLEL: pre-set cfg.brown_rng_seed = clock_base + k (unique per worker).
%     'default'   — alias for 'timestamp' (backward-compatible pre-rename value).
brown_rng_seed_mode = 'fixed';   % 'timestamp' | 'fixed'   ('default' = alias for 'timestamp')

% --- Used only when brown_rng_seed_mode == 'fixed'
rng_seed_base = 202607021;

% --- Parallel execution (Lucien, 2026-05-06; HPC_Linux hardening 2026-07-03) ---
%     use_parallel         : true  -> parfor branch (auto-sized pool, default)
%                            false -> legacy serial for-loop branch (same physics as pre-parallel)
%     mem_per_worker_gb    : conservative per-worker memory budget [GB].
%                            IMPORTANT: this default was NOT empirically measured
%                            for this exact configuration (N_body=300, animation +
%                            flowfield enabled for sampled repeats).  Before trusting
%                            the auto-sized pool on a real cluster allocation, run ONE
%                            standalone repeat and measure its real peak RSS (e.g.
%                            `/usr/bin/time -v matlab -batch ...` or `top`/`ps` on
%                            Linux), then set this to that measured value + ~30-50%
%                            margin.  Under-budgeting this is the most likely root
%                            cause of the "workers vanish mid-run, batch hangs
%                            forever" failure mode (Linux OOM-killer; parfor has no
%                            timeout / dead-worker recovery -- see file header).
%     cpu_resource_frac    : cap on physical-core count only (0.90 -> never request
%                            more than 90% of detected physical cores).
%     mem_resource_frac    : cap on available memory only, kept separate from and
%                            more conservative than cpu_resource_frac by default --
%                            per-run memory is harder to predict than CPU usage, and
%                            running out of RAM is the failure mode that hangs the
%                            whole batch, whereas running out of "spare" CPU just
%                            makes things slower.
%     worker_blas_threads  : per-worker maxNumCompThreads.  1 prevents the
%                            N_workers x N_BLAS_threads "thread storm" that would
%                            otherwise oversubscribe the CPU.  On Linux, ALSO set
%                            the equivalent OS environment variables
%                            (OMP_NUM_THREADS=1, MKL_NUM_THREADS=1,
%                            OPENBLAS_NUM_THREADS=1) in your job submission script
%                            BEFORE launching MATLAB: maxNumCompThreads alone does
%                            not always fully constrain every underlying threaded
%                            library on Linux.
use_parallel         = true;
mem_per_worker_gb    = 2.5;
cpu_resource_frac    = 0.90;
mem_resource_frac    = 0.65;
worker_blas_threads  = 1;

% --- Heavy-output sampling (HPC_Linux hardening) ---
%     generate_animation / compute_flowfield hold significant extra per-worker
%     memory (frame buffers, video encoding, flow-field grids) and are usually
%     only needed for a visual sanity check, not for every single repeat.  To
%     shrink peak per-worker memory footprint (reducing OOM risk -- see file
%     header), only the repeats listed in heavy_output_sample_idx (computed
%     below from n_repeats; default = first, middle, last) get the full
%     animation/flowfield output; every other repeat gets numeric-only output
%     regardless of what cfg_base / the JSON config specifies.  Set
%     heavy_output_sample_idx = 1:n_repeats to disable sampling (full output
%     every repeat, original behaviour) once memory headroom is confirmed safe.
heavy_output_sample_idx = [];   % placeholder; computed after n_repeats validation below

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
pili_cfg(1).site        = 2;
pili_cfg(2).type        = 'linear';
pili_cfg(2).length_um   = 3;
pili_cfg(2).site        = 9;
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
cfg_base.generate_animation = true;   % template value; per-repeat override below
                                       % forces this OFF except for
                                       % heavy_output_sample_idx repeats (see
                                       % USER SETTINGS -> Heavy-output sampling)
cfg_base.animation_fps      = 30;
cfg_base.animation_skip     = 5;

cfg_base.auto_confirm = true;

cfg_base.flowfield_grid_spacing = 0.05;
cfg_base.flowfield_grid_range   = 10.0;
cfg_base.flowfield_quiver_skip  = 5;
cfg_base.flowfield_chunk_size   = 2000;

%% ---- Load configuration ----
if isempty(strtrim(case_name_base))
    error('batch_repeat_brownian_HPC_Linux:EmptyCaseName', ...
        'case_name_base must be a non-empty string (folder name under output/).');
end

if n_repeats < 1 || floor(n_repeats) ~= n_repeats
    error('batch_repeat_brownian_HPC_Linux:InvalidNRepeats', ...
        'n_repeats must be a positive integer (>= 1).');
end

% --- Resolve heavy-output sampling indices (first, middle, last by default) ---
heavy_output_sample_idx = unique([1, ceil(n_repeats / 2), n_repeats]);
heavy_output_sample_idx = heavy_output_sample_idx( ...
    heavy_output_sample_idx >= 1 & heavy_output_sample_idx <= n_repeats);

valid_seed_modes = {'fixed', 'timestamp', 'default'};
if ~any(strcmpi(brown_rng_seed_mode, valid_seed_modes))
    error('batch_repeat_brownian_HPC_Linux:InvalidSeedMode', ...
        'brown_rng_seed_mode must be ''timestamp'', ''fixed'', or ''default'' (alias for ''timestamp'').');
end
brown_rng_seed_mode = lower(brown_rng_seed_mode);
if strcmp(brown_rng_seed_mode, 'default')
    brown_rng_seed_mode = 'timestamp';   % normalise backward-compatible alias
end

if ~isempty(strtrim(config_json))
    if ~isfile(config_json)
        error('batch_repeat_brownian_HPC_Linux:ConfigNotFound', ...
            'Config file not found:\n  %s', config_json);
    end
    cfg = load_sim_config(config_json);
else
    cfg = cfg_base;
end

% --- Require Brownian ON ---
if ~isfield(cfg, 'enable_brown') || ~logical(cfg.enable_brown)
    error('batch_repeat_brownian_HPC_Linux:BrownianRequired', ...
        ['This script is only for Brownian simulations.  Set enable_brown = true ', ...
         '(inline) or brownian.enable = true in JSON.']);
end

% --- Batch defaults: no preview, headless ---
cfg.preview_initial = false;
cfg.auto_confirm    = true;

% --- Config source metadata (written to batch summary header) ---
config_source = struct( ...
    'type', 'inline', ...
    'json_path', '', ...
    'note', 'inline cfg_base block in batch_repeat_brownian_HPC_Linux.m');
if ~isempty(strtrim(config_json))
    config_source.type = 'json';
    if isfile(config_json)
        config_source.json_path = config_json;
    else
        config_source.json_path = fullfile(pwd, config_json);
    end
end

%% ========================================================================
%  PARALLEL SCHEDULING DECISION
%  ========================================================================
%  Decide whether to actually run in parallel and, if so, how many workers
%  to spin up.  CPU and memory are capped independently (cpu_resource_frac,
%  mem_resource_frac -- see USER SETTINGS for why they are separate).  All
%  branches set parallel_info so the header / summary log can report what
%  happened.

parallel_info = struct( ...
    'use_parallel',        use_parallel, ...
    'mem_per_worker_gb',   mem_per_worker_gb, ...
    'cpu_resource_frac',   cpu_resource_frac, ...
    'mem_resource_frac',   mem_resource_frac, ...
    'worker_blas_threads', worker_blas_threads, ...
    'n_physical_cores',    NaN, ...
    'available_mem_gb',    NaN, ...
    'cpu_cap',             NaN, ...
    'mem_cap',             NaN, ...
    'binding',             "n/a", ...
    'n_pool',              1, ...
    'fallback_reason',     "");

parallel_requested = use_parallel;

if use_parallel
    if ~license('test', 'Distrib_Computing_Toolbox')
        warning('batch_repeat_brownian_HPC_Linux:NoPCT', ...
            'Parallel Computing Toolbox is not licensed on this MATLAB. Falling back to serial.');
        use_parallel = false;
        parallel_info.fallback_reason = "PCT license unavailable";
    end
end

if use_parallel
    [n_pool_decided, parallel_info] = decide_pool_size( ...
        n_repeats, mem_per_worker_gb, cpu_resource_frac, mem_resource_frac, parallel_info, ...
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
%  run_simulation_HPC_Linux generates the timestamp seed itself), so 'fixed' serial
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
        warning('batch_repeat_brownian_HPC_Linux:PoolCreateFailed', ...
            'parpool creation failed. Falling back to serial.');
        use_parallel = false;
        parallel_info.fallback_reason = "parpool creation failed";
    else
        parallel_info.n_pool = pool_handle.NumWorkers;
    end
end
% Unconditionally (re)apply worker-side defaults whenever we have a pool.
% This covers the case where ensure_parpool REUSED an existing pool created
% by another script — without this, that pool's workers might still have
% default BLAS thread counts (= physical cores), causing an N_workers x
% N_cores thread storm.  parfevalOnAll is idempotent for both settings,
% so reapplying on a freshly-created pool is harmless.
if use_parallel && ~isempty(pool_handle)
    f_vis = parfevalOnAll(pool_handle, @() set(0,'DefaultFigureVisible','off'), 0);
    f_thr = parfevalOnAll(pool_handle, @apply_max_num_comp_threads, 0, worker_blas_threads);
    wait(f_vis);
    wait(f_thr);
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

batch_start_dt   = datetime('now');
% NOTE: char(datetimeObj, 'Format', fmt) is NOT a valid signature -- char()
% only accepts a positional format string, not a 'Format' name/value pair
% (that name/value pair is only valid on the datetime() constructor itself).
% Passing it here made MATLAB try to parse the literal text 'Format' as the
% format pattern, which fails on the unsupported token 'F' (confirmed
% reproducible on R2023a AND R2026a -- this was a latent bug, not a
% version-specific issue).  Fixed to use the positional-format signature.
log_ts           = char(batch_start_dt, 'uuuuMMdd''_''HHmmss');
started_human    = char(batch_start_dt, 'uuuu-MM-dd HH:mm:ss');
batch_snapshot_path = fullfile(out_root, sprintf('batch_config_snapshot_%s.json', log_ts));
batch_snapshot_log  = '';
cfg_snapshot = cfg;
cfg_snapshot.case_name = case_name_base;
try
    save_sim_config(cfg_snapshot, batch_snapshot_path);
    batch_snapshot_log = path_for_log(batch_snapshot_path);
catch ME_snap
    warning('batch_repeat_brownian_HPC_Linux:ConfigSnapshotFailed', ...
        'Could not save batch config snapshot: %s', ME_snap.message);
end

sum_path = fullfile(out_root, sprintf('batch_repeat_summary_%s.txt', log_ts));
fid_sum  = fopen(sum_path, 'w', 'n', 'UTF-8');
if fid_sum == -1
    warning('batch_repeat_brownian_HPC_Linux:SummaryOpenFailed', ...
        'Could not open summary file for writing:\n  %s\n  Per-run records will not be persisted.', sum_path);
    sum_path = '';
end

%% ---- Header (console) ----
sep = repmat('=', 1, 60);
fprintf('\n%s\n', sep);
fprintf('  Batch repeat (Brownian): %d independent runs\n', n_repeats);
fprintf('  Output root: %s\n', out_root);
if ~isempty(batch_snapshot_log)
    fprintf('  Config snapshot: %s\n', batch_snapshot_log);
end
if use_parallel
    fprintf('  Mode: PARALLEL  (n_pool = %d, %s-bound)\n', ...
        parallel_info.n_pool, char(parallel_info.binding));
    fprintf(['  Physical cores: %d  (cap %.0f%% -> %d) | ', ...
             'Available memory: %.2f GB  (per-worker %g GB, cap %.0f%% -> %d)\n'], ...
        parallel_info.n_physical_cores, cpu_resource_frac*100, parallel_info.cpu_cap, ...
        parallel_info.available_mem_gb, mem_per_worker_gb, mem_resource_frac*100, parallel_info.mem_cap);
    if isfield(parallel_info, 'core_probe') && isfield(parallel_info, 'mem_probe')
        fprintf('  Core probe: %s | Memory probe: %s\n', ...
            char(parallel_info.core_probe), char(parallel_info.mem_probe));
    end
    fprintf('  BLAS threads / worker: %d\n', worker_blas_threads);
    fprintf('  Heavy output (animation/flowfield) sampled at repeats: %s\n', ...
        strjoin(arrayfun(@(k) sprintf('%03d', k), heavy_output_sample_idx, 'UniformOutput', false), ', '));
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
    fprintf(fid_sum, 'batch_repeat_brownian_HPC_Linux\n');
    fprintf(fid_sum, 'batch_id: %s\n', log_ts);
    fprintf(fid_sum, 'started: %s\n', started_human);
    fprintf(fid_sum, 'case_name_base: %s\n', case_name_base);
    fprintf(fid_sum, 'n_repeats: %d\n', n_repeats);
    fprintf(fid_sum, 'config_source_type: %s\n', config_source.type);
    if strcmp(config_source.type, 'json')
        fprintf(fid_sum, 'config_json: %s\n', path_for_log(config_source.json_path));
    else
        fprintf(fid_sum, 'config_source_note: %s\n', config_source.note);
    end
    if ~isempty(batch_snapshot_log)
        fprintf(fid_sum, 'batch_config_snapshot: %s\n', batch_snapshot_log);
    end
    fprintf(fid_sum, 'project_root: %s\n', path_for_log(projectRoot));
    fprintf(fid_sum, 'script_path: %s\n', path_for_log(mfilename('fullpath')));
    fprintf(fid_sum, 'brown_rng_seed_mode: %s\n', brown_rng_seed_mode);
    if strcmpi(brown_rng_seed_mode, 'fixed')
        fprintf(fid_sum, 'rng_seed_base: %d\n', rng_seed_base);
    end
    fprintf(fid_sum, 'parallel_requested: %s\n', char(string(parallel_requested)));
    fprintf(fid_sum, 'execution_mode: %s\n', char(string(use_parallel)));
    fprintf(fid_sum, 'brown_rng_policy: %s\n', ...
        brown_rng_policy_line(brown_rng_seed_mode, use_parallel, rng_seed_base, n_repeats, seeds));
    if use_parallel
        fprintf(fid_sum, 'summary_record_order: completion\n');
    else
        fprintf(fid_sum, 'summary_record_order: sequential\n');
    end
    resource_probed = isfinite(parallel_info.n_physical_cores);
    seeds_precomputed_unused = ~isempty(seeds) && ~use_parallel;
    if resource_probed
        write_resource_probe_summary(fid_sum, parallel_info, mem_per_worker_gb, ...
            cpu_resource_frac, mem_resource_frac, worker_blas_threads);
    end
    if use_parallel
        fprintf(fid_sum, 'n_pool: %d\n', parallel_info.n_pool);
        fprintf(fid_sum, 'seed_origin: %s\n', seed_origin_label);
    end
    if strlength(parallel_info.fallback_reason) > 0
        fprintf(fid_sum, 'parallel_fallback_reason: %s\n', char(parallel_info.fallback_reason));
    end
    if seeds_precomputed_unused
        fprintf(fid_sum, 'seeds_precomputed_unused: true\n');
    end
    fprintf(fid_sum, 'heavy_output_sample_idx: %s\n', ...
        strjoin(arrayfun(@(k) sprintf('%03d', k), heavy_output_sample_idx, 'UniformOutput', false), ', '));
    fprintf(fid_sum, '\n');
    fclose(fid_sum);   % flush header to disk; per-run records appended incrementally
    fid_sum = -1;      % marker: header is closed; use append_run_summary for per-run writes
end

%% ========================================================================
%  SLICE-PRE-ALLOCATED RESULT ARRAYS (parfor requires this layout)
%  ========================================================================
run_status = repmat("PENDING", n_repeats, 1);
run_outdir = repmat("",        n_repeats, 1);
run_errmsg = repmat("",        n_repeats, 1);
run_errid  = repmat("",        n_repeats, 1);   % ME.identifier for failed runs
run_errloc = repmat("",        n_repeats, 1);   % 'file (line N)' for failed runs
run_wall   = nan(n_repeats, 1);
run_seed   = zeros(n_repeats, 1);   % 0 = unknown (timestamp serial branch)

t_batch = tic;

%% ========================================================================
%  EXECUTION
%  ========================================================================
%  Wrapped in an outer try/catch so that a fatal, uncaught escape (e.g. a
%  worker process crash, or an error thrown outside the per-iteration
%  try/catch below) still falls through to the FINAL REPORT section
%  instead of aborting the script silently before any report is printed
%  or the summary footer is written.  Per-iteration failures are already
%  handled locally and never reach here; this only catches whole-batch
%  aborts.  Sliced result arrays (run_status, run_outdir, ...) retain
%  whatever partial data had already been committed for iterations that
%  completed before the abort -- MATLAB's parfor commits sliced-output
%  writes incrementally as each iteration returns, it does not defer all
%  of them to loop completion (verified empirically).

batch_fatal_error = [];
try

if use_parallel
    %% ---- PARALLEL BRANCH (parfor + DataQueue progress) ----
    progress_q = parallel.pool.DataQueue;
    print_progress(struct('reset', true), n_repeats);   % reset persistent counter
    afterEach(progress_q, @(msg) dispatch_progress_msg(msg, n_repeats, sum_path, brown_rng_seed_mode));

    fprintf('Launching parfor with %d workers...\n\n', parallel_info.n_pool);

    parfor k = 1:n_repeats
        cfg_k                 = cfg;                                         %#ok<PFBNS>
        cfg_k.case_name       = fullfile(case_name_base, sprintf('%03d', k));
        cfg_k.brown_rng_seed  = seeds(k);
        cfg_k.preview_initial = false;
        cfg_k.auto_confirm    = true;
        cfg_k.silent          = true;   % suppress per-run stdout on workers (see run_simulation_HPC_Linux.dual_log)
        % Heavy-output sampling: force animation/flowfield OFF except for the
        % sampled repeats, to shrink per-worker peak memory (see USER SETTINGS).
        is_sample_k                = ismember(k, heavy_output_sample_idx);   %#ok<PFBNS>
        cfg_k.generate_animation   = is_sample_k;
        cfg_k.compute_flowfield    = is_sample_k;
        % Live sub-case progress (HPC_Linux addition): reuse the same
        % DataQueue used for completion messages -- run_simulation_HPC_Linux
        % sends throttled 'heartbeat' messages through it (15 s interval),
        % which dispatch_progress_msg below routes separately from the
        % 'done' completion messages.
        cfg_k.progress_queue       = progress_q;                             %#ok<PFBNS>
        cfg_k.progress_tag         = k;
        run_seed(k) = seeds(k);

        dir_k = '';
        t_run = tic;
        try
            [res_k, ~, dir_k] = run_simulation_HPC_Linux(cfg_k);
            wall_k = toc(t_run);
            run_wall(k) = wall_k;

            if isempty(res_k)
                run_status(k) = "FAILED";
                run_errmsg(k) = "run_simulation_HPC_Linux returned empty results";
                run_errid(k)  = "N/A (empty results)";
                run_errloc(k) = "";
                run_outdir(k) = string(dir_k);
                send(progress_q, struct('type', 'done', 'k', k, 'status', "FAILED", ...
                    'msg', char(run_errmsg(k)), 'wall', wall_k, ...
                    'seed', seeds(k), 'outdir', char(string(dir_k)), ...
                    'identifier', char(run_errid(k)), ...
                    'stack_name', '', 'stack_line', NaN));
            else
                run_status(k) = "OK";
                run_outdir(k) = string(dir_k);
                send(progress_q, struct('type', 'done', 'k', k, 'status', "OK", ...
                    'msg', '', 'wall', wall_k, ...
                    'seed', seeds(k), 'outdir', char(string(dir_k)), ...
                    'identifier', '', 'stack_name', '', 'stack_line', NaN));
            end
        catch ME
            wall_k = toc(t_run);
            run_wall(k)   = wall_k;
            run_status(k) = "FAILED";
            run_errmsg(k) = string(ME.message);
            run_errid(k)  = string(ME.identifier);
            expected_dir_k = fullfile(out_root, sprintf('%03d', k));
            dir_k_capture = capture_partial_outdir(dir_k, expected_dir_k);
            run_outdir(k) = string(dir_k_capture);
            if ~isempty(ME.stack)
                run_errloc(k) = sprintf('%s (line %d)', ME.stack(1).name, ME.stack(1).line);
                stack_name_k  = ME.stack(1).name;
                stack_line_k  = ME.stack(1).line;
            else
                run_errloc(k) = "";
                stack_name_k  = '';
                stack_line_k  = NaN;
            end
            send(progress_q, struct('type', 'done', 'k', k, 'status', "FAILED", ...
                'msg', ME.message, 'wall', wall_k, ...
                'seed', seeds(k), 'outdir', dir_k_capture, ...
                'identifier', ME.identifier, ...
                'stack_name', stack_name_k, 'stack_line', stack_line_k));
        end
    end

else
    %% ---- SERIAL BRANCH (legacy for-loop; summary lines written incrementally) ----
    for k = 1:n_repeats
        sub = sprintf('%03d', k);
        fprintf('%s\n  Repeat %d/%d  ->  output/%s/%s\n%s\n', ...
            repmat('-', 1, 56), k, n_repeats, case_name_base, sub, repmat('-', 1, 56));

        cfg_k = cfg;
        cfg_k.case_name = fullfile(case_name_base, sub);

        % Heavy-output sampling: force animation/flowfield OFF except for the
        % sampled repeats, to shrink peak memory (see USER SETTINGS).
        is_sample_k              = ismember(k, heavy_output_sample_idx);
        cfg_k.generate_animation = is_sample_k;
        cfg_k.compute_flowfield  = is_sample_k;

        if strcmp(brown_rng_seed_mode, 'fixed')
            cfg_k.brown_rng_seed = rng_seed_base + k;
            run_seed(k) = rng_seed_base + k;
        else
            % Timestamp policy: let run_simulation_HPC_Linux use brown_rng_seed_from_clock (MMddHHmmss).
            if isfield(cfg_k, 'brown_rng_seed')
                cfg_k = rmfield(cfg_k, 'brown_rng_seed');
            end
        end

        res_k   = [];
        dir_k   = '';
        t_run = tic;
        try
            [res_k, ~, dir_k] = run_simulation_HPC_Linux(cfg_k);
            run_wall(k) = toc(t_run);

            if isempty(res_k)
                run_status(k) = "FAILED";
                run_errmsg(k) = "run_simulation_HPC_Linux returned empty results";
                run_errid(k)  = "N/A (empty results)";
                run_errloc(k) = "";
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
            expected_dir_k = fullfile(out_root, sub);
            run_outdir(k) = string(capture_partial_outdir(dir_k, expected_dir_k));
            run_errid(k)  = string(ME.identifier);
            if ~isempty(ME.stack)
                run_errloc(k) = sprintf('%s (line %d)', ME.stack(1).name, ME.stack(1).line);
            else
                run_errloc(k) = "";
            end
            fprintf('  [ERR]  %s\n', ME.message);
            if ~isempty(ME.stack)
                fprintf('         Location: %s (line %d)\n', ...
                    ME.stack(1).name, ME.stack(1).line);
            end
        end

        % --- Incremental summary write (per-run, flushed immediately) ---
        append_run_summary(sum_path, k, run_seed(k), run_status(k), ...
            run_wall(k), run_errid(k), run_errmsg(k), run_errloc(k), ...
            run_outdir(k), brown_rng_seed_mode);
    end
end

catch ME_batch
    % Whole-batch abort (not caught by any per-iteration try/catch above).
    % Report it clearly, then continue to FINAL REPORT with whatever
    % partial run_status/run_wall/... data already made it back to the
    % client, instead of letting the script die with a bare stack trace
    % and no summary.
    batch_fatal_error = ME_batch;
    fprintf('\n%s\n', sep);
    fprintf('  [FATAL] Batch execution aborted: %s\n', ME_batch.message);
    if ~isempty(ME_batch.stack)
        fprintf('          Location: %s (line %d)\n', ...
            ME_batch.stack(1).name, ME_batch.stack(1).line);
    end
    fprintf('  Continuing to final report with partial results...\n');
    fprintf('%s\n', sep);
end

elapsed = toc(t_batch);

%% ========================================================================
%  FINAL REPORT
%  ========================================================================
%  Per-run summary records are written incrementally inside the loops
%  (serial) or via the dispatch_progress_msg -> on_run_complete DataQueue
%  callback chain (parallel) — see append_run_summary.  Only the footer
%  (aggregate stats) is written here, in append mode.

n_ok      = sum(run_status == "OK");
n_failed  = sum(run_status == "FAILED");
n_pending = sum(run_status == "PENDING");

fprintf('\n%s\n', sep);
if isempty(batch_fatal_error)
    fprintf('  Batch complete: %.1f s total\n', elapsed);
else
    fprintf('  Batch ABORTED after %.1f s (see [FATAL] above): %s\n', ...
        elapsed, batch_fatal_error.message);
end
fprintf('  Success: %d / %d\n', n_ok, n_repeats);
if n_failed > 0
    fprintf('  Failed:  %d\n', n_failed);
end
if n_pending > 0
    fprintf('  Pending: %d (not completed)\n', n_pending);
end
if use_parallel
    sum_wall = sum(run_wall(~isnan(run_wall)));
    speedup  = sum_wall / max(elapsed, eps);
    fprintf('  Parallel scheduling report:\n');
    fprintf('    n_pool                  : %d\n',          parallel_info.n_pool);
    fprintf('    physical cores detected : %d\n',          parallel_info.n_physical_cores);
    fprintf('    available memory (start): %.2f GB\n',     parallel_info.available_mem_gb);
    fprintf('    cap                     : cpu %.0f%% -> %d, mem %.0f%% -> %d (binding=%s)\n', ...
        cpu_resource_frac*100, parallel_info.cpu_cap, ...
        mem_resource_frac*100, parallel_info.mem_cap, char(parallel_info.binding));
    fprintf('    sum of per-case wall    : %.1f s\n',      sum_wall);
    fprintf('    parallel speed-up est.  : %.2fx\n',       speedup);
end
fprintf('%s\n', sep);

if ~isempty(sum_path)
    fid_foot = fopen(sum_path, 'a', 'n', 'UTF-8');
    if fid_foot ~= -1
        finished_human = char(datetime('now', 'Format', 'uuuu-MM-dd HH:mm:ss'));
        fprintf(fid_foot, '\nfinished: %s\n', finished_human);
        fprintf(fid_foot, 'elapsed_total_s: %.2f\n', elapsed);
        if isempty(batch_fatal_error)
            fprintf(fid_foot, 'batch_aborted: false\n');
        else
            fprintf(fid_foot, 'batch_aborted: true\n');
            fprintf(fid_foot, 'batch_abort_identifier: %s\n', batch_fatal_error.identifier);
            fprintf(fid_foot, 'batch_abort_message: %s\n', ...
                regexprep(batch_fatal_error.message, '[\r\n]+', ' '));
        end
        fprintf(fid_foot, 'success_count: %d / %d\n', n_ok, n_repeats);
        if use_parallel
            sum_wall = sum(run_wall(~isnan(run_wall)));
            fprintf(fid_foot, 'sum_per_case_wall_s: %.2f\n',       sum_wall);
            fprintf(fid_foot, 'parallel_speedup_estimate: %.2f\n', sum_wall / max(elapsed, eps));
        end
        fprintf(fid_foot, 'failed_count: %d\n', n_failed);
        if n_failed > 0
            failed_idx = find(run_status == "FAILED");
            idx_labels = arrayfun(@(k) sprintf('%03d', k), failed_idx, 'UniformOutput', false);
            fprintf(fid_foot, 'failed_indices: %s\n', strjoin(idx_labels, ', '));
            fprintf(fid_foot, 'failed_runs:\n');
            for ii = 1:n_failed
                k = failed_idx(ii);
                if run_seed(k) > 0
                    seed_str = sprintf('%d', run_seed(k));
                else
                    seed_str = 'n/a';
                end
                errid_str = char(run_errid(k));
                if strlength(errid_str) == 0
                    errid_str = 'n/a';
                end
                msg_str = regexprep(char(run_errmsg(k)), '[\r\n]+', ' ');
                if numel(msg_str) > 120
                    msg_str = [msg_str(1:117), '...'];
                end
                if strlength(run_outdir(k)) > 0
                    outdir_str = path_for_log(char(run_outdir(k)));
                else
                    outdir_str = 'n/a';
                end
                fprintf(fid_foot, '  k=%03d  seed=%s  id=%s  msg=%s  outdir=%s\n', ...
                    k, seed_str, errid_str, msg_str, outdir_str);
            end
        end
        fprintf(fid_foot, 'pending_count: %d\n', n_pending);
        if n_pending > 0
            pending_idx = find(run_status == "PENDING");
            idx_labels = arrayfun(@(k) sprintf('%03d', k), pending_idx, 'UniformOutput', false);
            fprintf(fid_foot, 'pending_indices: %s\n', strjoin(idx_labels, ', '));
        end
        fclose(fid_foot);
        fprintf('\n  Summary written: %s\n', sum_path);
    end
end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function apply_max_num_comp_threads(n)
%APPLY_MAX_NUM_COMP_THREADS  Worker-side wrapper for parfevalOnAll BLAS cap.
    maxNumCompThreads(n);
end


function [n_pool, info] = decide_pool_size(n_repeats, mem_per_worker_gb, cpu_frac, mem_frac, info, blas_threads)
%DECIDE_POOL_SIZE  Probe machine resources and compute capped worker count.
%
%  HPC_Linux variant: adds a Linux-native physical-core probe (lscpu) and a
%  Linux-native available-memory probe (/proc/meminfo MemAvailable), and
%  fixes a bug in the final feature('numCores') fallback -- feature('numCores')
%  already returns the PHYSICAL (not logical) core count (confirmed against
%  MATLAB community references and MathWorks' own "local pool default size"
%  documentation), so the previous "/2" halving under-counted real capacity
%  by ~2x.  This fallback is reached on EVERY Linux run (the preceding
%  wmic/CIM/NUMBER_OF_PROCESSORS probes are all Windows-only), so the fix
%  matters a great deal there even though it rarely triggered on Windows.
%
%  cpu_frac and mem_frac are now separate (previously a single max_frac
%  capped both) -- see USER SETTINGS in the caller for why memory is kept
%  more conservative by default.
%
%  blas_threads is accepted for API stability and future diagnostic use but
%  is not currently read inside this function (the console probe block was
%  removed; the main batch header prints BLAS threads directly).  Kept as a
%  parameter so callers do not need to change if diagnostics are restored.
%
%  Returns the smaller of:
%    - floor(n_physical_cores * cpu_frac)
%    - floor(available_mem_gb * mem_frac / mem_per_worker_gb)
%    - n_repeats
%  and labels which constraint binds (used for the diagnostic header).

    if nargin < 6 || isempty(blas_threads)
        blas_threads = 1;   %#ok<NASGU>  retained for API stability; not read below
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

    if (~isfinite(n_phys) || n_phys < 1) && isunix
        % Linux-native probe: count unique (Core,Socket) pairs from lscpu, the
        % standard portable way to get the PHYSICAL core count (excludes
        % hyper-threaded sibling duplicates).  lscpu (util-linux) is present
        % on essentially every Linux distribution, including HPC images.
        try
            [st_lscpu, out_lscpu] = system('lscpu -p=Core,Socket 2>/dev/null | grep -v "^#"');
            if st_lscpu == 0
                lines = strtrim(strsplit(strtrim(out_lscpu), newline));
                lines = lines(~cellfun(@isempty, lines));
                if ~isempty(lines)
                    n_phys = numel(unique(lines));
                    probe_note = "lscpu";
                end
            end
        catch
        end
    end

    if ~isfinite(n_phys) || n_phys < 1
        try
            % feature('numCores') already returns the PHYSICAL core count on
            % both Windows and Linux (verified against MATLAB community
            % references) -- do NOT halve it here.
            n_phys = max(1, round(feature('numCores')));
            probe_note = "feature(numCores)";
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

    %% --- Available physical memory ---
    avail_gb = NaN;
    mem_probe_note = "";
    try
        if ispc
            [~, sys] = memory;
            avail_gb = sys.PhysicalMemory.Available / 1e9;
            mem_probe_note = "memory()";
        end
    catch
        avail_gb = NaN;
    end

    if (~isfinite(avail_gb) || avail_gb <= 0) && isunix
        % Linux-native probe: memory() is Windows-only (confirmed against
        % MATLAB's own documentation) and throws on Linux.  MemAvailable in
        % /proc/meminfo is the kernel's own "usable without swapping"
        % estimate (more accurate than MemFree, standard since Linux 3.14),
        % present on every modern Linux including HPC compute-node images.
        try
            meminfo_txt = fileread('/proc/meminfo');
            tok = regexp(meminfo_txt, 'MemAvailable:\s*(\d+)\s*kB', 'tokens', 'once');
            if ~isempty(tok)
                avail_gb = str2double(tok{1}) * 1024 / 1e9;   % kB -> GB
                mem_probe_note = "/proc/meminfo:MemAvailable";
            end
        catch
            avail_gb = NaN;
        end
    end

    if ~isfinite(avail_gb) || avail_gb <= 0
        avail_gb = 4.0;
        mem_probe_note = "fallback=4GB";
        fprintf(['  [WARN] Memory probe unavailable on this platform; assuming 4 GB ', ...
            'available.  This will badly under-size the pool on a real HPC node -- ', ...
            'verify /proc/meminfo is readable, or set mem_per_worker_gb/mem_resource_frac ', ...
            'manually if this warning appears on your cluster.\n']);
    end
    info.available_mem_gb = avail_gb;
    info.mem_probe = mem_probe_note;

    %% --- Two independent caps + n_repeats ---
    n_cpu = max(1, floor(n_phys   * cpu_frac));
    n_mem = max(1, floor(avail_gb * mem_frac / mem_per_worker_gb));
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

    % Console probe output removed: the main batch header (printed after this
    % function returns) already shows physical cores, memory, caps and binding.
    % The probe_note source is recorded in the summary file (core_probe /
    % mem_probe fields).
end


function pool = ensure_parpool(n_target)
%ENSURE_PARPOOL  Reuse-or-create a process-based parpool of the requested size.
%
%  Tries 'Processes' (R2022a+) first, falls back to 'local' on older MATLAB.
%  Returns [] if pool creation fails.
%
%  HPC_Linux variant: on Linux, first tries relocating the local cluster's
%  JobStorageLocation (the small inter-process coordination/lock files
%  parpool uses internally) away from its MATLAB default -- which resolves
%  under the user's home directory and, on HPC clusters where the home
%  directory lives on distributed/network storage (Lustre/NFS/BeeGFS/...),
%  is a documented cause of slow pool startup and, when concurrent jobs
%  share that same home-directory path, race conditions or outright
%  pool-start failures (see MathWorks Answers #437544; Harvard FASRC and
%  other HPC-center docs recommend exactly this relocation).  The relocated
%  path uses MATLAB's own tempdir() (which already honours $TMPDIR --
%  typically a node-local scratch path set by the job scheduler -- falling
%  back to /tmp otherwise), in a job-unique subfolder so concurrent jobs on
%  the same or different nodes never collide.  If this specific attempt
%  fails for any reason, falls straight back to the original
%  'Processes'/'local' string-profile calls (MATLAB default
%  JobStorageLocation), so this change can never leave pool creation worse
%  off than before.  Windows behaviour is completely unchanged.

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

    if isunix
        try
            job_dir = local_job_storage_location();
            if ~isfolder(job_dir)
                mkdir(job_dir);
            end
            c = parcluster('local');
            c.NumWorkers = n_target;
            c.JobStorageLocation = job_dir;
            pool = c.parpool(n_target);
            fprintf('  [INFO] parpool JobStorageLocation (node-local): %s\n', job_dir);
            return;
        catch ME
            fprintf(['  [WARN] Relocated-JobStorageLocation parpool start failed (%s); ', ...
                'falling back to default JobStorageLocation.\n'], ME.message);
            pool = [];
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


function dir_path = local_job_storage_location()
%LOCAL_JOB_STORAGE_LOCATION  Node-local, job-unique scratch dir for the
%  local parallel pool's JobStorageLocation on Linux (see ensure_parpool).
%  Prefers a scheduler job ID (SLURM_JOB_ID, then PBS_JOBID) for
%  traceability; falls back to a random UUID if neither is set (e.g. an
%  interactive, non-scheduled Linux session).
    job_id = strtrim(getenv('SLURM_JOB_ID'));
    if isempty(job_id)
        job_id = strtrim(getenv('PBS_JOBID'));
    end
    if isempty(job_id)
        job_id = char(java.util.UUID.randomUUID);
    end
    user_tag = strtrim(getenv('USER'));
    if isempty(user_tag)
        user_tag = strtrim(getenv('USERNAME'));
    end
    if isempty(user_tag)
        user_tag = 'user';
    end
    dir_path = fullfile(tempdir(), sprintf('matlab_pool_%s_%s', user_tag, job_id));
end


function append_run_summary(sum_path, k, seed, status, wall, errid, errmsg, errloc, outdir, seed_mode)
%APPEND_RUN_SUMMARY  Append one per-run record to summary file with forced flush.
%
%  Opens the summary file in append mode, writes a single per-run record
%  (seed, status, wall time, and on failure the identifier/message/location,
%  outdir, and log_hint when present, or on success the output directory),
%  then immediately closes the file.
%  The close-then-reopen pattern guarantees the record reaches disk before
%  the next run starts, so a mid-batch crash never loses already-completed
%  records.  No-op when sum_path is empty (summary disabled).
%
%  Inputs:
%    sum_path  - path to summary file (string/char); '' disables writing
%    k         - 1-based run index
%    seed      - Brownian RNG seed used (>0); 0 means unknown/unset
%    status    - "OK" | "FAILED" | "PENDING" (string)
%    wall      - wall time in seconds (NaN if not run)
%    errid     - ME.identifier on failure (string); "" on success
%    errmsg    - ME.message on failure (string); "" on success
%    errloc    - 'name (line N)' on failure (string); "" on success
%    outdir    - output directory of the run (string); "" on failure
%    seed_mode - 'fixed' | 'timestamp' (controls the '(timestamp)' tag)

    if isempty(sum_path)
        return;
    end
    fid = fopen(sum_path, 'a', 'n', 'UTF-8');
    if fid == -1
        warning('batch_repeat_brownian_HPC_Linux:SummaryAppendFailed', ...
            'Could not open summary file for append:\n  %s', sum_path);
        return;
    end

    if seed > 0
        seed_tag = sprintf('seed=%d', seed);
        if ~strcmp(seed_mode, 'fixed')
            seed_tag = [seed_tag, ' (timestamp)']; %#ok<AGROW>
        end
    else
        seed_tag = 'seed=n/a';
    end

    outdir_char = char(string(outdir));

    fprintf(fid, '%03d  %s  status=%s  wall=%.2f s\n', k, seed_tag, status, wall);
    if status == "FAILED" && strlength(errmsg) > 0
        fprintf(fid, '         Identifier : %s\n', errid);
        fprintf(fid, '         Message    : %s\n', errmsg);
        if strlength(errloc) > 0
            fprintf(fid, '         Location   : %s\n', errloc);
        end
    elseif status == "OK" && strlength(outdir_char) > 0
        fprintf(fid, '  out: %s\n', path_for_log(outdir_char));
    end
    if status == "FAILED" && strlength(outdir_char) > 0
        fprintf(fid, '         outdir     : %s\n', path_for_log(outdir_char));
        log_hint = fullfile(outdir_char, 'simulation.log');
        if isfile(log_hint)
            fprintf(fid, '         log_hint   : %s\n', path_for_log(log_hint));
        end
    end

    fclose(fid);   % forced flush — record is on disk before next run starts
end


function s = path_for_log(p)
%PATH_FOR_LOG  Prefer project-relative paths in summary files for portability.
    s = char(string(p));
    pw = char(pwd);
    sep_ch = filesep;
    if ~isempty(pw) && pw(end) ~= sep_ch
        pw = [pw, sep_ch]; %#ok<AGROW>
    end
    if strlength(s) >= strlength(pw) && strcmpi(s(1:strlength(pw)), pw)
        s = s(strlength(pw)+1:end);
    end
end


function s = brown_rng_policy_line(seed_mode, use_par, rng_base, n_rep, seeds_vec)
%BROWN_RNG_POLICY_LINE  One-line summary of the Brownian seed rule in effect.
    if strcmp(seed_mode, 'fixed')
        s = sprintf('seed(k) = rng_seed_base + k, rng_seed_base=%d, k=1..%d', ...
            rng_base, n_rep);
    elseif use_par && ~isempty(seeds_vec)
        clock_base = seeds_vec(1) - 1;
        s = sprintf('parallel: seed(k) = clock_base + k, clock_base=%d, k=1..%d', ...
            clock_base, n_rep);
    else
        s = ['serial: MMddHHmmss per run via brown_rng_seed_from_clock() ', ...
             'at run_simulation_HPC_Linux start'];
    end
end


function write_resource_probe_summary(fid, info, mem_per_worker_gb, cpu_frac, mem_frac, blas_threads)
%WRITE_RESOURCE_PROBE_SUMMARY  Resource cap fields for summary header/footer.
    fprintf(fid, 'n_physical_cores: %d\n',     info.n_physical_cores);
    if isfield(info, 'core_probe')
        fprintf(fid, 'core_probe: %s\n',       char(info.core_probe));
    end
    fprintf(fid, 'available_mem_gb: %.2f\n',   info.available_mem_gb);
    if isfield(info, 'mem_probe')
        fprintf(fid, 'mem_probe: %s\n',        char(info.mem_probe));
    end
    fprintf(fid, 'mem_per_worker_gb: %g\n',    mem_per_worker_gb);
    fprintf(fid, 'cpu_resource_frac: %g\n',    cpu_frac);
    fprintf(fid, 'mem_resource_frac: %g\n',    mem_frac);
    fprintf(fid, 'cpu_cap: %d\n',              info.cpu_cap);
    fprintf(fid, 'mem_cap: %d\n',              info.mem_cap);
    fprintf(fid, 'binding: %s\n',              char(info.binding));
    fprintf(fid, 'worker_blas_threads: %d\n',  blas_threads);
end


function outdir_char = capture_partial_outdir(dir_k, expected_dir)
%CAPTURE_PARTIAL_OUTDIR  Best-effort recovery of a failed run's output folder.
%
%  MATLAB does not propagate a function's output arguments back to the
%  caller when that function throws an uncaught error -- even though
%  run_simulation_HPC_Linux computes output_dir (and creates the folder, and may
%  have written a partial simulation.log with a [FATAL] block) well
%  before the point where the error is raised, [res_k, ~, dir_k] =
%  run_simulation_HPC_Linux(cfg_k) leaves dir_k at whatever value it held before
%  the call (verified empirically: a function erroring after assigning
%  an output argument still returns nothing to the caller).  So dir_k
%  here is always '' on failure and cannot recover the real path.
%
%  Instead, since cfg_k.case_name is deterministic per repeat
%  (case_name_base/NNN) and run_simulation_HPC_Linux only appends a timestamp
%  suffix when the target folder already contains prior run data (which
%  a fresh repeat folder normally does not), the caller passes the
%  expected path as expected_dir; we confirm it actually exists on disk
%  before reporting it, so a batch_repeat_brownian_HPC_Linux bug or a genuinely
%  never-created folder still correctly yields ''.
    outdir_char = '';
    if nargin >= 1 && ~isempty(dir_k)
        dir_str = strtrim(char(string(dir_k)));
        if ~isempty(dir_str)
            outdir_char = dir_str;
            return;
        end
    end
    if nargin >= 2 && ~isempty(expected_dir) && isfolder(expected_dir)
        outdir_char = expected_dir;
    end
end


function dispatch_progress_msg(msg, n_total, sum_path, seed_mode)
%DISPATCH_PROGRESS_MSG  DataQueue callback: routes 'heartbeat' vs 'done' messages.
%
%  run_simulation_HPC_Linux sends 'heartbeat' messages (intra-run progress,
%  at most once per 15 s per repeat) through the same DataQueue used for
%  the existing per-repeat 'done' completion messages.  Both share one
%  queue/one afterEach callback (simplest way to guarantee they interleave
%  correctly and never race on the client thread), so this dispatcher is
%  the single entry point that tells them apart via msg.type and routes
%  to the appropriate handler.  Older-style messages with no 'type' field
%  are treated as 'done' for backward compatibility (should not occur in
%  practice since every send() site now sets 'type' explicitly).
    if isfield(msg, 'type') && strcmp(msg.type, 'heartbeat')
        print_heartbeat(msg);
    else
        on_run_complete(msg, n_total, sum_path, seed_mode);
    end
end


function print_heartbeat(msg)
%PRINT_HEARTBEAT  Console-only intra-run progress line for a still-running repeat.
%
%  Deliberately NOT written to the persistent summary file (append_run_summary)
%  -- the summary's job is to record final per-repeat outcomes, not the
%  dozens of intermediate heartbeats a single long repeat may emit; keeping
%  them console-only avoids bloating that file with noise.  Printed as a
%  simple appended line (no cursor/redraw tricks): HPC batch output is
%  typically redirected to a log file, not a live terminal, so append-only
%  is the robust choice here.
    fprintf('  [HB] k=%03d running  %5.1f%%  elapsed=%s\n', ...
        msg.tag, msg.pct, format_duration_hms(msg.elapsed));
end


function s = format_duration_hms(seconds_in)
%FORMAT_DURATION_HMS  Compact H:MM:SS (or M:SS) formatting for heartbeat elapsed time.
    seconds_in = max(0, seconds_in);
    h = floor(seconds_in / 3600);
    m = floor(mod(seconds_in, 3600) / 60);
    s_val = mod(seconds_in, 60);
    if h > 0
        s = sprintf('%d:%02d:%02.0f', h, m, s_val);
    else
        s = sprintf('%d:%02.0f', m, s_val);
    end
end


function on_run_complete(msg, n_total, sum_path, seed_mode)
%ON_RUN_COMPLETE  DataQueue callback: print console progress + append summary.
%
%  Called by dispatch_progress_msg on the main thread whenever a parfor
%  worker sends a completion notice.  Delegates the console line (with ETA)
%  to print_progress and the persistent summary append to
%  append_run_summary.  Runs on the main thread, so there is no
%  concurrent-file-write race.

    print_progress(msg, n_total);

    errloc = '';
    if isfield(msg, 'stack_name') && ~isempty(msg.stack_name) ...
            && isfield(msg, 'stack_line') && isfinite(msg.stack_line)
        errloc = sprintf('%s (line %d)', msg.stack_name, msg.stack_line);
    end
    outdir = '';
    if isfield(msg, 'outdir')
        outdir = char(msg.outdir);
    end
    seed = 0;
    if isfield(msg, 'seed')
        seed = msg.seed;
    end

    append_run_summary(sum_path, msg.k, seed, msg.status, msg.wall, ...
        msg.identifier, msg.msg, errloc, outdir, seed_mode);
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
        % Extended error details (identifier + location), mirroring the
        % serial branch's console output and batch_run_from_json's parallel
        % branch.  Fields are always present (sent by both catch paths).
        if isfield(msg, 'identifier') && ~isempty(msg.identifier)
            fprintf('         Identifier : %s\n', msg.identifier);
        end
        if isfield(msg, 'stack_name') && ~isempty(msg.stack_name) ...
                && isfield(msg, 'stack_line') && isfinite(msg.stack_line)
            fprintf('         Location   : %s (line %d)\n', ...
                msg.stack_name, msg.stack_line);
        end
    end
end


