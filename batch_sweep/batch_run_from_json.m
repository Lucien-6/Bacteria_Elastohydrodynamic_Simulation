%% batch_run_from_json.m
%  Automated Batch Runner: Execute All JSON Config Files in a Folder
%
%  Description:
%    Interactively prompts the user to select a folder, then discovers
%    all *.json configuration files within it.  Each file is loaded via
%    load_sim_config(), forced into headless mode, and executed via
%    run_simulation().  If a run fails for any reason (bad config,
%    numerical divergence, MATLAB error), the error is caught and logged,
%    and execution continues with the next config file.  A full summary
%    table is printed upon completion.
%
%    Workflow:
%      1. User selects the target folder via a GUI dialog (or supplies
%         BATCH_JSON_DIR at the top of this file to bypass the dialog).
%      2. All *.json files are sorted alphabetically and queued.
%      3. Pre-scan: resolve case_name; skip duplicate non-empty case_name
%         (first file in sort order wins; later duplicates -> SKIPPED).
%      4. For each queued file (parallel or serial):
%           a. Load config with load_sim_config()
%           b. Override auto_confirm = true  (suppress all dialogs)
%              Override preview_initial = false  (skip pose preview)
%           c. Call run_simulation(cfg)  inside try-catch
%           d. On success: record V_sink and output directory
%           e. On failure: record error message and continue
%      5. Print a formatted summary table of all runs.
%
%    Key points:
%      - Skips failed runs; never halts the entire batch
%      - Headless flags are always enforced regardless of JSON settings
%      - Brownian RNG, animation, flowfield follow each JSON (no overrides)
%      - Optional parallel mode (use_parallel, default true) via parfor
%      - Resource cap: max_resource_frac = 0.75 (reserve 25% CPU and RAM)
%      - mem_per_worker_gb = 3.0 (typical per-case memory budget)
%      - Empty case_name files are never deduplicated (run_simulation assigns
%        unique run-<timestamp> output directories)
%      - Each run writes its own output/ subdirectory with full logs
%      - Real-time batch log: batch_summary_<timestamp>.log in JSON folder
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-18
%  Updated: 2026-05-25 — Lucien — optional parfor parallel branch, case_name
%                        deduplication, real-time DataQueue logging, unified
%                        run header/itemized summary (P0+P2 blog consistency).

close all force;
clear; clc;

%% ---- Project root (robust to current directory) ----
this_file = mfilename('fullpath');
if isempty(this_file)
    proj_root = pwd;
else
    proj_root = fileparts(fileparts(this_file));
end
if ~isfolder(proj_root) || ~isfile(fullfile(proj_root, 'run_simulation.m'))
    error('batch_run_from_json:BadRoot', ...
        ['Cannot resolve project root (run_simulation.m not found). ', ...
         'Run from Coupled_Body_Pili/batch_sweep/ or cd to project root.']);
end
cd(proj_root);
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));

%% ---- Optional: hard-code folder path to skip the GUI dialog ----
%  Leave empty ('') to use the interactive folder browser.
BATCH_JSON_DIR = '';

%% ---- USER SETTINGS (parallel scheduling) ----
%  use_parallel        : true  -> parfor branch (auto-sized pool, default)
%                        false -> serial for-loop branch
%  mem_per_worker_gb   : conservative per-worker memory budget [GB]
%  max_resource_frac   : cap on BOTH physical-core count and available memory
%                        (0.75 reserves 25% of CPU and RAM for the host)
%  worker_blas_threads : per-worker maxNumCompThreads (1 avoids thread storm)
use_parallel        = true;
mem_per_worker_gb   = 3.0;
max_resource_frac   = 0.75;
worker_blas_threads = 1;

%% ---- Folder selection ----
if isempty(BATCH_JSON_DIR)
    fprintf('Select the folder containing JSON configuration files...\n');
    selected = uigetdir(pwd, 'Select JSON Config Folder');
    if isequal(selected, 0)
        fprintf('No folder selected. Batch run cancelled.\n');
        return;
    end
    json_folder = selected;
else
    json_folder = BATCH_JSON_DIR;
end

if ~isfolder(json_folder)
    error('batch_run_from_json:InvalidFolder', ...
        'Target folder does not exist:\n  %s', json_folder);
end

%% ---- Discover JSON files ----
file_list = dir(fullfile(json_folder, '*.json'));
file_list = file_list(~[file_list.isdir]);

if isempty(file_list)
    fprintf('No *.json files found in:\n  %s\nBatch run cancelled.\n', json_folder);
    return;
end

% Sort alphabetically by filename
[~, sort_idx] = sort({file_list.name});
file_list     = file_list(sort_idx);
n_total       = numel(file_list);

%% ---- Open batch summary log ----
log_timestamp  = char(datetime('now', 'Format', 'uuuuMMdd''_''HHmmss'));
batch_log_path = fullfile(json_folder, ...
    sprintf('batch_summary_%s.log', log_timestamp));
blog_fid = fopen(batch_log_path, 'w', 'n', 'UTF-8');
if blog_fid == -1
    warning('batch_run_from_json:LogFailed', ...
        'Cannot create batch log file. Continuing without log.');
end

blog = @(varargin) batch_log(blog_fid, varargin{:});

%% ---- Print header ----
sep = repmat('=', 1, 56);
fprintf('\n%s\n', sep);
fprintf('  Batch Runner: JSON Config Files\n');
fprintf('  Folder : %s\n', json_folder);
fprintf('  Found  : %d config file(s)\n', n_total);
fprintf('  Start  : %s\n', log_timestamp);
fprintf('%s\n\n', sep);

blog('%s\n', sep);
blog('  Batch Runner: JSON Config Files\n');
blog('  Folder : %s\n', json_folder);
blog('  Found  : %d config file(s)\n', n_total);
blog('  Start  : %s\n', log_timestamp);
blog('%s\n\n', sep);

for ii = 1:n_total
    fprintf('  [%02d/%02d]  %s\n', ii, n_total, file_list(ii).name);
    blog('  [%02d/%02d]  %s\n', ii, n_total, file_list(ii).name);
end
fprintf('\n');
blog('\n');

%% ---- Pre-allocate result tracking ----
run_status   = repmat("PENDING",  n_total, 1);  % OK | FAILED | SKIPPED | PENDING
run_Vsink    = nan(n_total, 1);                 % Steady V_sink [um/s]
run_walltime = nan(n_total, 1);                 % Wall time [s]
run_outdir   = repmat("", n_total, 1);          % Output directory
run_errmsg   = repmat("", n_total, 1);          % Error / skip reason

%% ---- Pre-scan and build run queue ----
seen_case          = containers.Map('KeyType', 'char', 'ValueType', 'char');
run_ii             = zeros(0, 1);
run_json_path      = cell(0, 1);
n_skipped_prescan  = 0;
n_prescan_fail     = 0;

fprintf('%s\n', repmat('-', 1, 56));
fprintf('  Pre-scan: case_name deduplication\n');
fprintf('%s\n', repmat('-', 1, 56));
blog('%s\n', repmat('-', 1, 56));
blog('  Pre-scan: case_name deduplication\n');
blog('%s\n', repmat('-', 1, 56));

for ii = 1:n_total
    json_path = fullfile(json_folder, file_list(ii).name);
    fname     = file_list(ii).name;

    try
        cfg_probe = load_sim_config(json_path);
        cn        = strtrim(cfg_probe.case_name);

        if ~isempty(cn) && isKey(seen_case, cn)
            run_status(ii)  = "SKIPPED";
            run_errmsg(ii)  = sprintf( ...
                'Duplicate case_name "%s"; first file: %s', cn, seen_case(cn));
            n_skipped_prescan = n_skipped_prescan + 1;
            fprintf('  [SKIPPED] [%02d/%02d] %s  |  %s\n', ...
                ii, n_total, fname, run_errmsg(ii));
            blog('  [SKIPPED] [%02d/%02d] %s  |  %s\n', ...
                ii, n_total, fname, run_errmsg(ii));
            batch_log_flush(blog_fid);
        else
            if ~isempty(cn)
                seen_case(cn) = fname;
            end
            run_ii(end+1, 1)        = ii;           %#ok<AGROW>
            run_json_path{end+1, 1} = json_path;     %#ok<AGROW>
        end

    catch ME
        run_status(ii)   = "FAILED";
        run_errmsg(ii)   = string(ME.message);
        n_prescan_fail   = n_prescan_fail + 1;
        fprintf('  [PRESCAN FAIL] [%02d/%02d] %s  |  %s\n', ...
            ii, n_total, fname, ME.message);
        blog('  [PRESCAN FAIL] [%02d/%02d] %s  |  %s\n', ...
            ii, n_total, fname, ME.message);
        batch_log_flush(blog_fid);
    end
end

n_run = numel(run_ii);
fprintf('\n  Total JSON: %d  |  To run: %d  |  Skipped (duplicate): %d  |  Prescan failed: %d\n\n', ...
    n_total, n_run, n_skipped_prescan, n_prescan_fail);
blog('\n  Total JSON: %d  |  To run: %d  |  Skipped (duplicate): %d  |  Prescan failed: %d\n\n', ...
    n_total, n_run, n_skipped_prescan, n_prescan_fail);
batch_log_flush(blog_fid);

%% ---- Parallel scheduling decision ----
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

pool_handle = [];

if n_run > 0
    if use_parallel
        if ~license('test', 'Distrib_Computing_Toolbox')
            warning('batch_run_from_json:NoPCT', ...
                'Parallel Computing Toolbox is not licensed. Falling back to serial.');
            use_parallel = false;
            parallel_info.fallback_reason = "PCT license unavailable";
        end
    end

    if use_parallel
        [n_pool_decided, parallel_info] = decide_pool_size( ...
            n_run, mem_per_worker_gb, max_resource_frac, parallel_info, ...
            worker_blas_threads);
        parallel_info.n_pool = n_pool_decided;
    end

    if use_parallel
        pool_handle = ensure_parpool(parallel_info.n_pool);
        if isempty(pool_handle)
            warning('batch_run_from_json:PoolCreateFailed', ...
                'parpool creation failed. Falling back to serial.');
            use_parallel = false;
            parallel_info.fallback_reason = "parpool creation failed";
        else
            parallel_info.n_pool = pool_handle.NumWorkers;
            f_vis = parfevalOnAll(pool_handle, @() set(0, 'DefaultFigureVisible', 'off'), 0);
            f_thr = parfevalOnAll(pool_handle, @apply_max_num_comp_threads, 0, worker_blas_threads);
            wait(f_vis);
            wait(f_thr);
            % Worker defaults (figure visibility, BLAS threads) persist for the
            % lifetime of this parpool — intentional for batch use.  To restore
            % MATLAB defaults after the batch, run:  delete(gcp('nocreate'))
        end
    end
end

%% ---- Execute batch ----
t_batch_start = tic;

if n_run > 0
    exec_sep = repmat('-', 1, 56);
    fprintf('%s\n', exec_sep);
    blog('%s\n', exec_sep);
    if use_parallel && ~isempty(pool_handle)
        fprintf('  Execution mode: PARALLEL  (n_pool = %d, %s-bound)\n', ...
            parallel_info.n_pool, char(parallel_info.binding));
        blog('  Execution mode: PARALLEL  (n_pool = %d, %s-bound)\n', ...
            parallel_info.n_pool, char(parallel_info.binding));
    else
        fprintf('  Execution mode: SERIAL\n');
        blog('  Execution mode: SERIAL\n');
        if strlength(parallel_info.fallback_reason) > 0
            fprintf('  Parallel fallback: %s\n', char(parallel_info.fallback_reason));
            blog('  Parallel fallback: %s\n', char(parallel_info.fallback_reason));
        end
    end
    fprintf('%s\n\n', exec_sep);
    blog('%s\n\n', exec_sep);
    batch_log_flush(blog_fid);

    if use_parallel && ~isempty(pool_handle)
        %% ---- PARALLEL BRANCH (parfor + DataQueue progress) ----
        run_fname = cell(n_run, 1);
        for k = 1:n_run
            run_fname{k} = file_list(run_ii(k)).name;
        end

        pq_status   = repmat("PENDING", n_run, 1);
        pq_Vsink    = nan(n_run, 1);
        pq_walltime = nan(n_run, 1);
        pq_outdir   = repmat("", n_run, 1);
        pq_errmsg   = repmat("", n_run, 1);

        progress_q = parallel.pool.DataQueue;
        afterEach(progress_q, @(msg) on_progress_impl(msg, n_run, n_total, blog, blog_fid));
        print_progress_json(struct('reset', true), n_run, n_total);

        parfor k = 1:n_run
            json_path = run_json_path{k};
            fname     = run_fname{k};

            cfg_k                    = load_sim_config(json_path);
            cfg_k.auto_confirm       = true;
            cfg_k.preview_initial    = false;

            t_run        = tic;
            status_k     = "FAILED";
            V_steady_k   = NaN;
            errmsg_k     = "";
            outdir_k     = "";
            identifier_k = "";
            stack_name_k = "";
            stack_line_k = NaN;

            try
                [res_i, ~, dir_i] = run_simulation(cfg_k);
                wall_k = toc(t_run);

                [status_k, is_complete_k] = simulation_result_status(res_i);
                if ~is_complete_k
                    if isempty(res_i)
                        errmsg_k = "run_simulation returned empty results";
                    else
                        errmsg_k = "run_simulation returned partial results";
                    end
                    outdir_k  = char(string(dir_i));
                elseif isfield(res_i, 'V_sink') && ~isempty(res_i.V_sink)
                    Vs            = res_i.V_sink;
                    V_steady_k    = mean(Vs(max(1, round(0.95 * numel(Vs))):end));
                    outdir_k      = char(string(dir_i));
                else
                    status_k = "FAILED";
                    errmsg_k = "results missing V_sink field";
                    outdir_k = char(string(dir_i));
                end

                pq_walltime(k) = wall_k;
                pq_status(k)   = status_k;
                pq_Vsink(k)    = V_steady_k;
                pq_outdir(k)   = string(outdir_k);
                pq_errmsg(k)   = string(errmsg_k);

                send(progress_q, struct( ...
                    'ii', run_ii(k), 'k', k, 'fname', fname, 'status', status_k, ...
                    'V_steady', V_steady_k, 'wall', wall_k, 'errmsg', char(errmsg_k), ...
                    'outdir', outdir_k, 'identifier', char(identifier_k), ...
                    'stack_name', char(stack_name_k), 'stack_line', stack_line_k));

            catch ME
                wall_k = toc(t_run);
                pq_walltime(k) = wall_k;
                pq_status(k)   = "FAILED";
                pq_errmsg(k)   = string(ME.message);
                pq_outdir(k)   = string(capture_partial_outdir('', cfg_k.case_name));

                identifier_k = ME.identifier;
                if ~isempty(ME.stack)
                    stack_name_k = ME.stack(1).name;
                    stack_line_k = ME.stack(1).line;
                end

                send(progress_q, struct( ...
                    'ii', run_ii(k), 'k', k, 'fname', fname, 'status', "FAILED", ...
                    'V_steady', NaN, 'wall', wall_k, 'errmsg', ME.message, ...
                    'outdir', char(pq_outdir(k)), 'identifier', identifier_k, ...
                    'stack_name', stack_name_k, 'stack_line', stack_line_k));
            end
        end

        for k = 1:n_run
            ii = run_ii(k);
            run_status(ii)   = pq_status(k);
            run_Vsink(ii)    = pq_Vsink(k);
            run_walltime(ii) = pq_walltime(k);
            run_outdir(ii)   = pq_outdir(k);
            run_errmsg(ii)   = pq_errmsg(k);
        end

    else
        %% ---- SERIAL BRANCH ----
        n_serial_done = 0;
        for k = 1:n_run
            ii        = run_ii(k);
            json_path = run_json_path{k};

            run_line = format_run_queue_line(k, n_run, ii, n_total, file_list(ii).name);
            fprintf('\n%s\n', repmat('-', 1, 56));
            fprintf('%s\n', run_line);
            fprintf('%s\n', repmat('-', 1, 56));
            blog('\n%s\n', repmat('-', 1, 56));
            blog('%s\n', run_line);
            blog('%s\n', repmat('-', 1, 56));
            batch_log_flush(blog_fid);

            t_run_start = tic;

            try
                cfg                    = load_sim_config(json_path);
                cfg.auto_confirm       = true;
                cfg.preview_initial    = false;

                [res_i, ~, dir_i] = run_simulation(cfg);
                run_walltime(ii) = toc(t_run_start);

                [run_status(ii), is_complete_i] = simulation_result_status(res_i);
                if ~is_complete_i
                    if isempty(res_i)
                        run_errmsg(ii) = "run_simulation returned empty results";
                    else
                        run_errmsg(ii) = "run_simulation returned partial results";
                    end
                    run_outdir(ii)  = string(dir_i);
                    fprintf('  [%s] Run %d returned partial/empty results.\n', run_status(ii), ii);
                    blog('  [%s] Run %d returned partial/empty results.\n', run_status(ii), ii);
                elseif isfield(res_i, 'V_sink') && ~isempty(res_i.V_sink)
                    Vs              = res_i.V_sink;
                    V_steady        = mean(Vs(max(1, round(0.95 * numel(Vs))):end));
                    run_Vsink(ii)   = V_steady;
                    run_outdir(ii)  = string(dir_i);
                    fprintf('  [OK]   Run %d complete  |  V_sink = %.4e um/s  |  %.1f s\n', ...
                        ii, V_steady, run_walltime(ii));
                    blog('  [OK]   Run %d complete  |  V_sink = %.4e um/s  |  %.1f s\n', ...
                        ii, V_steady, run_walltime(ii));
                else
                    run_status(ii)  = "FAILED";
                    run_errmsg(ii)  = "results missing V_sink field";
                    run_outdir(ii)  = string(dir_i);
                    fprintf('  [FAIL] Run %d: results missing V_sink field.\n', ii);
                    blog('  [FAIL] Run %d: results missing V_sink field.\n', ii);
                end

            catch ME
                run_walltime(ii) = toc(t_run_start);
                run_status(ii)   = "FAILED";
                run_errmsg(ii)   = string(ME.message);
                if exist('cfg', 'var') && isstruct(cfg) && isfield(cfg, 'case_name')
                    run_outdir(ii) = string(capture_partial_outdir('', cfg.case_name));
                end

                fprintf('  [FAIL] Run %d FAILED after %.1f s\n', ii, run_walltime(ii));
                fprintf('         Error: %s\n', ME.message);
                blog('  [FAIL] Run %d FAILED after %.1f s\n', ii, run_walltime(ii));
                blog('         Identifier : %s\n', ME.identifier);
                blog('         Message    : %s\n', ME.message);
                if ~isempty(ME.stack)
                    blog('         Location   : %s (line %d)\n', ...
                        ME.stack(1).name, ME.stack(1).line);
                end
                if strlength(run_outdir(ii)) > 0
                    outdir_hint = char(run_outdir(ii));
                    fprintf('         Outdir: %s\n', outdir_hint);
                    blog('         Outdir     : %s\n', outdir_hint);
                    log_hint = fullfile(outdir_hint, 'simulation.log');
                    if isfile(log_hint)
                        fprintf('         Log hint: %s\n', log_hint);
                        blog('         Log hint   : %s\n', log_hint);
                    end
                end
            end

            batch_log_flush(blog_fid);

            n_serial_done = n_serial_done + 1;
            if n_serial_done > 0 && n_serial_done < n_run
                done_mask = run_status(run_ii(1:k)) ~= "PENDING";
                avg_t     = mean(run_walltime(run_ii(done_mask)), 'omitnan');
                eta_s     = avg_t * (n_run - n_serial_done);
                fprintf('  ETA remaining: %.0f s  (~%.1f min)\n', eta_s, eta_s / 60);
            end
        end
    end
end

t_batch_total = toc(t_batch_start);

%% ---- Final Summary ----
n_ok      = sum(run_status == "OK");
n_failed  = sum(run_status == "FAILED");
n_partial = sum(run_status == "TIMEOUT" | run_status == "ABORTED");
n_skipped = sum(run_status == "SKIPPED");

fprintf('\n%s\n', sep);
fprintf('  BATCH COMPLETE\n');
fprintf('  Total JSON files : %d\n', n_total);
fprintf('  Queued to run    : %d\n', n_run);
fprintf('  Succeeded        : %d\n', n_ok);
fprintf('  Failed           : %d\n', n_failed);
fprintf('  Partial          : %d (excluded from steady statistics)\n', n_partial);
fprintf('  Skipped          : %d\n', n_skipped);
fprintf('  Total wall time  : %.1f s  (%.1f min)\n', t_batch_total, t_batch_total / 60);
fprintf('%s\n', sep);

blog('\n%s\n', sep);
blog('  BATCH COMPLETE\n');
blog('  Total JSON files : %d\n', n_total);
blog('  Queued to run    : %d\n', n_run);
blog('  Succeeded        : %d\n', n_ok);
blog('  Failed           : %d\n', n_failed);
blog('  Partial          : %d (excluded from steady statistics)\n', n_partial);
blog('  Skipped          : %d\n', n_skipped);
blog('  Total wall time  : %.1f s  (%.1f min)\n', t_batch_total, t_batch_total / 60);
blog('%s\n', sep);

if n_run > 0 && use_parallel && ~isempty(pool_handle)
    sum_wall = sum(run_walltime(run_ii), 'omitnan');
    speedup  = sum_wall / max(t_batch_total, eps);
    fprintf('  sum_per_case_wall_s         : %.2f\n', sum_wall);
    fprintf('  parallel_speedup_estimate   : %.2f\n', speedup);
    blog('  sum_per_case_wall_s         : %.2f\n', sum_wall);
    blog('  parallel_speedup_estimate   : %.2f\n', speedup);
end

% --- Per-run result table ---
fprintf('\n  %-4s  %-32s  %-8s  %14s  %10s\n', ...
    '#', 'Config File', 'Status', 'V_sink (um/s)', 'Time (s)');
fprintf('  %s\n', repmat('-', 1, 74));
blog('\n  %-4s  %-32s  %-8s  %14s  %10s\n', ...
    '#', 'Config File', 'Status', 'V_sink (um/s)', 'Time (s)');
blog('  %s\n', repmat('-', 1, 74));

for ii = 1:n_total
    fname_short = file_list(ii).name;
    if length(fname_short) > 32
        fname_short = ['...' fname_short(end-28:end)];
    end

    if run_status(ii) == "OK"
        fprintf('  %-4d  %-32s  %-8s  %14.4e  %10.1f\n', ...
            ii, fname_short, 'OK', run_Vsink(ii), run_walltime(ii));
        blog('  %-4d  %-32s  %-8s  %14.4e  %10.1f\n', ...
            ii, fname_short, 'OK', run_Vsink(ii), run_walltime(ii));
    elseif run_status(ii) == "FAILED"
        fprintf('  %-4d  %-32s  %-8s  %14s  %10.1f\n', ...
            ii, fname_short, 'FAILED', '---', run_walltime(ii));
        blog('  %-4d  %-32s  %-8s  %14s  %10.1f\n', ...
            ii, fname_short, 'FAILED', '---', run_walltime(ii));
        fprintf('         Error: %s\n', run_errmsg(ii));
        blog('         Error: %s\n', run_errmsg(ii));
        if strlength(run_outdir(ii)) > 0
            outdir_hint = char(run_outdir(ii));
            fprintf('         Outdir: %s\n', outdir_hint);
            blog('         Outdir: %s\n', outdir_hint);
            log_hint = fullfile(outdir_hint, 'simulation.log');
            if isfile(log_hint)
                fprintf('         Log hint: %s\n', log_hint);
                blog('         Log hint: %s\n', log_hint);
            end
        end
    elseif run_status(ii) == "SKIPPED"
        fprintf('  %-4d  %-32s  %-8s  %14s  %10s\n', ...
            ii, fname_short, 'SKIPPED', '---', '---');
        blog('  %-4d  %-32s  %-8s  %14s  %10s\n', ...
            ii, fname_short, 'SKIPPED', '---', '---');
        fprintf('         Reason: %s\n', run_errmsg(ii));
        blog('         Reason: %s\n', run_errmsg(ii));
    else
        fprintf('  %-4d  %-32s  %-8s\n', ii, fname_short, char(run_status(ii)));
        blog('  %-4d  %-32s  %-8s\n', ii, fname_short, char(run_status(ii)));
    end
end

fprintf('  %s\n', repmat('-', 1, 74));
blog('  %s\n', repmat('-', 1, 74));

% --- Log file location ---
if blog_fid ~= -1
    fclose(blog_fid);
    fprintf('\n  Batch log saved: %s\n', batch_log_path);
end
fprintf('%s\n\n', sep);


%% ==================== Local Functions ====================

function batch_log(fid, varargin)
%BATCH_LOG  Write formatted message to batch log file (if open)
    if fid == -1
        return;
    end
    fprintf(fid, varargin{:});
end


function batch_log_flush(fid)
%BATCH_LOG_FLUSH  Flush batch log file handle for real-time logging
    if fid ~= -1
        fflush(fid);
    end
end


function on_progress_impl(msg, n_run, n_total, blog_fn, blog_fid)
%ON_PROGRESS_IMPL  DataQueue callback: console progress + real-time batch log
    print_progress_json(msg, n_run, n_total);

    blog_fn('\n%s\n', repmat('-', 1, 56));
    blog_fn('%s\n', format_run_queue_line(msg.k, n_run, msg.ii, n_total, msg.fname));
    blog_fn('%s\n', repmat('-', 1, 56));

    if msg.status == "OK"
        blog_fn('  [OK]   JSON %02d complete  |  V_sink = %.4e um/s  |  %.1f s\n', ...
            msg.ii, msg.V_steady, msg.wall);
    else
        blog_fn('  [FAIL] JSON %02d FAILED after %.1f s\n', msg.ii, msg.wall);
        if isfield(msg, 'identifier') && ~isempty(msg.identifier)
            blog_fn('         Identifier : %s\n', msg.identifier);
        end
        if isfield(msg, 'errmsg') && ~isempty(msg.errmsg)
            blog_fn('         Message    : %s\n', msg.errmsg);
        end
        if isfield(msg, 'stack_name') && ~isempty(msg.stack_name) ...
                && isfield(msg, 'stack_line') && isfinite(msg.stack_line)
            blog_fn('         Location   : %s (line %d)\n', ...
                msg.stack_name, msg.stack_line);
        end
        if isfield(msg, 'outdir') && strlength(string(msg.outdir)) > 0
            outdir_hint = char(msg.outdir);
            blog_fn('         Outdir     : %s\n', outdir_hint);
            log_hint = fullfile(outdir_hint, 'simulation.log');
            if isfile(log_hint)
                blog_fn('         Log hint   : %s\n', log_hint);
            end
        end
    end

    batch_log_flush(blog_fid);
end


function line = format_run_queue_line(k, n_run, ii, n_total, fname)
%FORMAT_RUN_QUEUE_LINE  Unified run header: queue index + original JSON index
    line = sprintf('  Run queue %02d/%02d  [JSON %02d/%02d]: %s', ...
        k, n_run, ii, n_total, fname);
end


function apply_max_num_comp_threads(n)
%APPLY_MAX_NUM_COMP_THREADS  Worker-side wrapper for parfevalOnAll BLAS cap.
    maxNumCompThreads(n);
end


function outdir_char = capture_partial_outdir(dir_i, case_name)
%CAPTURE_PARTIAL_OUTDIR  Best-effort recovery of a failed run's output folder.
%
%  MATLAB does not propagate a function's output arguments back to the
%  caller when that function throws an uncaught error -- even though
%  run_simulation computes output_dir (and creates the folder, and may
%  have written a partial simulation.log with a [FATAL] block) well
%  before the point where the error is raised, [res_i, ~, dir_i] =
%  run_simulation(cfg) leaves dir_i at whatever value it held before the
%  call.  So dir_i is always '' on failure and cannot recover the real
%  path directly.
%
%  When case_name is non-empty, run_simulation derives output_dir
%  deterministically as output/<case_name> (only appending a timestamp
%  suffix if that folder already held prior run data), so we reconstruct
%  the expected path and confirm it exists on disk.  Empty case_name is
%  not recoverable this way, because run_simulation then generates a
%  wall-clock-based run-<timestamp> folder name that cannot be predicted
%  in advance.
    outdir_char = '';
    if nargin >= 1 && ~isempty(dir_i)
        dir_str = strtrim(char(string(dir_i)));
        if ~isempty(dir_str)
            outdir_char = dir_str;
            return;
        end
    end
    if nargin >= 2 && ~isempty(case_name)
        cn = strtrim(char(string(case_name)));
        if ~isempty(cn)
            expected_dir = fullfile(pwd, 'output', cn);
            if isfolder(expected_dir)
                outdir_char = expected_dir;
            end
        end
    end
end


function [n_pool, info] = decide_pool_size(n_cases, mem_per_worker_gb, max_frac, info, blas_threads)
%DECIDE_POOL_SIZE  Probe machine resources and compute capped worker count.

    if nargin < 5 || isempty(blas_threads)
        blas_threads = 1;
    end

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

    n_cpu = max(1, floor(n_phys   * max_frac));
    n_mem = max(1, floor(avail_gb * max_frac / mem_per_worker_gb));
    info.cpu_cap = n_cpu;
    info.mem_cap = n_mem;

    n_pool = min([n_cpu, n_mem, n_cases]);

    if n_pool == n_cases && n_cases <= n_cpu && n_cases <= n_mem
        info.binding = "n_cases";
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
    fprintf('         n_cases        : %d\n', n_cases);
    fprintf('         BLAS threads   : %d / worker (after parfevalOnAll init)\n', blas_threads);
    fprintf('         => n_pool = %d  (binding: %s)\n\n', n_pool, char(info.binding));
end


function pool = ensure_parpool(n_target)
%ENSURE_PARPOOL  Reuse-or-create a process-based parpool of the requested size.

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


function print_progress_json(msg, n_run, n_total)
%PRINT_PROGRESS_JSON  Console streaming progress for JSON batch parfor runs.

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
        return;
    end

    n_done  = n_done + 1;
    elapsed = toc(t0);
    eta_s   = elapsed * max(0, n_run - n_done) / max(n_done, 1);
    run_tag = format_run_queue_tag(msg.k, n_run, msg.ii, n_total, msg.fname);

    if msg.status == "OK"
        fprintf('  [%3d/%3d] OK     %s  wall=%6.1f s  | elapsed=%6.1f s  ETA=%6.1f s\n', ...
            n_done, n_run, run_tag, msg.wall, elapsed, eta_s);
    else
        err_one_line = regexprep(char(msg.errmsg), '[\r\n]+', ' ');
        if numel(err_one_line) > 80
            err_one_line = [err_one_line(1:77), '...'];
        end
        fprintf('  [%3d/%3d] FAILED %s  wall=%6.1f s  | err: %s\n', ...
            n_done, n_run, run_tag, msg.wall, err_one_line);
    end
end


function tag = format_run_queue_tag(k, n_run, ii, n_total, fname)
%FORMAT_RUN_QUEUE_TAG  Compact run tag for inline console progress lines
    tag = sprintf('Run queue %02d/%02d  [JSON %02d/%02d]: %s', ...
        k, n_run, ii, n_total, fname);
end
