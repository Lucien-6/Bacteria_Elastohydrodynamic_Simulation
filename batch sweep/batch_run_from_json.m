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
%      3. For each file:
%           a. Load config with load_sim_config()
%           b. Override auto_confirm = true  (suppress all dialogs)
%              Override preview_initial = false  (skip pose preview)
%           c. Call run_simulation(cfg)  inside try-catch
%           d. On success: record V_sink and output directory
%           e. On failure: record error message and continue
%      4. Print a formatted summary table of all runs.
%
%    Key points:
%      - Skips failed runs; never halts the entire batch
%      - Headless flags are always enforced regardless of JSON settings
%      - Each run writes its own output/ subdirectory with full logs
%      - The batch itself produces a batch_summary.log in the JSON folder
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-18

close all force;
clear; clc;

%% ---- Optional: hard-code folder path to skip the GUI dialog ----
%  Leave empty ('') to use the interactive folder browser.
BATCH_JSON_DIR = '';

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
run_status   = repmat("PENDING",  n_total, 1);  % "OK" | "FAILED" | "PENDING"
run_Vsink    = nan(n_total, 1);                 % Steady V_sink [um/s]
run_walltime = nan(n_total, 1);                 % Wall time [s]
run_outdir   = repmat("", n_total, 1);          % Output directory
run_errmsg   = repmat("", n_total, 1);          % Error message (if failed)

%% ---- Execute batch ----
t_batch_start = tic;

for ii = 1:n_total
    json_path = fullfile(json_folder, file_list(ii).name);

    fprintf('\n%s\n', repmat('-', 1, 56));
    fprintf('  Run %d/%d: %s\n', ii, n_total, file_list(ii).name);
    fprintf('%s\n', repmat('-', 1, 56));
    blog('\n%s\n', repmat('-', 1, 56));
    blog('  Run %d/%d: %s\n', ii, n_total, file_list(ii).name);
    blog('%s\n', repmat('-', 1, 56));

    t_run_start = tic;

    try
        % --- Load and patch config (enforce headless mode) ---
        cfg                    = load_sim_config(json_path);
        cfg.auto_confirm    = true;    % Suppress all interactive dialogs
        cfg.preview_initial = false;   % Skip initial pose preview
        % generate_animation: intentionally left as loaded from JSON

        % --- Execute simulation ---
        [res_i, ~, dir_i] = run_simulation(cfg);

        run_walltime(ii) = toc(t_run_start);

        if isempty(res_i)
            % run_simulation returned gracefully but produced no results
            % (e.g. user-abort path — should not happen in headless mode)
            run_status(ii)  = "FAILED";
            run_errmsg(ii)  = "run_simulation returned empty results (possible abort)";
            run_outdir(ii)  = string(dir_i);
            fprintf('  [WARN] Run %d returned empty results.\n', ii);
            blog('  [WARN] Run %d returned empty results.\n', ii);
        else
            Vs              = res_i.V_sink;
            V_steady        = mean(Vs(max(1, round(0.95 * numel(Vs))):end));
            run_status(ii)  = "OK";
            run_Vsink(ii)   = V_steady;
            run_outdir(ii)  = string(dir_i);
            fprintf('  [OK]   Run %d complete  |  V_sink = %.4e um/s  |  %.1f s\n', ...
                ii, V_steady, run_walltime(ii));
            blog('  [OK]   Run %d complete  |  V_sink = %.4e um/s  |  %.1f s\n', ...
                ii, V_steady, run_walltime(ii));
        end

    catch ME
        run_walltime(ii) = toc(t_run_start);
        run_status(ii)   = "FAILED";
        run_errmsg(ii)   = string(ME.message);

        fprintf('  [FAIL] Run %d FAILED after %.1f s\n', ii, run_walltime(ii));
        fprintf('         Error: %s\n', ME.message);
        blog('  [FAIL] Run %d FAILED after %.1f s\n', ii, run_walltime(ii));
        blog('         Identifier : %s\n', ME.identifier);
        blog('         Message    : %s\n', ME.message);
        if ~isempty(ME.stack)
            blog('         Location   : %s (line %d)\n', ...
                ME.stack(1).name, ME.stack(1).line);
        end
    end

    % ETA estimate based on average so far
    completed = sum(run_status ~= "PENDING");
    if completed > 0 && completed < n_total
        avg_t    = mean(run_walltime(1:ii), 'omitnan');
        eta_s    = avg_t * (n_total - completed);
        fprintf('  ETA remaining: %.0f s  (~%.1f min)\n', eta_s, eta_s / 60);
    end
end

t_batch_total = toc(t_batch_start);

%% ---- Final Summary ----
n_ok     = sum(run_status == "OK");
n_failed = sum(run_status == "FAILED");

fprintf('\n%s\n', sep);
fprintf('  BATCH COMPLETE\n');
fprintf('  Total runs  : %d\n', n_total);
fprintf('  Succeeded   : %d\n', n_ok);
fprintf('  Failed      : %d\n', n_failed);
fprintf('  Total time  : %.1f s  (%.1f min)\n', t_batch_total, t_batch_total / 60);
fprintf('%s\n', sep);

blog('\n%s\n', sep);
blog('  BATCH COMPLETE\n');
blog('  Total runs  : %d\n', n_total);
blog('  Succeeded   : %d\n', n_ok);
blog('  Failed      : %d\n', n_failed);
blog('  Total time  : %.1f s  (%.1f min)\n', t_batch_total, t_batch_total / 60);
blog('%s\n', sep);

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
    else
        fprintf('  %-4d  %-32s  %-8s\n', ii, fname_short, 'PENDING');
        blog('  %-4d  %-32s  %-8s\n', ii, fname_short, 'PENDING');
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


%% ==================== Local Function ====================

function batch_log(fid, varargin)
%BATCH_LOG  Write formatted message to batch log file (if open)
    if fid == -1
        return;
    end
    fprintf(fid, varargin{:});
end
