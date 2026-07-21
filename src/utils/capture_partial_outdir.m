function outdir_char = capture_partial_outdir(dir_k, expected_dir)
%CAPTURE_PARTIAL_OUTDIR  Best-effort recovery of a failed run's output folder.
%
%  Syntax:
%    outdir_char = capture_partial_outdir(dir_k, expected_dir)
%
%  MATLAB does not propagate a function's output arguments back to the
%  caller when that function throws an uncaught error -- even though
%  run_simulation computes output_dir (and creates the folder, and may
%  have written a partial simulation.log with a [FATAL] block) well
%  before the point where the error is raised, [res_k, ~, dir_k] =
%  run_simulation(cfg_k) leaves dir_k at whatever value it held before
%  the call (verified empirically: a function erroring after assigning
%  an output argument still returns nothing to the caller).  So dir_k
%  here is always '' on failure and cannot recover the real path.
%
%  Instead, since cfg_k.case_name is deterministic per repeat
%  (case_name_base/NNN) and run_simulation only appends a timestamp
%  suffix when the target folder already contains prior run data (which
%  a fresh repeat folder normally does not), the caller passes the
%  expected path as expected_dir; we confirm it actually exists on disk
%  before reporting it, so a batch driver bug or a genuinely
%  never-created folder still correctly yields ''.
%
%  This helper lives under src/utils (not as a script-local function)
%  so parfor workers can resolve it via the MATLAB path.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-03
%  Updated: 2026-07-16 — Lucien — moved from batch_repeat_brownian.m
%                        local function so parfor can resolve it.

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
