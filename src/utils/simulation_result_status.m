function [status, is_complete] = simulation_result_status(results)
%SIMULATION_RESULT_STATUS  Classify complete and partial simulation returns
%
%  Syntax:
%    [status, is_complete] = simulation_result_status(results)
%
%  Returns "OK" only for a completed, non-empty result.  A timeout is a
%  normal partial return from run_simulation and is classified as "TIMEOUT"
%  so batch drivers do not accidentally include it in steady-state metrics.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

if isempty(results)
    status = "ABORTED";
    is_complete = false;
    return;
end

if isfield(results, 'completed') && ~logical(results.completed)
    is_complete = false;
    reason = "aborted";
    if isfield(results, 'termination_reason') && ...
            ~isempty(results.termination_reason)
        reason = lower(string(results.termination_reason));
    end
    if reason == "timeout"
        status = "TIMEOUT";
    else
        status = "ABORTED";
    end
    return;
end

status = "OK";
is_complete = true;

end
