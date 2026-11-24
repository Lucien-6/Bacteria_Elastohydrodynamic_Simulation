function status = odetpbar(t, ~, flag, ctx)
%ODETPBAR  Adaptive ODE integration progress bar (OutputFcn for ode15s)
%
%  Description:
%    Displays a real-time text progress bar during ode15s integration with
%    percentage, elapsed time, and adaptive cost-rate EWMA ETA (plus
%    estimated wall-clock finish time for runs longer than 1 min).
%
%    Supports both:
%      - Legacy 3-argument OutputFcn call: odetpbar(t, y, flag)
%        Correctly anchors progress at segment start t(1) (fixing multi-
%        segment t0 != 0 distortion) and uses an internal ProgressEstimator.
%      - Global-context 4-argument call: odetpbar(t, y, flag, ctx)
%        Displays true global progress across multi-segment integrations
%        (Brownian intervals, Euler-basis rotations) using a shared
%        ProgressEstimator handle passed in ctx.estimator.
%
%    In headless / redirected sessions (~usejava('desktop')), automatically
%    switches from in-place '\b' redraws to throttled newline-delimited
%    status lines so log files never accumulate backspace garbage.
%
%  Syntax:
%    status = odetpbar(t, y, flag)
%    status = odetpbar(t, y, flag, ctx)
%
%  Inputs:
%    t    - Scalar or vector of current time values (provided by ode15s)
%    y    - Current state vector (unused; required by OutputFcn interface)
%    flag - Integration stage string: '' (step), 'init', or 'done'
%    ctx  - (Optional) struct with fields:
%             .estimator      ProgressEstimator handle covering [0, T_total]
%             .T_total        Global target dimensionless time
%             .throttle_init  Logical: defer initial draw for short segments
%                             (default false)
%
%  Outputs:
%    status - Always 0 (continue integration)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-11-24

status = 0;
persistent local_est prev_len last_upd last_t_cur is_desktop_cached;

BAR_W           = 30;     % Progress bar width (characters)
THROTTLE_TTY    = 0.25;   % Min seconds between in-place desktop refreshes
THROTTLE_BATCH  = 15.0;   % Min seconds between appended lines in non-desktop mode

if isempty(is_desktop_cached)
    is_desktop_cached = usejava('desktop');
end
if isempty(prev_len)
    prev_len = 0;
end
is_tty   = is_desktop_cached;
throttle = THROTTLE_TTY;
if ~is_tty
    throttle = THROTTLE_BATCH;
end

if nargin < 4
    ctx = [];
end
has_ctx = isstruct(ctx) && isfield(ctx, 'estimator') && ~isempty(ctx.estimator);
throttle_init = has_ctx && isfield(ctx, 'throttle_init') && logical(ctx.throttle_init);

if isempty(flag)
    % ===== Normal integration step =====
    if isempty(t)
        return;
    end
    t_cur = t(end);
    last_t_cur = t_cur;

    est = active_estimator(has_ctx, ctx, local_est);
    if isempty(est)
        return;
    end
    now_s = toc(est.tic_start);
    est.update(t_cur, now_s);

    if ~isempty(last_upd) && (now_s - last_upd) < throttle
        return;
    end

    info = est.estimate(now_s);
    msg  = build_msg(info, BAR_W, false);
    prev_len = emit_line(msg, prev_len, is_tty, false);
    last_upd = now_s;
    drawnow limitrate;

else
    switch flag
        case 'init'
            t0 = 0;
            tf = 1;
            if ~isempty(t)
                t0 = t(1);
                tf = t(end);
            end
            last_t_cur = t0;
            prev_len   = 0;

            if has_ctx
                est = ctx.estimator;
                now_s = toc(est.tic_start);
                est.update(t0, now_s);
            else
                local_est = ProgressEstimator(tf, 'StartWork', t0, ...
                    'TauMin', 10.0, 'TauMax', 300.0, 'Gamma', 0.25, ...
                    'WarmupSec', 0.5, 'WarmupPct', 0.5);
                est   = local_est;
                now_s = 0.0;
            end

            if throttle_init
                % Short intermediate segment: do not draw unless it runs past throttle
                last_upd = now_s;
            else
                info = est.estimate(now_s);
                msg  = build_msg(info, BAR_W, false);
                prev_len = emit_line(msg, 0, is_tty, false);
                last_upd = now_s;
                drawnow limitrate;
            end

        case 'done'
            est = active_estimator(has_ctx, ctx, local_est);
            if ~isempty(est)
                now_s = toc(est.tic_start);
                if ~isempty(last_t_cur) && isfinite(last_t_cur)
                    est.update(last_t_cur, now_s);
                end
                info = est.estimate(now_s);
                is_complete = (info.pct >= 99.999);
                if prev_len > 0 || ~throttle_init || is_complete
                    msg = build_msg(info, BAR_W, is_complete);
                    emit_line(msg, prev_len, is_tty, true);
                    drawnow;
                end
            elseif prev_len > 0 && is_tty
                fprintf('\n');
                drawnow;
            end

            if ~has_ctx
                local_est = [];
            end
            prev_len   = 0;
            last_upd   = [];
            last_t_cur = [];

        otherwise
            error('odetpbar:UnknownFlag', 'Unknown flag: %s', flag);
    end
end

end


% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function est = active_estimator(has_ctx, ctx, local_est)
    if has_ctx
        est = ctx.estimator;
    else
        est = local_est;
    end
end


function new_len = emit_line(msg, prev_len, is_tty, end_with_newline)
    if is_tty
        if ~isempty(prev_len) && prev_len > 0
            fprintf(repmat('\b', 1, prev_len));
        end
        fprintf('%s', msg);
        extra = 0;
        if ~isempty(prev_len)
            extra = prev_len - numel(msg);
        end
        if extra > 0
            fprintf('%s', repmat(' ', 1, extra));
            if ~end_with_newline
                fprintf(repmat('\b', 1, extra));
            end
        end
        if end_with_newline
            fprintf('\n');
            new_len = 0;
        else
            new_len = numel(msg);
        end
    else
        fprintf('%s\n', msg);
        new_len = 0;
    end
end


function msg = build_msg(info, bar_w, is_complete)
%BUILD_MSG  Construct formatted progress bar string
    pct     = min(100.0, max(0.0, info.pct));
    filled  = min(bar_w, max(0, round(pct / 100.0 * bar_w)));
    bar_str = [repmat('=', 1, filled), repmat('.', 1, bar_w - filled)];
    elap_str = format_duration(info.elapsed_s, 'compact');

    if is_complete
        msg = sprintf('[ODE]   [%s] 100.0%%  Elapsed: %s', bar_str, elap_str);
    elseif isfield(info, 'eta_s') && isfinite(info.eta_s)
        eta_str = format_duration(info.eta_s, 'compact');
        if info.eta_s >= 60 && isfield(info, 'finish_str') && ~isempty(info.finish_str)
            msg = sprintf('[ODE]   [%s] %5.1f%%  Elapsed: %s  ETA: %s (~%s)', ...
                bar_str, pct, elap_str, eta_str, info.finish_str);
        else
            msg = sprintf('[ODE]   [%s] %5.1f%%  Elapsed: %s  ETA: %s', ...
                bar_str, pct, elap_str, eta_str);
        end
    else
        msg = sprintf('[ODE]   [%s] %5.1f%%  Elapsed: %s  ETA: --', ...
            bar_str, pct, elap_str);
    end
end
