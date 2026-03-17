function render_progress_bar(tag, current, total, tstart)
%RENDER_PROGRESS_BAR  Generic in-place text progress bar for long loops
%
%  Description:
%    Displays a real-time text progress bar in the Command Window, similar
%    to odetpbar but for arbitrary loops (animation rendering, flow field
%    evaluation, etc.).  Updates are throttled to avoid excessive I/O.
%    Uses persistent variables for cross-call state (same pattern as odetpbar).
%
%  Syntax:
%    render_progress_bar(tag, current, total, tstart)
%
%  Inputs:
%    tag     - Prefix label string, e.g. '[ANIM]  ' or '[FLOW-XZ] '
%    current - Current step number (1-based)
%    total   - Total number of steps
%    tstart  - tic handle for elapsed / ETA computation
%
%  Display format:
%    [ANIM]   [==============................]  46.7%  Elapsed: 12.3 s  ETA: 14.0 s
%
%  Note:
%    This function only writes to the console (fprintf).  It does NOT write
%    to any log file.  Call logf() separately for a final completion line.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17
%  Updated: 2026-03-17 — v3.2: Fixed pass-by-value bug; replaced prev_len_container
%                         parameter with persistent prev_msg_len (same pattern as odetpbar)

persistent prev_msg_len last_upd_time;

BAR_W    = 30;
THROTTLE = 0.25;   % Min seconds between display refreshes

if current == 1
    prev_msg_len  = 0;
    last_upd_time = -Inf;
end

now_s = toc(tstart);

if (now_s - last_upd_time) < THROTTLE && current < total
    return;
end

pct    = 100 * current / total;
filled = round(pct / 100 * BAR_W);
bar_str = [repmat('=', 1, filled), repmat('.', 1, BAR_W - filled)];

if pct > 0.5 && now_s > 0.5
    eta_s = now_s * (total - current) / current;
    msg = sprintf('%s[%s] %5.1f%%  Elapsed: %s  ETA: %s', ...
        tag, bar_str, pct, fmt_dur_rpb(now_s), fmt_dur_rpb(eta_s));
else
    msg = sprintf('%s[%s] %5.1f%%  Elapsed: %s', ...
        tag, bar_str, pct, fmt_dur_rpb(now_s));
end

% In-place overwrite via \b (backspace) — same mechanism as odetpbar.
% State is tracked via persistent prev_msg_len (NOT via a parameter,
% because MATLAB passes numeric arrays by value, not by reference).
if prev_msg_len > 0
    fprintf(repmat('\b', 1, prev_msg_len));
end
fprintf('%s', msg);
prev_msg_len  = numel(msg);
last_upd_time = now_s;

if current == total
    fprintf('\n');
    prev_msg_len = 0;
end

drawnow;

end


function s = fmt_dur_rpb(t)
%FMT_DUR_RPB  Convert seconds to human-readable duration string
    if t < 0
        s = '-- s';
    elseif t < 60
        s = sprintf('%.1f s', t);
    elseif t < 3600
        s = sprintf('%d min %.0f s', floor(t/60), mod(t, 60));
    else
        s = sprintf('%d h %d min', floor(t/3600), floor(mod(t, 3600)/60));
    end
end
