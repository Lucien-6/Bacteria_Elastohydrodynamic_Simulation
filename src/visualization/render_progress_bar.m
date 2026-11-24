function render_progress_bar(tag, current, total, tstart, detail, eta_opts)
%RENDER_PROGRESS_BAR  Generic adaptive in-place text progress bar for long loops
%
%  Description:
%    Displays a real-time text progress bar in the Command Window for
%    arbitrary loops (animation rendering, flow-field evaluation, etc.).
%    Uses an internal ProgressEstimator (cost-rate EWMA) by default and
%    supports optional heterogeneous / multi-class ETA overrides via
%    eta_opts.  Automatically degrades to throttled newline-delimited
%    output in headless / redirected sessions (~usejava('desktop')).
%
%  Syntax:
%    render_progress_bar(tag, current, total, tstart)
%    render_progress_bar(tag, current, total, tstart, detail)
%    render_progress_bar(tag, current, total, tstart, detail, eta_opts)
%
%  Inputs:
%    tag      - Prefix label string, e.g. '[ANIM]   ' or '[FLOW] Compute '
%    current  - Current step number (1-based)
%    total    - Total number of steps
%    tstart   - tic handle for elapsed / ETA computation
%    detail   - (Optional) short status inserted after the tag, e.g. '3/54'
%    eta_opts - (Optional) scalar ETA [s] override, or struct with optional
%               fields .eta_s (override ETA [s]) and .pct (override % in [0,100])
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17
%  Updated: 2026-11-24

persistent prev_msg_len last_upd_time est prev_tag is_desktop_cached;

BAR_W           = 30;
THROTTLE_TTY    = 0.25;   % Min seconds between desktop in-place refreshes
THROTTLE_BATCH  = 15.0;   % Min seconds between non-desktop appended lines

if isempty(is_desktop_cached)
    is_desktop_cached = usejava('desktop');
end
is_tty   = is_desktop_cached;
throttle = THROTTLE_TTY;
if ~is_tty
    throttle = THROTTLE_BATCH;
end

if nargin < 5 || isempty(detail)
    detail = '';
elseif ~endsWith(detail, ' ')
    detail = [detail, ' '];
end

total_safe = max(double(total), 1.0);
cur_safe   = min(max(double(current), 0.0), total_safe);

is_new_loop = (cur_safe <= 1) || isempty(est) || isempty(prev_tag) ...
    || ~strcmp(prev_tag, char(tag)) || (est.total_work ~= total_safe);
if is_new_loop
    prev_msg_len  = 0;
    last_upd_time = -Inf;
    prev_tag      = char(tag);
    est = ProgressEstimator(total_safe, 'TicStart', tstart, ...
        'TauMin', 5.0, 'TauMax', 120.0, 'Gamma', 0.25, ...
        'MinDtSample', 0.10, 'WarmupSec', 0.5, 'WarmupPct', 0.5);
end

now_s = toc(tstart);
est.update(cur_safe, now_s);

if (now_s - last_upd_time) < throttle && cur_safe < total_safe
    return;
end

info = est.estimate(now_s);
pct   = 100.0 * cur_safe / total_safe;
eta_s = info.eta_s;
fin_s = info.finish_str;

if nargin >= 6 && ~isempty(eta_opts)
    if isnumeric(eta_opts) && isscalar(eta_opts)
        eta_s = double(eta_opts);
        fin_s = ProgressEstimator.formatFinishTime(eta_s);
    elseif isstruct(eta_opts)
        if isfield(eta_opts, 'pct') && ~isempty(eta_opts.pct) && isfinite(eta_opts.pct)
            pct = min(100.0, max(0.0, double(eta_opts.pct)));
        end
        if isfield(eta_opts, 'eta_s') && ~isempty(eta_opts.eta_s)
            eta_s = double(eta_opts.eta_s);
            fin_s = ProgressEstimator.formatFinishTime(eta_s);
        end
    end
end

if cur_safe >= total_safe
    pct   = 100.0;
    eta_s = 0.0;
end

filled  = min(BAR_W, max(0, round(pct / 100.0 * BAR_W)));
bar_str = [repmat('=', 1, filled), repmat('.', 1, BAR_W - filled)];
elap_s  = format_duration(now_s, 'compact');

if cur_safe >= total_safe
    msg = sprintf('%s%s[%s] 100.0%%  Elapsed: %s', ...
        tag, detail, bar_str, elap_s);
elseif isfinite(eta_s) && (pct > 0.5 || now_s > 0.5)
    eta_txt = format_duration(eta_s, 'compact');
    if eta_s >= 60 && ~isempty(fin_s)
        msg = sprintf('%s%s[%s] %5.1f%%  Elapsed: %s  ETA: %s (~%s)', ...
            tag, detail, bar_str, pct, elap_s, eta_txt, fin_s);
    else
        msg = sprintf('%s%s[%s] %5.1f%%  Elapsed: %s  ETA: %s', ...
            tag, detail, bar_str, pct, elap_s, eta_txt);
    end
else
    msg = sprintf('%s%s[%s] %5.1f%%  Elapsed: %s', ...
        tag, detail, bar_str, pct, elap_s);
end

if is_tty
    if prev_msg_len > 0
        fprintf(repmat('\b', 1, prev_msg_len));
    end
    fprintf('%s', msg);
    extra = prev_msg_len - numel(msg);
    if extra > 0
        fprintf('%s', repmat(' ', 1, extra));
        if cur_safe < total_safe
            fprintf(repmat('\b', 1, extra));
        end
    end
    prev_msg_len  = numel(msg);
    last_upd_time = now_s;

    if cur_safe >= total_safe
        fprintf('\n');
        prev_msg_len = 0;
        prev_tag     = '';
    end
else
    fprintf('%s\n', msg);
    prev_msg_len  = 0;
    last_upd_time = now_s;
    if cur_safe >= total_safe
        prev_tag = '';
    end
end

drawnow limitrate;

end
