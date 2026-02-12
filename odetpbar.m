function status = odetpbar(t, y, flag)
%ODETPBAR  Professional ODE integration progress bar (OutputFcn for ode15s)
%
%  Description:
%    Displays a real-time text progress bar during ode15s integration.
%    Shows a visual bar, percentage, elapsed time, and estimated time
%    remaining (ETA).  Updates are throttled to avoid excessive console
%    output.  Self-contained — no external dependencies.
%
%  Usage:
%    opts = odeset('OutputFcn', @odetpbar);
%    sol  = ode15s(odefun, tspan, y0, opts);
%
%  Display format:
%    [ODE]   [==============................]  46.7%  Elapsed: 12.3 s  ETA: 14.0 s
%
%  Note:
%    drawnow is required so that fprintf output is flushed to the Command
%    Window during ode15s execution.  Without it MATLAB buffers the output
%    and nothing is visible until the solver returns.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-12 — Professional formatting with ETA; self-contained

persistent tf tstart prev_len last_upd;

BAR_W    = 30;       % Progress bar width (characters)
THROTTLE = 0.25;     % Min seconds between display updates

if isempty(flag)
    % ===== Normal integration step =====
    now_s = toc(tstart);
    if ~isempty(last_upd) && (now_s - last_upd) < THROTTLE
        status = 0;
        return;
    end

    pct = min(100, 100 * mean(t) / tf);
    msg = build_msg(pct, now_s, BAR_W);

    % Overwrite previous line
    if ~isempty(prev_len) && prev_len > 0
        fprintf(repmat('\b', 1, prev_len));
    end
    fprintf('%s', msg);
    prev_len = numel(msg);
    last_upd = now_s;
    drawnow;

    status = 0;

else
    switch flag
        case 'init'
            tstart   = tic;
            tf       = max(t);
            last_upd = [];

            msg = build_msg(0, 0, BAR_W);
            fprintf('%s', msg);
            prev_len = numel(msg);
            drawnow;

        case 'done'
            elapsed = toc(tstart);

            % Overwrite with final 100% bar
            if ~isempty(prev_len) && prev_len > 0
                fprintf(repmat('\b', 1, prev_len));
            end
            bar_done = repmat('=', 1, BAR_W);
            msg = sprintf('[ODE]   [%s] 100.0%%  Elapsed: %s', ...
                bar_done, fmt_dur(elapsed));
            fprintf('%s\n', msg);
            drawnow;

            tf = []; tstart = []; prev_len = []; last_upd = [];

        otherwise
            error('odetpbar:UnknownFlag', 'Unknown flag: %s', flag);
    end
end

end


% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function msg = build_msg(pct, elapsed_s, bar_w)
%BUILD_MSG  Construct formatted progress bar string
    filled  = round(pct / 100 * bar_w);
    bar_str = [repmat('=', 1, filled), repmat('.', 1, bar_w - filled)];

    if pct > 0.5 && elapsed_s > 0.5
        eta_s = elapsed_s * (100 - pct) / pct;
        msg = sprintf('[ODE]   [%s] %5.1f%%  Elapsed: %s  ETA: %s', ...
            bar_str, pct, fmt_dur(elapsed_s), fmt_dur(eta_s));
    else
        msg = sprintf('[ODE]   [%s] %5.1f%%  Elapsed: %s', ...
            bar_str, pct, fmt_dur(elapsed_s));
    end
end


function s = fmt_dur(t)
%FMT_DUR  Convert seconds to human-readable duration string
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
