function status = ode_progress_output(t, y, flag, tstart, tlim, silent, prog_ctx)
%ODE_PROGRESS_OUTPUT  Unified ODE OutputFcn: timeout guard, progress bar, file & queue heartbeats
%
%  Description:
%    Installed as ode15s' OutputFcn by run_simulation.m and
%    run_simulation_HPC_Linux.m.  Combines four responsibilities without any
%    persistent state (all telemetry state lives in prog_ctx.estimator):
%      1. Wall-clock timeout enforcement (returns status = 1 when elapsed
%         segment wall time reaches tlim, while still allowing flag='done'
%         to finalise the progress bar and telemetry cleanly).
%      2. Adaptive cost-rate EWMA updates on the shared ProgressEstimator.
%      3. Throttled atomic progress.txt writes (temp file + movefile) with
%         backward-compatible key=value fields plus [P10, P90], finish
%         clock, and ISO update timestamp.
%      4. Optional throttled parallel.pool.DataQueue heartbeat messages
%         carrying live pct, elapsed_s, smoothed eta_s, and finish_str
%         even during long single-segment non-Brownian integrations.
%
%  Syntax:
%    status = ode_progress_output(t, y, flag, tstart, tlim, silent, prog_ctx)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

status = 0;

has_ctx = (nargin >= 7) && isstruct(prog_ctx) && isfield(prog_ctx, 'estimator') && ~isempty(prog_ctx.estimator);

if has_ctx
    est   = prog_ctx.estimator;
    now_s = toc(est.tic_start);

    % Extract current dimensionless integration time
    t_cur = NaN;
    if strcmp(flag, 'init')
        if ~isempty(t), t_cur = t(1); else, t_cur = est.current_work; end
    elseif isempty(flag)
        if ~isempty(t), t_cur = t(end); end
    elseif strcmp(flag, 'done')
        t_cur = est.current_work;
    end

    if isfinite(t_cur)
        est.update(t_cur, now_s);
    end

    % --- Throttled atomic progress.txt heartbeat ---
    if isfield(prog_ctx, 'file') && ~isempty(prog_ctx.file)
        interval_s = 10.0;
        if isfield(prog_ctx, 'interval_s') && ~isempty(prog_ctx.interval_s)
            interval_s = prog_ctx.interval_s;
        end
        seg_num = 1;
        if isfield(prog_ctx, 'seg') && ~isempty(prog_ctx.seg)
            seg_num = prog_ctx.seg;
        end
        T_total = est.total_work;
        if isfield(prog_ctx, 'T_total') && ~isempty(prog_ctx.T_total)
            T_total = prog_ctx.T_total;
        end

        is_first_init   = strcmp(flag, 'init') && (seg_num == 1) && isinf(est.last_file_wall);
        is_final_done   = strcmp(flag, 'done') && ((est.current_work >= T_total - 1e-12) || (toc(tstart) >= tlim));
        if is_first_init || is_final_done
            do_write = true;
            est.markFileWritten(now_s);
        else
            do_write = est.shouldWriteFile(interval_s, now_s);
        end

        if do_write
            info = est.estimate(now_s);
            write_progress_file_atomic(prog_ctx.file, info, est.current_work, T_total, seg_num);
        end
    end

    % --- Optional throttled DataQueue heartbeat (HPC parallel batches) ---
    if isfield(prog_ctx, 'queue') && ~isempty(prog_ctx.queue)
        hb_interval_s = 30.0;
        if isfield(prog_ctx, 'hb_interval_s') && ~isempty(prog_ctx.hb_interval_s)
            hb_interval_s = prog_ctx.hb_interval_s;
        end
        if est.shouldSendHeartbeat(hb_interval_s, now_s)
            info = est.estimate(now_s);
            tag  = [];
            if isfield(prog_ctx, 'tag'), tag = prog_ctx.tag; end
            try
                send(prog_ctx.queue, struct( ...
                    'type',    'heartbeat', ...
                    'tag',     tag, ...
                    'pct',     info.pct, ...
                    'elapsed', info.elapsed_s, ...
                    'eta',     info.eta_s, ...
                    'finish',  info.finish_str));
            catch
                % Never abort an integration because of a telemetry send failure
            end
        end
    end
end

% Always let 'done' finalise the console progress bar even if timeout triggered
if strcmp(flag, 'done')
    if ~silent
        if has_ctx
            odetpbar(t, y, 'done', prog_ctx);
        else
            odetpbar(t, y, 'done');
        end
    end
    status = 0;
    return;
end

% Wall-clock segment timeout check
if toc(tstart) >= tlim
    status = 1;
    return;
end

if ~silent
    if has_ctx
        status = odetpbar(t, y, flag, prog_ctx);
    else
        status = odetpbar(t, y, flag);
    end
end

end


function write_progress_file_atomic(filepath, info, t_cur, T_total, seg_num)
%WRITE_PROGRESS_FILE_ATOMIC  Write progress.txt via temp file + atomic rename.
    try
        eta_s   = info.eta_s;
        p10_s   = info.eta_p10_s;
        p90_s   = info.eta_p90_s;
        fin_str = info.finish_str;
        if isempty(fin_str), fin_str = '--'; end
        % Replace space in multi-day finish string ('MM-dd HH:mm') with '_' so
        % whitespace-delimited key=value parsers keep every token intact.
        fin_tok = strrep(fin_str, ' ', '_');
        upd_tok = char(datetime('now', 'Format', 'uuuu-MM-dd''T''HH:mm:ss'));

        line = sprintf(['pct=%.2f t=%.6e T=%.6e elapsed=%.1f eta=%.1f ', ...
            'seg=%d eta_p10=%.1f eta_p90=%.1f finish=%s updated=%s\n'], ...
            info.pct, t_cur, T_total, info.elapsed_s, eta_s, ...
            seg_num, p10_s, p90_s, fin_tok, upd_tok);

        tmp_path = [filepath, '.tmp'];
        fid = fopen(tmp_path, 'w');
        if fid ~= -1
            fprintf(fid, '%s', line);
            fclose(fid);
            [ok, ~] = movefile(tmp_path, filepath, 'f');
            if ~ok
                % Fallback to direct overwrite if movefile fails (e.g. antivirus lock)
                fid2 = fopen(filepath, 'w');
                if fid2 ~= -1
                    fprintf(fid2, '%s', line);
                    fclose(fid2);
                end
                if isfile(tmp_path), delete(tmp_path); end
            end
        else
            fid2 = fopen(filepath, 'w');
            if fid2 ~= -1
                fprintf(fid2, '%s', line);
                fclose(fid2);
            end
        end
    catch
        % Never fail a simulation because of a progress file write error
    end
end
