classdef ProgressEstimator < handle
%PROGRESSESTIMATOR  Adaptive cost-rate EWMA progress and ETA estimator
%
%  Description:
%    Estimates remaining wall-clock time (ETA) and uncertainty intervals
%    for both continuous-time adaptive ODE integrations (ode15s) and
%    discrete loops (animation frames, flowfield chunks, post-processing).
%
%    Key design choices over naive cumulative linear extrapolation
%    (elapsed * (1 - p) / p):
%      1. Cost-rate space smoothing: smooths local slowness
%         c_k = d_wall / d_work [s / work_unit] rather than speed
%         v_k = d_work / d_wall.  Remaining time ETA = c_hat * rem_work is
%         linear in c, avoiding Jensen's inequality downward bias.
%      2. Continuous-time adaptive EWMA decay:
%           tau_eff = clamp(gamma * elapsed, tau_min, tau_max)
%           alpha_k = 1 - exp(-d_wall / tau_eff)
%         Fast startup forgetting discards initial stiffness transients and
%         JIT warmup; larger mid-run tau smooths periodic beating or
%         stochastic fluctuations while still adapting to regime shifts.
%      3. Huber / Winsorized outlier guard with regime-shift unlock:
%         Single spikes (basis rotations, OS stalls, isolated Jacobian
%         rebuilds) are clipped to [0.25*c_hat, 4*c_hat]; 3 consecutive
%         deviations unlock rapid adaptation to a true speed change.
%      4. Optional two-component Brownian/LaggedJacobian amortization:
%         Tracks fresh-Jacobian segments and normal lagged-Jacobian
%         segments separately so periodic Jacobian rebuilds do not cause
%         sawtooth ETA jitter.
%      5. Optional log-normal [P10, P90] uncertainty bounds and wall-clock
%         finish timestamp formatting.
%
%  Usage:
%    est = ProgressEstimator(T_total);
%    est.update(t_cur);              % reads internal tic
%    info = est.estimate();          % struct: pct, elapsed_s, eta_s, eta_p10_s, eta_p90_s, ...
%
%    % Deterministic / unit-test mode with explicit wall seconds:
%    est.update(work_done, wall_s);
%    info = est.estimate(wall_s);
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

    properties
        total_work       (1,1) double = 1.0
        start_work       (1,1) double = 0.0
        tau_min          (1,1) double = 15.0    % Minimum EWMA horizon [s]
        tau_max          (1,1) double = 600.0   % Maximum EWMA horizon [s]
        gamma            (1,1) double = 0.25    % Horizon fraction of elapsed time
        min_dt_sample    (1,1) double = 0.15    % Minimum wall window per EWMA update [s]
        warmup_s         (1,1) double = 0.5     % Minimum elapsed time before reporting ETA [s]
        warmup_pct       (1,1) double = 0.2     % Minimum percent before reporting ETA [%]
        outlier_lo       (1,1) double = 0.25    % Lower Winsorization factor relative to c_hat
        outlier_hi       (1,1) double = 4.0     % Upper Winsorization factor relative to c_hat
        regime_unlock_n  (1,1) double = 3       % Consecutive outliers to trigger regime shift
    end

    properties (SetAccess = private)
        tic_start                       % uint64 tic handle (used when wall_s not supplied)
        current_work     (1,1) double = 0.0
        last_sample_work (1,1) double = 0.0
        last_sample_wall (1,1) double = 0.0
        last_wall_s      (1,1) double = 0.0

        c_hat            (1,1) double = NaN     % Smoothed cost rate [s / work_unit]
        c_var            (1,1) double = 0.0     % EWMA variance of c_k
        n_updates        (1,1) double = 0
        n_eff            (1,1) double = 0.0     % Effective sample count for uncertainty
        n_consec_outlier (1,1) double = 0
        outlier_buf      (1,:) double = []      % Recent consecutive outlier c_obs for clean re-seed
        regime_start_wall (1,1) double = 0.0    % Wall time of last detected regime shift [s]

        % Optional two-component segment model (normal vs Jacobian-refresh segments)
        seg_w_norm       (1,1) double = NaN
        seg_var_norm     (1,1) double = 0.0
        seg_n_norm       (1,1) double = 0
        seg_w_jac        (1,1) double = NaN
        seg_var_jac      (1,1) double = 0.0
        seg_n_jac        (1,1) double = 0

        % Lightweight telemetry throttle state (replaces fragile persistent vars)
        last_file_wall   (1,1) double = -Inf
        last_hb_wall     (1,1) double = -Inf
    end

    methods
        function obj = ProgressEstimator(total_work, varargin)
            if nargin >= 1 && ~isempty(total_work)
                obj.total_work = max(double(total_work), eps);
            end
            obj.tic_start = tic;

            % Parse optional name-value pairs
            if ~isempty(varargin)
                if mod(numel(varargin), 2) ~= 0
                    error('ProgressEstimator:InvalidArgs', ...
                        'Optional arguments must be name-value pairs.');
                end
                for k = 1:2:numel(varargin)
                    key = lower(char(varargin{k}));
                    val = varargin{k + 1};
                    switch key
                        case 'startwork'
                            obj.start_work = double(val);
                            obj.current_work = obj.start_work;
                            obj.last_sample_work = obj.start_work;
                        case 'ticstart'
                            if ~isempty(val), obj.tic_start = val; end
                        case 'taumin'
                            obj.tau_min = max(double(val), 0.1);
                        case 'taumax'
                            obj.tau_max = max(double(val), obj.tau_min);
                        case 'gamma'
                            obj.gamma = max(double(val), 0.01);
                        case 'mindtsample'
                            obj.min_dt_sample = max(double(val), 0.0);
                        case 'warmups'
                            obj.warmup_s = max(double(val), 0.0);
                        case 'warmupsec'
                            obj.warmup_s = max(double(val), 0.0);
                        case 'warmuppct'
                            obj.warmup_pct = max(double(val), 0.0);
                        otherwise
                            error('ProgressEstimator:UnknownOption', ...
                                'Unknown option ''%s''.', char(varargin{k}));
                    end
                end
            end
        end

        function reset(obj, total_work, start_work, tic_handle)
            %RESET  Re-initialise state for a new run or loop.
            if nargin >= 2 && ~isempty(total_work)
                obj.total_work = max(double(total_work), eps);
            end
            if nargin >= 3 && ~isempty(start_work)
                obj.start_work = double(start_work);
            else
                obj.start_work = 0.0;
            end
            if nargin >= 4 && ~isempty(tic_handle)
                obj.tic_start = tic_handle;
            else
                obj.tic_start = tic;
            end
            obj.current_work     = obj.start_work;
            obj.last_sample_work = obj.start_work;
            obj.last_sample_wall = 0.0;
            obj.last_wall_s      = 0.0;
            obj.c_hat            = NaN;
            obj.c_var            = 0.0;
            obj.n_updates        = 0;
            obj.n_eff            = 0.0;
            obj.n_consec_outlier = 0;
            obj.outlier_buf      = [];
            obj.regime_start_wall = 0.0;
            obj.seg_w_norm       = NaN;
            obj.seg_var_norm     = 0.0;
            obj.seg_n_norm       = 0;
            obj.seg_w_jac        = NaN;
            obj.seg_var_jac      = 0.0;
            obj.seg_n_jac        = 0;
            obj.last_file_wall   = -Inf;
            obj.last_hb_wall     = -Inf;
        end

        function update(obj, work_val, wall_s)
            %UPDATE  Feed current cumulative work position and optional wall clock [s].
            if nargin < 3 || isempty(wall_s)
                wall_s = toc(obj.tic_start);
            end
            work_val = double(work_val);
            wall_s   = max(double(wall_s), obj.last_wall_s);

            if ~isfinite(work_val) || ~isfinite(wall_s)
                return;
            end

            % Keep monotonic frontier
            if work_val > obj.current_work
                obj.current_work = work_val;
            end
            obj.last_wall_s = wall_s;

            d_work = obj.current_work - obj.last_sample_work;
            d_wall = wall_s - obj.last_sample_wall;
            span   = max(obj.total_work - obj.start_work, eps);

            % Guard against stiff-ODE startup step discovery (e.g. h0 ~ 1e-17):
            % require meaningful work advancement before opening or updating the
            % cost-rate filter so d_wall / 1e-17 never poisons c_hat.
            min_dwork = 1e-4 * span;
            if obj.n_updates == 0
                min_dwork = max(min_dwork, (obj.warmup_pct / 100.0) * span);
            end

            if d_work < min_dwork || d_wall < obj.min_dt_sample
                return;
            end

            c_obs = d_wall / d_work;
            if ~isfinite(c_obs) || c_obs <= 0
                return;
            end

            if obj.n_updates == 0 || ~isfinite(obj.c_hat)
                obj.c_hat = c_obs;
                obj.c_var = 0.0;
                obj.n_updates = 1;
                obj.n_eff = 1.0;
                obj.n_consec_outlier = 0;
                obj.outlier_buf = [];
            else
                age_in_regime = max(wall_s - obj.regime_start_wall, 0.0);
                tau_eff = min(obj.tau_max, max(obj.tau_min, obj.gamma * age_in_regime));
                alpha = 1.0 - exp(-d_wall / tau_eff);
                alpha = min(max(alpha, 1e-4), 0.95);

                % Huber / Winsorized outlier protection with clean regime-shift re-seed
                lo = obj.outlier_lo * obj.c_hat;
                hi = obj.outlier_hi * obj.c_hat;
                is_outlier = (obj.n_updates >= 3) && (c_obs < lo || c_obs > hi);

                if is_outlier
                    obj.n_consec_outlier = obj.n_consec_outlier + 1;
                    obj.outlier_buf(end + 1) = c_obs;
                    if obj.n_consec_outlier < obj.regime_unlock_n
                        c_eff = min(max(c_obs, lo), hi);
                        % Prevent an isolated long stall from wiping out history via huge d_wall
                        alpha = min(alpha, 0.20);
                        diff_c = c_eff - obj.c_hat;
                        obj.c_hat = obj.c_hat + alpha * diff_c;
                        obj.c_var = (1.0 - alpha) * (obj.c_var + alpha * diff_c * diff_c);
                        obj.n_updates = obj.n_updates + 1;
                        obj.n_eff = (1.0 - alpha) * obj.n_eff + 1.0;
                    else
                        % Sustained regime shift confirmed: re-seed cleanly from new regime
                        obj.c_hat             = median(obj.outlier_buf);
                        obj.c_var             = var(obj.outlier_buf, 1);
                        obj.n_updates         = obj.n_updates + 1;
                        obj.n_eff             = double(numel(obj.outlier_buf));
                        obj.n_consec_outlier  = 0;
                        obj.outlier_buf       = [];
                        obj.regime_start_wall = wall_s;
                    end
                else
                    obj.n_consec_outlier = 0;
                    obj.outlier_buf      = [];
                    diff_c = c_obs - obj.c_hat;
                    obj.c_hat = obj.c_hat + alpha * diff_c;
                    % EWMA variance update (West 1979 / Finch 2009 form)
                    obj.c_var = (1.0 - alpha) * (obj.c_var + alpha * diff_c * diff_c);
                    obj.n_updates = obj.n_updates + 1;
                    obj.n_eff = (1.0 - alpha) * obj.n_eff + 1.0;
                end
            end

            obj.last_sample_work = obj.current_work;
            obj.last_sample_wall = wall_s;
        end

        function recordSegment(obj, seg_wall_s, is_jac_refresh)
            %RECORDSEGMENT  Track separate EWMA wall times for normal vs Jacobian-refresh segments.
            if ~isfinite(seg_wall_s) || seg_wall_s <= 0
                return;
            end
            alpha_seg = 0.15;
            if nargin >= 3 && logical(is_jac_refresh)
                if obj.seg_n_jac == 0 || ~isfinite(obj.seg_w_jac)
                    obj.seg_w_jac   = seg_wall_s;
                    obj.seg_var_jac = 0.0;
                else
                    a = max(alpha_seg, 1.0 / (obj.seg_n_jac + 1));
                    d = seg_wall_s - obj.seg_w_jac;
                    obj.seg_w_jac   = obj.seg_w_jac + a * d;
                    obj.seg_var_jac = (1.0 - a) * (obj.seg_var_jac + a * d * d);
                end
                obj.seg_n_jac = obj.seg_n_jac + 1;
            else
                if obj.seg_n_norm == 0 || ~isfinite(obj.seg_w_norm)
                    obj.seg_w_norm   = seg_wall_s;
                    obj.seg_var_norm = 0.0;
                else
                    a = max(alpha_seg, 1.0 / (obj.seg_n_norm + 1));
                    d = seg_wall_s - obj.seg_w_norm;
                    obj.seg_w_norm   = obj.seg_w_norm + a * d;
                    obj.seg_var_norm = (1.0 - a) * (obj.seg_var_norm + a * d * d);
                end
                obj.seg_n_norm = obj.seg_n_norm + 1;
            end
        end

        function info = estimate(obj, wall_s, rem_norm_segs, rem_jac_segs)
            %ESTIMATE  Return progress percentage, smoothed ETA, [P10, P90] band, and finish time.
            if nargin < 2 || isempty(wall_s)
                wall_s = toc(obj.tic_start);
            end
            wall_s = max(double(wall_s), obj.last_wall_s);

            span = max(obj.total_work - obj.start_work, eps);
            done = max(0.0, obj.current_work - obj.start_work);
            rem  = max(0.0, obj.total_work - obj.current_work);
            pct  = min(100.0, max(0.0, 100.0 * done / span));

            info = struct( ...
                'pct',         pct, ...
                'elapsed_s',   wall_s, ...
                'eta_s',       NaN, ...
                'eta_p10_s',   NaN, ...
                'eta_p90_s',   NaN, ...
                'c_hat',       obj.c_hat, ...
                'finish_str',  '', ...
                'status',      'warming_up', ...
                'method',      'none');

            if pct >= 100.0 - 1e-12 || rem <= 0
                info.pct       = 100.0;
                info.eta_s     = 0.0;
                info.eta_p10_s = 0.0;
                info.eta_p90_s = 0.0;
                info.status    = 'done';
                info.method    = 'complete';
                return;
            end

            % Optional two-component Brownian+LaggedJacobian segment estimate
            use_two_comp = (nargin >= 4) && ~isempty(rem_norm_segs) && ~isempty(rem_jac_segs) ...
                && isfinite(rem_norm_segs) && isfinite(rem_jac_segs) ...
                && (rem_norm_segs + rem_jac_segs > 0) ...
                && obj.seg_n_norm >= 2 && isfinite(obj.seg_w_norm);

            if use_two_comp
                w_norm = obj.seg_w_norm;
                if obj.seg_n_jac >= 1 && isfinite(obj.seg_w_jac)
                    w_jac   = obj.seg_w_jac;
                    var_jac = obj.seg_var_jac;
                else
                    w_jac   = w_norm;
                    var_jac = obj.seg_var_norm;
                end
                eta_s = max(0.0, rem_norm_segs * w_norm + rem_jac_segs * w_jac);
                % Propagate segment-level uncertainty (scale by sqrt(rem) for IID + rate drift)
                sig_norm = sqrt(max(obj.seg_var_norm, 0.0));
                sig_jac  = sqrt(max(var_jac, 0.0));
                % Blend independent noise sqrt(K) and correlated rate uncertainty K/sqrt(n)
                sigma_eta = sqrt( ...
                    rem_norm_segs * (sig_norm^2) + ...
                    (rem_norm_segs^2) * (sig_norm^2) / max(obj.seg_n_norm, 1) + ...
                    rem_jac_segs  * (sig_jac^2)  + ...
                    (rem_jac_segs^2)  * (sig_jac^2)  / max(obj.seg_n_jac, 1));
                info.eta_s  = eta_s;
                info.method = 'two_component_seg';
                [info.eta_p10_s, info.eta_p90_s] = ProgressEstimator.logNormalBounds(eta_s, sigma_eta);
                info.status     = 'ok';
                info.finish_str = ProgressEstimator.formatFinishTime(eta_s);
                return;
            end

            % Standard continuous-time cost-rate EWMA estimate
            if wall_s < obj.warmup_s || pct < obj.warmup_pct
                return;
            end

            if obj.n_updates >= 1 && isfinite(obj.c_hat) && obj.c_hat > 0
                c_use = obj.c_hat;
                % Include any un-flushed partial progress since last_sample_wall
                d_wall_tail = wall_s - obj.last_sample_wall;
                d_work_tail = obj.current_work - obj.last_sample_work;
                if d_wall_tail > 2.0 * max(obj.min_dt_sample, 0.5) && d_work_tail > 0
                    % Smooth preview including current partial window
                    tau_eff = min(obj.tau_max, max(obj.tau_min, obj.gamma * max(wall_s, obj.tau_min)));
                    a_tail = 1.0 - exp(-d_wall_tail / tau_eff);
                    c_tail = d_wall_tail / d_work_tail;
                    c_tail = min(max(c_tail, obj.outlier_lo * c_use), obj.outlier_hi * c_use);
                    c_use  = c_use + a_tail * (c_tail - c_use);
                elseif d_wall_tail > 30.0 && d_work_tail <= 0
                    % Stalled without advancing work: do not under-report remaining time
                    c_use = max(c_use, wall_s / max(done, eps));
                end
                eta_s = max(0.0, c_use * rem);
                info.method = 'cost_rate_ewma';

                % Uncertainty from EWMA variance of c_k
                if obj.n_updates >= 3 && obj.c_var > 0
                    sigma_c   = sqrt(obj.c_var);
                    sigma_eta = sigma_c * rem / sqrt(max(obj.n_eff / 2.0, 1.0));
                else
                    sigma_eta = 0.25 * eta_s;
                end
            elseif done > 0
                % Fallback before first EWMA window closes
                eta_s     = max(0.0, wall_s * rem / done);
                sigma_eta = 0.35 * eta_s;
                info.method = 'cumulative_fallback';
            else
                return;
            end

            info.eta_s = eta_s;
            [info.eta_p10_s, info.eta_p90_s] = ProgressEstimator.logNormalBounds(eta_s, sigma_eta);
            info.status     = 'ok';
            info.finish_str = ProgressEstimator.formatFinishTime(eta_s);
        end

        function tf = shouldWriteFile(obj, interval_s, wall_s)
            if nargin < 3 || isempty(wall_s)
                wall_s = toc(obj.tic_start);
            end
            tf = (wall_s - obj.last_file_wall) >= interval_s;
            if tf
                obj.last_file_wall = wall_s;
            end
        end

        function markFileWritten(obj, wall_s)
            if nargin < 2 || isempty(wall_s)
                wall_s = toc(obj.tic_start);
            end
            obj.last_file_wall = wall_s;
        end

        function tf = shouldSendHeartbeat(obj, interval_s, wall_s)
            if nargin < 3 || isempty(wall_s)
                wall_s = toc(obj.tic_start);
            end
            tf = (wall_s - obj.last_hb_wall) >= interval_s;
            if tf
                obj.last_hb_wall = wall_s;
            end
        end
    end

    methods (Static)
        function [p10, p90] = logNormalBounds(eta_s, sigma_eta)
            %LOGNORMALBOUNDS  Strictly positive [P10, P90] interval for ETA.
            if ~isfinite(eta_s) || eta_s <= 0
                p10 = eta_s;
                p90 = eta_s;
                return;
            end
            if ~isfinite(sigma_eta) || sigma_eta <= 0
                p10 = eta_s;
                p90 = eta_s;
                return;
            end
            cv = min(max(sigma_eta / max(eta_s, eps), 0.0), 1.5);
            sigma_ln = sqrt(log(1.0 + cv * cv));
            z80 = 1.28155;   % 10th / 90th percentile of N(0, 1)
            factor = exp(z80 * sigma_ln);
            p10 = eta_s / factor;
            p90 = eta_s * factor;
        end

        function s = formatFinishTime(eta_s)
            %FORMATFINISHTIME  Wall-clock completion timestamp ('HH:MM' or 'MM-DD HH:MM').
            if ~isfinite(eta_s) || eta_s < 0 || eta_s > 30 * 86400
                s = '';
                return;
            end
            now_dt = datetime('now');
            fin_dt = now_dt + seconds(eta_s);
            if year(fin_dt) == year(now_dt) && month(fin_dt) == month(now_dt) && day(fin_dt) == day(now_dt)
                s = char(datetime(fin_dt, 'Format', 'HH:mm'));
            else
                s = char(datetime(fin_dt, 'Format', 'MM-dd HH:mm'));
            end
        end

        function s = formatEtaWithRange(info)
            %FORMATETAWITHRANGE  Human-readable ETA string with optional finish clock & range.
            if ~isstruct(info) || ~isfield(info, 'eta_s') || ~isfinite(info.eta_s)
                s = 'warming up';
                return;
            end
            base = format_duration(info.eta_s, 'compact');
            % Show [P10-P90] range only when ETA >= 120 s and spread is meaningful (> 15%)
            show_range = info.eta_s >= 120 && isfield(info, 'eta_p10_s') && isfield(info, 'eta_p90_s') ...
                && isfinite(info.eta_p10_s) && isfinite(info.eta_p90_s) ...
                && ((info.eta_p90_s - info.eta_p10_s) / max(info.eta_s, eps) > 0.15);
            if show_range
                base = sprintf('%s [%s..%s]', base, ...
                    format_duration(info.eta_p10_s, 'compact'), ...
                    format_duration(info.eta_p90_s, 'compact'));
            end
            if info.eta_s >= 60 && isfield(info, 'finish_str') && ~isempty(info.finish_str)
                s = sprintf('%s (~%s)', base, info.finish_str);
            else
                s = base;
            end
        end
    end
end
