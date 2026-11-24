classdef BatchEtaEstimator < handle
%BATCHETAESTIMATOR  Parallel-aware batch ETA via Empirical Bayes + list scheduling
%
%  Description:
%    Computes accurate remaining wall-clock time (makespan ETA) for both
%    serial and parallel (parfor / parfeval) simulation batches.
%
%    Resolves four core failure modes of naive elapsed*(N-n)/n:
%      1. Wave quantization distortion: tracks active worker start times
%         and simulates remaining queue assignment across W workers with
%         greedy min-slot list scheduling, preventing the 10-30x ETA spike
%         right after the first job in a parallel wave completes.
%      2. In-flight worker progress: incorporates live per-job heartbeats
%         (pct, elapsed, job_eta_s) so accurate batch ETA is available even
%         before the first job finishes and throughout the tail wave.
%      3. Heterogeneous queues & LPT scheduling (e.g. paper_run): calibrates
%         per-group prior runtime estimates via Empirical Bayes shrinkage:
%           rho_g = (m0 + sum_{OK in g} actual) / (m0 + sum_{OK in g} prior)
%      4. Failure isolation: FAILED/TIMEOUT/ABORTED jobs are removed from
%         the pending queue without biasing the calibrated runtime model of
%         healthy remaining jobs.
%
%  Usage:
%    % Homogeneous parallel batch (e.g. batch_repeat_brownian):
%    be = BatchEtaEstimator(n_repeats, n_workers);
%    be.updateHeartbeat(k, pct, job_elapsed_s, job_eta_s, batch_wall_s);
%    be.markComplete(k, wall_s, is_ok, batch_wall_s);
%    est = be.estimate(batch_wall_s);   % struct: eta_s, finish_str, n_done, n_ok, ...
%
%    % Heterogeneous prioritised batch (e.g. paper_run):
%    be = BatchEtaEstimator(numel(J), n_workers, ...
%             'PriorSec', [J.est]*60, 'Groups', group_cell, 'JobIds', 1:numel(J));
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

    properties (SetAccess = private)
        n_total          (1,1) double = 0
        n_workers        (1,1) double = 1
        tic_start                       % uint64 tic handle

        % Per-job metadata (indexed 1..n_total)
        job_ids          (:,1) double   % Caller-visible job id (default 1:n_total)
        prior_s          (:,1) double   % Prior wall time [s] or relative weight
        has_prior        (1,1) logical = false
        prior_is_weights (1,1) logical = false   % True if prior_s is dimensionless relative weight
        group_idx        (:,1) double   % Group index 1..n_groups
        group_names      (:,1) cell

        % Job states: 0 = queued, 1 = running (known/inferred), 2 = completed
        state            (:,1) double
        start_wall       (:,1) double   % Batch wall time when job started running [s]
        actual_wall      (:,1) double   % Completed job wall time [s]
        ok_flag          (:,1) logical

        % Latest live heartbeat per job (index 1..n_total)
        hb_pct           (:,1) double
        hb_elapsed       (:,1) double
        hb_eta           (:,1) double
        hb_batch_wall    (:,1) double

        % Virtual worker slots (tracks active jobs when explicit start events are absent)
        slot_job         (1,:) double   % Job index currently occupying slot (0 = idle)
        slot_start       (1,:) double   % Batch wall time when current slot task started

        % Summary counters
        n_done           (1,1) double = 0
        n_ok             (1,1) double = 0
        n_failed         (1,1) double = 0
        next_queue_ptr   (1,1) double = 1

        % Empirical Bayes shrinkage strength (in multiples of median prior duration)
        eb_prior_jobs    (1,1) double = 2.0
    end

    methods
        function obj = BatchEtaEstimator(n_total, n_workers, varargin)
            obj.n_total   = max(0, round(double(n_total)));
            obj.n_workers = max(1, round(double(n_workers)));
            obj.tic_start = tic;

            obj.job_ids     = (1:obj.n_total)';
            obj.prior_s     = nan(obj.n_total, 1);
            obj.has_prior   = false;
            obj.group_idx   = ones(obj.n_total, 1);
            obj.group_names = {'default'};

            % Parse optional name-value pairs
            if ~isempty(varargin)
                if mod(numel(varargin), 2) ~= 0
                    error('BatchEtaEstimator:InvalidArgs', ...
                        'Optional arguments must be name-value pairs.');
                end
                for k = 1:2:numel(varargin)
                    key = lower(char(varargin{k}));
                    val = varargin{k + 1};
                    switch key
                        case 'priorsec'
                            p = double(val(:));
                            if numel(p) == obj.n_total && all(isfinite(p) & p > 0)
                                obj.prior_s          = p;
                                obj.has_prior        = true;
                                obj.prior_is_weights = false;
                            end
                        case 'priorweights'
                            w_p = double(val(:));
                            if numel(w_p) == obj.n_total && all(isfinite(w_p) & w_p > 0)
                                obj.prior_s          = w_p;
                                obj.has_prior        = true;
                                obj.prior_is_weights = true;
                            end
                        case 'groups'
                            g_cell = cellstr(val(:));
                            if numel(g_cell) == obj.n_total
                                [u_grp, ~, g_idx] = unique(g_cell, 'stable');
                                obj.group_names = u_grp(:);
                                obj.group_idx   = g_idx(:);
                            end
                        case 'jobids'
                            ids = double(val(:));
                            if numel(ids) == obj.n_total
                                obj.job_ids = ids;
                            end
                        case 'ticstart'
                            if ~isempty(val), obj.tic_start = val; end
                        case 'ebpriorjobs'
                            obj.eb_prior_jobs = max(double(val), 0.1);
                        otherwise
                            error('BatchEtaEstimator:UnknownOption', ...
                                'Unknown option ''%s''.', char(varargin{k}));
                    end
                end
            end

            obj.state         = zeros(obj.n_total, 1);
            obj.start_wall    = nan(obj.n_total, 1);
            obj.actual_wall   = nan(obj.n_total, 1);
            obj.ok_flag       = false(obj.n_total, 1);
            obj.hb_pct        = nan(obj.n_total, 1);
            obj.hb_elapsed    = nan(obj.n_total, 1);
            obj.hb_eta        = nan(obj.n_total, 1);
            obj.hb_batch_wall = nan(obj.n_total, 1);

            % Assume the first min(W, N) jobs enter execution at t = 0
            n_init = min(obj.n_workers, obj.n_total);
            obj.slot_job   = zeros(1, obj.n_workers);
            obj.slot_start = zeros(1, obj.n_workers);
            for w = 1:n_init
                obj.slot_job(w)   = w;
                obj.slot_start(w) = 0.0;
                obj.state(w)      = 1;
                obj.start_wall(w) = 0.0;
            end
            obj.next_queue_ptr = n_init + 1;
        end

        function updateHeartbeat(obj, job_id, pct, job_elapsed_s, job_eta_s, batch_wall_s)
            %UPDATEHEARTBEAT  Record live intra-job progress from a running worker.
            if nargin < 6 || isempty(batch_wall_s)
                batch_wall_s = toc(obj.tic_start);
            end
            idx = obj.findJobIndex(job_id);
            if isempty(idx) || obj.state(idx) == 2
                return;
            end

            obj.hb_pct(idx)        = double(pct);
            obj.hb_elapsed(idx)    = max(0.0, double(job_elapsed_s));
            obj.hb_batch_wall(idx) = max(0.0, double(batch_wall_s));
            if nargin >= 5 && ~isempty(job_eta_s) && isfinite(job_eta_s) && job_eta_s >= 0
                obj.hb_eta(idx) = double(job_eta_s);
            elseif pct > 0.2 && pct < 100 && job_elapsed_s > 0
                obj.hb_eta(idx) = job_elapsed_s * (100.0 - pct) / pct;
            end

            % Infer exact start time of this job in batch wall-clock coordinates
            inferred_start = max(0.0, batch_wall_s - obj.hb_elapsed(idx));
            if obj.state(idx) == 0
                % Job was dispatched out of nominal order (common in parfor):
                % place it in an active slot while preserving slot/state invariant.
                obj.attachRunningJobToSlot(idx, inferred_start);
            else
                obj.start_wall(idx) = inferred_start;
                w = find(obj.slot_job == idx, 1);
                if ~isempty(w)
                    obj.slot_start(w) = inferred_start;
                end
            end
        end

        function markComplete(obj, job_id, job_wall_s, is_ok, batch_wall_s)
            %MARKCOMPLETE  Record completion of one job and advance virtual worker slots.
            if nargin < 5 || isempty(batch_wall_s)
                batch_wall_s = toc(obj.tic_start);
            end
            idx = obj.findJobIndex(job_id);
            if isempty(idx)
                % Defensive: still advance counts if caller passed an out-of-range tag
                return;
            end
            if obj.state(idx) == 2
                return;   % already recorded
            end

            obj.state(idx)       = 2;
            obj.actual_wall(idx) = max(0.0, double(job_wall_s));
            obj.ok_flag(idx)     = logical(is_ok);
            obj.n_done           = obj.n_done + 1;
            if obj.ok_flag(idx)
                obj.n_ok = obj.n_ok + 1;
            else
                obj.n_failed = obj.n_failed + 1;
            end

            % Free the slot occupied by this job (or the slot whose start best matches)
            w = find(obj.slot_job == idx, 1);
            if isempty(w)
                active_slots = find(obj.slot_job > 0);
                if ~isempty(active_slots)
                    % Pick the active slot without a recent heartbeat that started earliest
                    no_hb = active_slots(isnan(obj.hb_batch_wall(obj.slot_job(active_slots))));
                    if ~isempty(no_hb)
                        [~, pos] = min(obj.slot_start(no_hb));
                        w = no_hb(pos);
                    else
                        [~, pos] = min(obj.slot_start(active_slots));
                        w = active_slots(pos);
                    end
                    % Move the displaced job back to queued state since it hasn't finished
                    displaced = obj.slot_job(w);
                    if displaced > 0 && obj.state(displaced) == 1
                        obj.state(displaced) = 0;
                    end
                end
            end

            if ~isempty(w)
                next_j = obj.popNextQueuedJob();
                if next_j > 0
                    obj.slot_job(w)         = next_j;
                    obj.slot_start(w)       = batch_wall_s;
                    obj.state(next_j)       = 1;
                    obj.start_wall(next_j)  = batch_wall_s;
                else
                    obj.slot_job(w)   = 0;
                    obj.slot_start(w) = batch_wall_s;
                end
            end
        end

        function est = estimate(obj, batch_wall_s)
            %ESTIMATE  Compute parallel makespan ETA via calibrated greedy list scheduling.
            if nargin < 2 || isempty(batch_wall_s)
                batch_wall_s = toc(obj.tic_start);
            end
            batch_wall_s = max(0.0, double(batch_wall_s));

            est = struct( ...
                'eta_s',        NaN, ...
                'elapsed_s',    batch_wall_s, ...
                'finish_str',   '', ...
                'n_done',       obj.n_done, ...
                'n_ok',         obj.n_ok, ...
                'n_failed',     obj.n_failed, ...
                'n_total',      obj.n_total, ...
                'eff_pct',      0.0, ...
                'mean_ok_s',    NaN, ...
                'scale_global', 1.0);

            if obj.n_total <= 0 || obj.n_done >= obj.n_total
                est.eta_s   = 0.0;
                est.eff_pct = 100.0;
                return;
            end

            % --- Step 1: Calibrate predicted total duration T_hat(1..N) ---
            [T_hat, mean_ok_s, global_scale] = obj.predictDurations();
            est.mean_ok_s    = mean_ok_s;
            est.scale_global = global_scale;

            if any(~isfinite(T_hat(obj.state ~= 2)))
                % No OK completions, no priors, and no usable heartbeats yet
                est.eff_pct = 100.0 * obj.n_done / obj.n_total;
                return;
            end

            % --- Step 2: Estimate remaining wall time on every active worker slot ---
            worker_free = zeros(1, obj.n_workers);
            inflight_frac_sum = 0.0;

            for w = 1:obj.n_workers
                j = obj.slot_job(w);
                if j > 0 && obj.state(j) == 1
                    [rem_j, frac_j] = obj.estimateActiveRemaining(j, T_hat(j), batch_wall_s);
                    worker_free(w) = rem_j;
                    inflight_frac_sum = inflight_frac_sum + frac_j;
                end
            end

            % Ensure every running job (in case more heartbeats arrived than slots) is accounted for
            queued_idx = find(obj.state == 0);

            % --- Step 3: Greedy min-slot list scheduling for all still-queued jobs ---
            for qi = 1:numel(queued_idx)
                j = queued_idx(qi);
                [~, w_min] = min(worker_free);
                worker_free(w_min) = worker_free(w_min) + T_hat(j);
            end

            eta_s = max(0.0, max(worker_free));
            est.eta_s      = eta_s;
            est.finish_str = ProgressEstimator.formatFinishTime(eta_s);

            % Effective completion percentage (work-weighted when priors exist)
            if obj.has_prior
                total_w = sum(T_hat);
                done_w  = sum(T_hat(obj.state == 2));
                for w = 1:obj.n_workers
                    j = obj.slot_job(w);
                    if j > 0 && obj.state(j) == 1
                        [~, frac_j] = obj.estimateActiveRemaining(j, T_hat(j), batch_wall_s);
                        done_w = done_w + frac_j * T_hat(j);
                    end
                end
                est.eff_pct = min(99.9, max(0.0, 100.0 * done_w / max(total_w, eps)));
            else
                est.eff_pct = min(99.9, max(0.0, 100.0 * (obj.n_done + inflight_frac_sum) / obj.n_total));
            end
        end
    end

    methods (Access = private)
        function idx = findJobIndex(obj, job_id)
            idx = find(obj.job_ids == double(job_id), 1);
        end

        function next_j = popNextQueuedJob(obj)
            next_j = 0;
            while obj.next_queue_ptr <= obj.n_total
                cand = obj.next_queue_ptr;
                obj.next_queue_ptr = obj.next_queue_ptr + 1;
                if obj.state(cand) == 0
                    next_j = cand;
                    return;
                end
            end
            % Fallback scan in case a job was returned to queued state
            rem_q = find(obj.state == 0, 1);
            if ~isempty(rem_q)
                next_j = rem_q;
            end
        end

        function attachRunningJobToSlot(obj, idx, inferred_start)
            % Prefer an idle slot first
            w = find(obj.slot_job == 0, 1);
            if ~isempty(w)
                obj.slot_job(w)     = idx;
                obj.slot_start(w)   = inferred_start;
                obj.state(idx)      = 1;
                obj.start_wall(idx) = inferred_start;
                return;
            end
            % Otherwise swap with an assumed-running slot that has no heartbeat yet
            active_slots = find(obj.slot_job > 0);
            no_hb = active_slots(isnan(obj.hb_batch_wall(obj.slot_job(active_slots))));
            if ~isempty(no_hb)
                w = no_hb(end);
                displaced = obj.slot_job(w);
                if displaced > 0 && obj.state(displaced) == 1
                    obj.state(displaced) = 0;
                end
                obj.slot_job(w)     = idx;
                obj.slot_start(w)   = inferred_start;
                obj.state(idx)      = 1;
                obj.start_wall(idx) = inferred_start;
            end
        end

        function [T_hat, mean_ok_s, global_scale] = predictDurations(obj)
            ok_mask = (obj.state == 2) & obj.ok_flag & isfinite(obj.actual_wall) & (obj.actual_wall > 0);
            ok_walls = obj.actual_wall(ok_mask);
            global_scale = 1.0;

            if ~isempty(ok_walls)
                mean_ok_s = BatchEtaEstimator.robustMean(ok_walls);
            else
                mean_ok_s = NaN;
            end

            if obj.has_prior
                if obj.prior_is_weights
                    % Dimensionless relative weights: calibrate absolute scale from OK runs
                    if any(ok_mask)
                        global_scale = sum(obj.actual_wall(ok_mask)) / max(sum(obj.prior_s(ok_mask)), eps);
                        T_hat = obj.prior_s * global_scale;
                    else
                        T_hat = nan(obj.n_total, 1);
                    end
                else
                    % Absolute prior in seconds: Empirical Bayes shrinkage per group
                    med_prior = median(obj.prior_s);
                    m0_global = obj.eb_prior_jobs * med_prior;
                    if any(ok_mask)
                        sum_act_all = sum(obj.actual_wall(ok_mask));
                        sum_pri_all = sum(obj.prior_s(ok_mask));
                        global_scale = (m0_global + sum_act_all) / (m0_global + sum_pri_all);
                        global_scale = min(max(global_scale, 0.05), 20.0);
                    end

                    T_hat = obj.prior_s * global_scale;
                    n_groups = numel(obj.group_names);
                    for g = 1:n_groups
                        g_ok = ok_mask & (obj.group_idx == g);
                        if any(g_ok)
                            g_priors = obj.prior_s(obj.group_idx == g);
                            m0_g = obj.eb_prior_jobs * median(g_priors);
                            sum_act_g = sum(obj.actual_wall(g_ok));
                            sum_pri_g = sum(obj.prior_s(g_ok));
                            % Shrink group factor toward global_scale
                            rho_g = (m0_g * global_scale + sum_act_g) / (m0_g + sum_pri_g);
                            rho_g = min(max(rho_g, 0.05), 20.0);
                            T_hat(obj.group_idx == g) = obj.prior_s(obj.group_idx == g) * rho_g;
                        end
                    end
                end
            else
                % Homogeneous batch (no prior): use OK completions or infer from active heartbeats
                if isfinite(mean_ok_s)
                    T_hat = repmat(mean_ok_s, obj.n_total, 1);
                else
                    hb_totals = [];
                    active_idx = find(obj.state == 1);
                    for k = 1:numel(active_idx)
                        j = active_idx(k);
                        if isfinite(obj.hb_elapsed(j)) && isfinite(obj.hb_eta(j)) && obj.hb_elapsed(j) > 0
                            hb_totals(end + 1) = obj.hb_elapsed(j) + obj.hb_eta(j); %#ok<AGROW>
                        elseif isfinite(obj.hb_pct(j)) && obj.hb_pct(j) >= 1.0 && isfinite(obj.hb_elapsed(j))
                            hb_totals(end + 1) = obj.hb_elapsed(j) * 100.0 / obj.hb_pct(j); %#ok<AGROW>
                        end
                    end
                    if ~isempty(hb_totals)
                        est_dur = median(hb_totals);
                        T_hat = repmat(est_dur, obj.n_total, 1);
                        mean_ok_s = est_dur;
                    else
                        T_hat = nan(obj.n_total, 1);
                    end
                end
            end

            % Keep completed jobs at their actual recorded durations for weighted progress
            if any(ok_mask)
                T_hat(ok_mask) = obj.actual_wall(ok_mask);
            end
        end

        function [rem_j, frac_j] = estimateActiveRemaining(obj, j, T_hat_j, batch_wall_s)
            % Prefer fresh live heartbeat when available
            if isfinite(obj.hb_batch_wall(j)) && isfinite(obj.hb_eta(j))
                age_s = max(0.0, batch_wall_s - obj.hb_batch_wall(j));
                % Decay reported job ETA by time elapsed since the heartbeat arrived,
                % with a smooth floor so a slow update never falsely claims 0 s remaining.
                raw_rem = obj.hb_eta(j) - age_s;
                if raw_rem > 0
                    rem_j = raw_rem;
                else
                    rem_j = max(5.0, 0.1 * obj.hb_eta(j));
                end
                tot_j = max(obj.hb_elapsed(j) + age_s + rem_j, eps);
                frac_j = min(0.99, max(0.0, (obj.hb_elapsed(j) + age_s) / tot_j));
                return;
            end

            % No heartbeat: use slot elapsed time and smooth heavy-tail overrun model
            e_j = max(0.0, batch_wall_s - obj.start_wall(j));
            if e_j < 0.85 * T_hat_j
                rem_j = T_hat_j - e_j;
            else
                % Overrunning job: smooth positive tail rather than clamping to 0
                rem_j = 0.15 * T_hat_j + 0.25 * (e_j - 0.85 * T_hat_j);
            end
            frac_j = min(0.95, max(0.0, e_j / max(e_j + rem_j, eps)));
        end
    end

    methods (Static, Access = private)
        function m = robustMean(vals)
            %ROBUSTMEAN  10% trimmed mean when n >= 5, plain mean otherwise.
            vals = sort(vals(:));
            n = numel(vals);
            if n >= 5
                k = floor(0.1 * n);
                m = mean(vals(k + 1 : n - k));
            else
                m = mean(vals);
            end
        end
    end
end
