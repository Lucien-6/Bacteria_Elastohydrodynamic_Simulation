function summary = cpc_parallel_and_flowfield_benchmark(varargin)
%CPC_PARALLEL_AND_FLOWFIELD_BENCHMARK  Multi-worker parfor & GPU/ThreadPool flow-field benchmark
%
%  Description:
%    Complements validation/cpc_performance_benchmark.m for Section 5 of the
%    BacEHD CPC paper by measuring:
%      Part A — Multi-worker parfor strong scaling for independent Brownian
%               ensembles (N_runs = 16, N_body = 306, 1 pilus / 30 segs,
%               T_obs = 0.4 s, dt_B = 0.05 s, 1 BLAS thread per worker),
%               reporting compute wall time, speedup, parallel efficiency,
%               pool startup overhead, and cross-worker trajectory parity.
%      Part B — 6-panel 3-D/2-D flow-field reconstruction benchmark
%               (force extraction + 3-plane velocity kernels + 2-D planar
%               streamlines) comparing serial CPU, MATLAB ThreadPool CPU,
%               and hybrid CUDA FP64 GPU + ThreadPool on a real saved
%               3-D microswimmer state (E. coli single polar flagellum),
%               with exact FP64 GPU-vs-CPU velocity accuracy verification.
%
%  Syntax:
%    summary = cpc_parallel_and_flowfield_benchmark()
%    summary = cpc_parallel_and_flowfield_benchmark('OutDir', outdir)
%
%  Outputs:
%    summary - struct with .parallel (table), .flowfield (table),
%              .kernel_fp64_rel_err, .hardware, .passed, .elapsed_s.
%              Tables written to:
%                <OutDir>/parallel_scaling_benchmark.csv
%                <OutDir>/flowfield_gpu_benchmark.csv
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-25

ip = inputParser;
addParameter(ip, 'OutDir', '');
addParameter(ip, 'NRuns', 32, @(x) isnumeric(x) && isscalar(x) && x >= 4);
addParameter(ip, 'NFrames', 12, @(x) isnumeric(x) && isscalar(x) && x >= 3);
addParameter(ip, 'SkipParallel', false, @(x) islogical(x) || isnumeric(x));
parse(ip, varargin{:});

t_all = tic;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(genpath(fullfile(root, 'src')));
addpath(fullfile(root, 'chlamydomonas'));
addpath(fullfile(root, 'chlamydomonas', 'core'));
addpath(fullfile(root, 'MATLAB Add-On', ...
    'Suite of functions to perform uniform sampling of a sphere'));

if isempty(ip.Results.OutDir)
    outdir = fullfile(root, 'output', 'cpc_checks');
else
    outdir = ip.Results.OutDir;
end
if ~isfolder(outdir)
    mkdir(outdir);
end

n_runs   = ip.Results.NRuns;
n_frames = ip.Results.NFrames;
n_cores  = feature('numcores');
gpu_name = 'none';
has_gpu  = false;
try
    gdev = gpuDevice();
    has_gpu  = true;
    gpu_name = char(gdev.Name);
catch
    has_gpu = false;
end
has_pct = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));

fprintf('================================================================\n');
fprintf('  BacEHD v%s — Parallel & Flow-Field Scaling Benchmark\n', project_version());
fprintf('  Platform: MATLAB %s on %s | Physical cores: %d | GPU: %s\n', ...
    version, computer, n_cores, gpu_name);
fprintf('================================================================\n');

%% ---- Part A: Multi-worker parfor strong scaling (Brownian ensemble) ----
nolog = @(varargin) [];
par_csv = fullfile(outdir, 'parallel_scaling_benchmark.csv');
if ip.Results.SkipParallel && isfile(par_csv)
    fprintf('\n--- Part A: Loading existing %s ---\n', par_csv);
    T_par = readtable(par_csv);
else
fprintf('\n--- Part A: Brownian ensemble strong scaling (N_runs = %d) ---\n', n_runs);
[params_base, Z0_base] = build_brownian_case(nolog);
T_obs_s = 0.25;
dt_B_s  = 0.05;
seeds   = 2026112500 + (1:n_runs);

% Baseline 1: serial with 1 BLAS thread (apples-to-apples reference for parfor)
prev_blas = maxNumCompThreads(1);
cleanup_blas = onCleanup(@() maxNumCompThreads(prev_blas));

% Warm up JIT on 1 short step
[~] = run_one_brownian_repeat(params_base, Z0_base, dt_B_s, dt_B_s, 99999);

t0 = tic;
end_pos_ref = zeros(n_runs, 3);
for k = 1:n_runs
    end_pos_ref(k, :) = run_one_brownian_repeat(params_base, Z0_base, T_obs_s, dt_B_s, seeds(k));
end
wall_1_1t = toc(t0);

% Also record serial with default multi-threaded BLAS (for full transparency)
maxNumCompThreads(prev_blas);
t0 = tic;
end_pos_def = zeros(n_runs, 3);
for k = 1:n_runs
    end_pos_def(k, :) = run_one_brownian_repeat(params_base, Z0_base, T_obs_s, dt_B_s, seeds(k));
end
wall_1_def = toc(t0);
maxNumCompThreads(1);

p_rows = {};
p_rows{end+1} = table("serial_default_blas", 1, prev_blas, n_runs, 0, ...
    wall_1_def, wall_1_1t / wall_1_def, 100 * (wall_1_1t / wall_1_def), ...
    max(vecnorm(end_pos_def - end_pos_ref, 2, 2)), ...
    'VariableNames', {'mode', 'n_workers', 'blas_threads_per_worker', 'n_runs', ...
    'pool_start_s', 'compute_wall_s', 'speedup_vs_1w1t', 'efficiency_pct', 'max_pos_diff_L0'});
p_rows{end+1} = table("serial_1thread", 1, 1, n_runs, 0, ...
    wall_1_1t, 1.0, 100.0, 0.0, ...
    'VariableNames', p_rows{1}.Properties.VariableNames);

fprintf('  %-20s | w=%2d (%2d BLAS) | compute = %6.2f s | speedup = %5.2fx | eff = %5.1f%%\n', ...
    'serial_default_blas', 1, prev_blas, wall_1_def, wall_1_1t / wall_1_def, 100 * (wall_1_1t / wall_1_def));
fprintf('  %-20s | w=%2d (%2d BLAS) | compute = %6.2f s | speedup = %5.2fx | eff = %5.1f%%\n', ...
    'serial_1thread', 1, 1, wall_1_1t, 1.0, 100.0);

if has_pct
    cand_w = [2, 4, 8, 16];
    cand_w = cand_w(cand_w <= n_cores);
    delete(gcp('nocreate'));
    cl = parcluster('Processes');
    if cl.NumWorkers < max(cand_w)
        cl.NumWorkers = max(cand_w);
    end
    for iw = 1:numel(cand_w)
        w = cand_w(iw);
        t_pool = tic;
        pool = parpool(cl, w);
        pool_start_s = toc(t_pool);

        % Match batch_run_from_json.m: 1 BLAS thread per worker + local body_precomp
        wait(parfevalOnAll(pool, @maxNumCompThreads, 0, 1));
        p_nopre = params_base;
        p_nopre.body_precomp = [];
        p_const = parallel.pool.Constant(@() rebuild_precomp_on_worker(p_nopre));
        % Apples-to-apples JIT warm-up on every worker (same 1-step warmup as serial)
        parfor k = 1:w
            pw = p_const.Value;
            [~] = run_one_brownian_repeat(pw, Z0_base, dt_B_s, dt_B_s, 99999);
        end

        t0 = tic;
        end_pos_w = zeros(n_runs, 3);
        parfor k = 1:n_runs
            pw = p_const.Value;
            end_pos_w(k, :) = run_one_brownian_repeat(pw, Z0_base, T_obs_s, dt_B_s, seeds(k));
        end
        wall_w = toc(t0);
        max_diff = max(vecnorm(end_pos_w - end_pos_ref, 2, 2));
        sp = wall_1_1t / wall_w;
        eff = 100 * sp / w;

        p_rows{end+1} = table(sprintf("parfor_%dw", w), w, 1, n_runs, pool_start_s, ...
            wall_w, sp, eff, max_diff, ...
            'VariableNames', p_rows{1}.Properties.VariableNames);
        fprintf('  %-20s | w=%2d ( 1 BLAS) | compute = %6.2f s | speedup = %5.2fx | eff = %5.1f%% | pool = %4.1f s | diff = %.1e\n', ...
            sprintf('parfor_%dw', w), w, wall_w, sp, eff, pool_start_s, max_diff);

        delete(pool);
    end
end
clear cleanup_blas;
maxNumCompThreads(prev_blas);

T_par = vertcat(p_rows{:});
writetable(T_par, par_csv);
end

%% ---- Part B: 6-panel 3-D/2-D Flow-Field Reconstruction Benchmark ----
fprintf('\n--- Part B: Flow-field reconstruction benchmark (N_frames = %d) ---\n', n_frames);
ecoli_mat = fullfile(root, 'output', 'Ecoli-Single-Polar-Flagellum', 'simulation_results.mat');
if ~isfile(ecoli_mat)
    error('cpc_parallel_and_flowfield_benchmark:MissingCase', ...
        'Expected saved E. coli simulation at %s', ecoli_mat);
end
ws = warning('off', 'MATLAB:decomposition:LoadNotSupported');
S_ecoli = load(ecoli_mat, 'results', 'params');
warning(ws);
res_e = S_ecoli.results;
par_e = S_ecoli.params;
par_e.body_precomp = build_body_precomp(par_e.body_local, par_e.eps_body);

n_avail = res_e.flow_state.n_times;
frame_idx = unique(round(linspace(max(2, round(0.1 * n_avail)), n_avail, n_frames)));
n_frames_actual = numel(frame_idx);

flow_opts_cpu = struct('grid_spacing_um', 0.05, 'grid_range_um', 10.0, ...
    'quiver_skip', 5, 'chunk_size', 2000, 'use_gpu', false);
flow_opts_gpu = flow_opts_cpu;
flow_opts_gpu.use_gpu = true;

% B1: Single-frame velocity kernel CPU vs GPU + exact FP64 accuracy check
k0 = frame_idx(end);
[Z_k0, par_k0, t_k0] = unpack_flow_frame(res_e, par_e, k0);
hydro_k0 = flowfield_from_saved_state(Z_k0, par_k0, t_k0, nolog);

% Warm-up CPU & GPU
pl_cpu = flowfield_grid_velocity(hydro_k0, par_k0, flow_opts_cpu, nolog, true);
n_ext_pts = sum(~pl_cpu.xz.mask(:)) + sum(~pl_cpu.yz.mask(:)) + sum(~pl_cpu.xy.mask(:));

n_rep_k = 3;
t0 = tic;
for r = 1:n_rep_k
    pl_cpu = flowfield_grid_velocity(hydro_k0, par_k0, flow_opts_cpu, nolog, true);
end
kern_cpu_ms = 1e3 * toc(t0) / n_rep_k;

kern_gpu_ms = NaN;
fp64_rel_err = NaN;
if has_gpu
    pl_gpu = flowfield_grid_velocity(hydro_k0, par_k0, flow_opts_gpu, nolog, true); % warm-up
    wait(gpuDevice());
    t0 = tic;
    for r = 1:n_rep_k
        pl_gpu = flowfield_grid_velocity(hydro_k0, par_k0, flow_opts_gpu, nolog, true);
        wait(gpuDevice());
    end
    kern_gpu_ms = 1e3 * toc(t0) / n_rep_k;
    fp64_rel_err = max_plane_rel_diff(pl_cpu, pl_gpu);
    fprintf('  Kernel 3-plane (%d exterior pts, Nb=%d, Nseg=%d): CPU = %.1f ms -> GPU = %.1f ms (%.2fx) | FP64 rel err = %.2e\n', ...
        n_ext_pts, par_e.N_body, sum(par_e.N_segments), kern_cpu_ms, kern_gpu_ms, ...
        kern_cpu_ms / kern_gpu_ms, fp64_rel_err);
end

% B2: Multi-frame batch reconstruction (force + 3-plane velocity + 2-D streamlines)
f_rows = {};

% Mode 1: Serial CPU
[t_force_s, t_vel_s, t_str_s, slim_mib] = bench_flowfield_serial( ...
    res_e, par_e, frame_idx, flow_opts_cpu, nolog);
t_tot_serial = t_force_s + t_vel_s + t_str_s;
f_rows{end+1} = table("cpu_serial", n_frames_actual, n_ext_pts, 1, false, ...
    kern_cpu_ms, t_force_s, t_vel_s, t_str_s, t_tot_serial, 1.0, slim_mib, 0.0, ...
    'VariableNames', {'mode', 'n_frames', 'n_exterior_pts', 'n_threads', 'use_gpu', ...
    'kernel_ms', 'force_s', 'velocity_s', 'streamlines_s', 'total_s', ...
    'speedup_vs_serial', 'slim_cache_MiB_per_frame', 'fp64_rel_err'});
fprintf('  %-20s | force=%5.2fs | vel=%5.2fs | stream=%5.2fs | total=%5.2fs (1.00x) | cache=%.2f MiB/frame\n', ...
    'cpu_serial', t_force_s, t_vel_s, t_str_s, t_tot_serial, slim_mib);

% Mode 2 & 3: ThreadPool CPU and Hybrid GPU + ThreadPool
if has_pct
    delete(gcp('nocreate'));
    tpool = parpool('Threads');
    n_th = tpool.NumWorkers;

    % Warm up thread pool
    parfor i = 1:n_th
        tmp = sin(i); %#ok<NASGU>
    end

    [t_force_tp, t_vel_tp, t_str_tp, ~] = bench_flowfield_threadpool_cpu( ...
        res_e, par_e, frame_idx, flow_opts_cpu, nolog);
    t_tot_tp = t_force_tp + t_vel_tp + t_str_tp;
    f_rows{end+1} = table("cpu_threadpool", n_frames_actual, n_ext_pts, n_th, false, ...
        kern_cpu_ms, t_force_tp, t_vel_tp, t_str_tp, t_tot_tp, t_tot_serial / t_tot_tp, slim_mib, 0.0, ...
        'VariableNames', f_rows{1}.Properties.VariableNames);
    fprintf('  %-20s | force=%5.2fs | vel=%5.2fs | stream=%5.2fs | total=%5.2fs (%.2fx) | threads=%d\n', ...
        'cpu_threadpool', t_force_tp, t_vel_tp, t_str_tp, t_tot_tp, t_tot_serial / t_tot_tp, n_th);

    if has_gpu
        [t_force_g, t_vel_g, t_str_g, ~] = bench_flowfield_gpu_threadpool( ...
            res_e, par_e, frame_idx, flow_opts_gpu, nolog);
        t_tot_g = t_force_g + t_vel_g + t_str_g;
        f_rows{end+1} = table("gpu_plus_threadpool", n_frames_actual, n_ext_pts, n_th, true, ...
            kern_gpu_ms, t_force_g, t_vel_g, t_str_g, t_tot_g, t_tot_serial / t_tot_g, slim_mib, fp64_rel_err, ...
            'VariableNames', f_rows{1}.Properties.VariableNames);
        fprintf('  %-20s | force=%5.2fs | vel=%5.2fs | stream=%5.2fs | total=%5.2fs (%.2fx) | FP64 err=%.1e\n', ...
            'gpu_plus_threadpool', t_force_g, t_vel_g, t_str_g, t_tot_g, t_tot_serial / t_tot_g, fp64_rel_err);
    end
    delete(tpool);
end

T_flow = vertcat(f_rows{:});
flow_csv = fullfile(outdir, 'flowfield_gpu_benchmark.csv');
writetable(T_flow, flow_csv);

passed = all(T_par.max_pos_diff_L0 < 1e-10) && (isnan(fp64_rel_err) || fp64_rel_err < 1e-13);
elapsed_s = toc(t_all);

summary = struct( ...
    'parallel',           T_par, ...
    'flowfield',          T_flow, ...
    'kernel_fp64_rel_err', fp64_rel_err, ...
    'hardware',           struct('matlab', version, 'computer', computer, ...
                                 'physical_cores', n_cores, 'gpu', gpu_name), ...
    'passed',             passed, ...
    'elapsed_s',          elapsed_s);

fprintf('----------------------------------------------------------------\n');
fprintf('Parallel & Flow-Field Benchmark %s in %.1f s.\n  Tables: %s\n          %s\n', ...
    tern(passed, 'PASSED', 'FAILED'), elapsed_s, par_csv, flow_csv);
end


function [params, Z0] = build_brownian_case(nolog)
pili_cfg = struct('type', 'linear', 'length_um', 3, 'site', 2);
[Pili_Matrix, ext] = parse_pili_config(pili_cfg);
params = init_coupled_params( ...
    'geometry_mode',       'parametric', ...
    'minor_axis',          0.4e-6, ...
    'major_axis',          1.25e-6, ...
    'N_body',              300, ...
    'alpha_0',             pi / 3, ...
    'E',                   3.279e8, ...
    'r_pili',              3e-9, ...
    'sigma',               0.2, ...
    'mu',                  1e-3, ...
    'rho_bac',             1180, ...
    'rho_fluid',           998.232, ...
    'Temper',              293.15, ...
    'hydro_radius_factor', sqrt(2), ...
    'seg_per_um',          10, ...
    'N_segments',          30, ...
    'Pili_Matrix',         Pili_Matrix, ...
    'root_coupling',       {ext.root_coupling}, ...
    'motor_torque',        {ext.motor_torque}, ...
    'logf',                nolog, ...
    'logf_file',           nolog);
Z0 = init_coupled_state(params, nolog);
[Z0, params] = adaptive_basis_check(Z0, params);
end


function p = rebuild_precomp_on_worker(p)
p.body_precomp = build_body_precomp(p.body_local, p.eps_body);
end


function pos_end = run_one_brownian_repeat(params, Z0, T_obs_s, dt_B_s, seed)
n_steps = round(T_obs_s / dt_B_s);
dt_B_star = dt_B_s / params.scales.T_0;
opts_base = odeset('AbsTol', 1e-4, 'RelTol', 1e-4, 'Stats', 'off');
rng_s = RandStream('mt19937ar', 'Seed', seed);
jc = LaggedJacobian(50);
Z = Z0;
p = params;
for s = 1:n_steps
    p.brown_force_dimless = sample_brownian_body_force(Z, p, dt_B_s, rng_s);
    rhs_f = @(t, ZZ) ode_coupled_rhs(t, ZZ, p);
    jac_f = @(t, ZZ) ode_coupled_jacobian(t, ZZ, p);
    jc.tick();
    jc.beginSegment(rhs_f, jac_f);
    opts = odeset(opts_base, 'Jacobian', @(t, ZZ) jc.evaluate(t, ZZ));
    sol = ode15s(rhs_f, [0, dt_B_star], Z, opts);
    Z = sol.y(:, end);
    Z(4:7) = Z(4:7) / norm(Z(4:7));
    Z = enforce_clamped_geometry(Z, p);
end
pos_end = Z(1:3).';
end


function [Z, p, t_dl] = unpack_flow_frame(results, params, k)
spec = flowfield_state_data(results, params, k, 'history');
Z = spec.Z;
p = spec.params;
t_dl = spec.t_dimless;
end


function err = max_plane_rel_diff(a, b)
diffs = [ ...
    plane_diff(a.xz.Ux, a.xz.Uz, b.xz.Ux, b.xz.Uz, a.xz.mask), ...
    plane_diff(a.yz.Uy, a.yz.Uz, b.yz.Uy, b.yz.Uz, a.yz.mask), ...
    plane_diff(a.xy.Ux, a.xy.Uy, b.xy.Ux, b.xy.Uy, a.xy.mask)];
err = max(diffs);
end


function d = plane_diff(u1, v1, u2, v2, mask)
ok = ~mask & isfinite(u1) & isfinite(v1);
num = max(hypot(u1(ok) - u2(ok), v1(ok) - v2(ok)));
den = max(hypot(u1(ok), v1(ok)));
d = num / max(den, eps);
end


function [t_force, t_vel, t_str, slim_mib] = bench_flowfield_serial(res, par, frame_idx, fopts, nolog)
nf = numel(frame_idx);
hydros = cell(nf, 1);
t0 = tic;
for i = 1:nf
    [Z, p, t_dl] = unpack_flow_frame(res, par, frame_idx(i));
    hydros{i} = flowfield_from_saved_state(Z, p, t_dl, nolog);
end
t_force = toc(t0);

planes = cell(nf, 1);
t0 = tic;
for i = 1:nf
    planes{i} = flowfield_grid_velocity(hydros{i}, par, fopts, nolog, true);
end
t_vel = toc(t0);

streams = cell(nf, 1);
t0 = tic;
for i = 1:nf
    streams{i} = compute_three_plane_streams(planes{i}, hydros{i});
end
t_str = toc(t0);

% Measure slim per-frame cache footprint (speed grids + streamlines, matching v10.2.0)
p1 = planes{1};
slim = struct('Umag_xz', single(p1.xz.Umag), 'Umag_yz', single(p1.yz.Umag), ...
    'Umag_xy', single(p1.xy.Umag), 'streams', streams{1}); %#ok<NASGU>
info_w = whos('slim');
slim_mib = info_w.bytes / 1024^2;
end


function [t_force, t_vel, t_str, slim_mib] = bench_flowfield_threadpool_cpu(res, par, frame_idx, fopts, nolog)
nf = numel(frame_idx);
Z_all = cell(nf, 1);
t_all = zeros(nf, 1);
rot_all = cell(nf, 1);
for i = 1:nf
    [Z_all{i}, pi_tmp, t_all(i)] = unpack_flow_frame(res, par, frame_idx(i));
    rot_all{i} = pi_tmp.rot;
end
p_nopre = par;
p_nopre.body_precomp = [];
p_const = parallel.pool.Constant(@() rebuild_precomp_on_worker(p_nopre));
% Warm up Constant and JIT on workers before timing
parfor i = 1:min(nf, 8)
    pw = p_const.Value;
    pw.rot = rot_all{1};
    [~] = flowfield_from_saved_state(Z_all{1}, pw, t_all(1), nolog);
end

hydros = cell(nf, 1);
t0 = tic;
parfor i = 1:nf
    pw = p_const.Value;
    pw.rot = rot_all{i};
    hydros{i} = flowfield_from_saved_state(Z_all{i}, pw, t_all(i), nolog);
end
t_force = toc(t0);

planes = cell(nf, 1);
t0 = tic;
parfor i = 1:nf
    pw = p_const.Value;
    planes{i} = flowfield_grid_velocity(hydros{i}, pw, fopts, nolog, true);
end
t_vel = toc(t0);

streams = cell(nf, 1);
t0 = tic;
parfor i = 1:nf
    streams{i} = compute_three_plane_streams(planes{i}, hydros{i});
end
t_str = toc(t0);
slim_mib = NaN;
end


function [t_force, t_vel, t_str, slim_mib] = bench_flowfield_gpu_threadpool(res, par, frame_idx, fopts_gpu, nolog)
nf = numel(frame_idx);
Z_all = cell(nf, 1);
t_all = zeros(nf, 1);
rot_all = cell(nf, 1);
for i = 1:nf
    [Z_all{i}, pi_tmp, t_all(i)] = unpack_flow_frame(res, par, frame_idx(i));
    rot_all{i} = pi_tmp.rot;
end
p_nopre = par;
p_nopre.body_precomp = [];
p_const = parallel.pool.Constant(@() rebuild_precomp_on_worker(p_nopre));
parfor i = 1:min(nf, 8)
    pw = p_const.Value;
    pw.rot = rot_all{1};
    [~] = flowfield_from_saved_state(Z_all{1}, pw, t_all(1), nolog);
end

% Force extraction in parallel on ThreadPool
hydros = cell(nf, 1);
t0 = tic;
parfor i = 1:nf
    pw = p_const.Value;
    pw.rot = rot_all{i};
    hydros{i} = flowfield_from_saved_state(Z_all{i}, pw, t_all(i), nolog);
end
t_force = toc(t0);

% Batched 3-plane FP64 GPU velocity on main host (avoids GPU context contention)
planes = cell(nf, 1);
wait(gpuDevice());
t0 = tic;
for i = 1:nf
    planes{i} = flowfield_grid_velocity(hydros{i}, par, fopts_gpu, nolog, true);
end
wait(gpuDevice());
t_vel = toc(t0);

% 2-D streamlines in parallel on ThreadPool
streams = cell(nf, 1);
t0 = tic;
parfor i = 1:nf
    streams{i} = compute_three_plane_streams(planes{i}, hydros{i});
end
t_str = toc(t0);
slim_mib = NaN;
end


function st = compute_three_plane_streams(pl, hydro)
orig = hydro.X_cm * hydro.L_to_um;
sp = 0.6;
ln = 4.0;
st = struct();
st.xz = flowfield_plane_streamlines(pl.xz.h_vec, pl.xz.v_vec, pl.xz.Ux, pl.xz.Uz, ...
    pl.xz.mask, [orig(1), orig(3)], sp, ln, 'XZ', true);
st.yz = flowfield_plane_streamlines(pl.yz.h_vec, pl.yz.v_vec, pl.yz.Uy, pl.yz.Uz, ...
    pl.yz.mask, [orig(2), orig(3)], sp, ln, 'YZ', true);
st.xy = flowfield_plane_streamlines(pl.xy.h_vec, pl.xy.v_vec, pl.xy.Ux, pl.xy.Uy, ...
    pl.xy.mask, [orig(1), orig(2)], sp, ln, 'XY', true);
end


function s = tern(flag, a, b)
if flag, s = a; else, s = b; end
end
