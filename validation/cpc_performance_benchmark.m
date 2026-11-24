function summary = cpc_performance_benchmark(varargin)
%CPC_PERFORMANCE_BENCHMARK  Operator & end-to-end acceleration benchmark for BacEHD
%
%  Description:
%    Measures the wall-clock speedup and numerical equivalence of BacEHD
%    v10.0.0 solver accelerations relative to the legacy dense LU / default
%    finite-difference Jacobian baseline:
%      Part A (Operator & Jacobian scaling across N_body and M_pili):
%        - Dense LU RHS (body_precomp disabled) vs body-frame Block Schur RHS
%        - Full finite-difference Jacobian vs 3-colour frozen-hydrodynamics
%          Jacobian (ode_coupled_jacobian)
%        - Exact numerical equivalence: max |dZ_schur - dZ_dense| / max(|dZ_dense|)
%      Part B (End-to-end integration speedup & trajectory equivalence):
%        - Deterministic settling (N_body=300, 1 pilus, 0.5 s):
%          dense+FD vs Block Schur + frozen_hydro Jacobian
%        - Brownian motion (N_body=300, 1 pilus, 0.2 s, dt_B=0.05 s):
%          dense+FD (no reuse) vs Block Schur + frozen_hydro + LaggedJacobian
%
%  Syntax:
%    summary = cpc_performance_benchmark()
%    summary = cpc_performance_benchmark('OutDir', outdir)
%
%  Outputs:
%    summary - struct with .operators (table), .e2e (table), .passed,
%              .matlab_version, .computer, .elapsed_s.
%              Tables are written to:
%                <OutDir>/performance_benchmark_operators.csv
%                <OutDir>/performance_benchmark_e2e.csv
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

ip = inputParser;
addParameter(ip, 'OutDir', '');
parse(ip, varargin{:});

t_all = tic;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(genpath(fullfile(root, 'src')));
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

nolog = @(varargin) [];

fprintf('================================================================\n');
fprintf('  BacEHD v%s — Performance & Equivalence Benchmark\n', project_version());
fprintf('  Platform: MATLAB %s on %s\n', version, computer);
fprintf('================================================================\n');

%% ---- Part A: Operator & Jacobian Scaling ----
configs_a = {
    'bare_N120',   120, [];
    'bare_N300',   300, [];
    'pili1_N300',  300, 2;
    'pili4_N300',  300, [2, 4, 8, 10];
    'pili4_N600',  600, [2, 4, 8, 10];
};

n_cfg = size(configs_a, 1);
op_rows = cell(n_cfg, 1);

for i = 1:n_cfg
    label  = configs_a{i, 1};
    n_body = configs_a{i, 2};
    sites  = configs_a{i, 3};

    [p_schur, Z0] = build_case_params(n_body, sites, nolog);
    p_dense = p_schur;
    p_dense.body_precomp = [];

    % Warm up both paths and check exact algebraic equivalence
    dZ_dense = ode_coupled_rhs(0, Z0, p_dense);
    dZ_schur = ode_coupled_rhs(0, Z0, p_schur);
    rel_diff = max(abs(dZ_schur - dZ_dense)) / max(max(abs(dZ_dense)), eps);

    % Repeat count keeps timing accurate without excessive runtime
    n_rep_rhs = 5;
    if numel(sites) >= 4 && n_body >= 600
        n_rep_rhs = 3;
    end
    t0 = tic;
    for r = 1:n_rep_rhs
        ode_coupled_rhs(0, Z0, p_dense);
    end
    rhs_dense_ms = 1e3 * toc(t0) / n_rep_rhs;

    t0 = tic;
    for r = 1:n_rep_rhs
        ode_coupled_rhs(0, Z0, p_schur);
    end
    rhs_schur_ms = 1e3 * toc(t0) / n_rep_rhs;
    rhs_speedup  = rhs_dense_ms / rhs_schur_ms;

    % Jacobian benchmark: FD Jacobian vs 3-colour frozen-hydro Jacobian
    % For high-DOF cases, time a 10-column sample of FD and scale to n_ode
    % so the benchmark finishes within seconds while remaining exact.
    n_ode = numel(Z0);
    t0 = tic;
    J_fast = ode_coupled_jacobian(0, Z0, p_schur); %#ok<NASGU>
    jac_fast_ms = 1e3 * toc(t0);

    n_sample = min(n_ode, 10);
    t0 = tic;
    for j = 1:n_sample
        Zj = Z0;
        Zj(j) = Zj(j) + 1e-7 * max(abs(Z0(j)), 1);
        ode_coupled_rhs(0, Zj, p_dense);
    end
    jac_fd_est_ms = (1e3 * toc(t0) / n_sample) * (n_ode + 1);
    jac_speedup   = jac_fd_est_ms / jac_fast_ms;

    n_seg_tot = sum(p_schur.N_segments);
    op_rows{i} = table(string(label), p_schur.N_body, p_schur.M_pili, n_seg_tot, n_ode, ...
        rhs_dense_ms, rhs_schur_ms, rhs_speedup, rel_diff, ...
        jac_fd_est_ms, jac_fast_ms, jac_speedup, ...
        'VariableNames', {'case_id', 'N_body_actual', 'M_pili', 'N_seg_total', 'n_ode', ...
        'rhs_dense_ms', 'rhs_schur_ms', 'rhs_speedup', 'rhs_rel_diff', ...
        'jac_fd_est_ms', 'jac_fast_ms', 'jac_speedup'});

    fprintf(['  %-11s (Nb=%3d, M=%d, DOF=%3d) | RHS: %6.2f -> %5.2f ms (%5.2fx, rel=%.1e) | ', ...
             'Jac: %7.1f -> %5.1f ms (%6.1fx)\n'], ...
        label, p_schur.N_body, p_schur.M_pili, n_ode, ...
        rhs_dense_ms, rhs_schur_ms, rhs_speedup, rel_diff, ...
        jac_fd_est_ms, jac_fast_ms, jac_speedup);
end
T_ops = vertcat(op_rows{:});

%% ---- Part B: End-to-End Integration Speedup & Trajectory Equivalence ----
fprintf('\n  --- End-to-end integration speedup & trajectory equivalence ---\n');
e2e_rows = cell(2, 1);

% B1: Deterministic settling (N_body=300, 1 pilus at site 2, T_obs=0.5 s)
[p_fast, Z0] = build_case_params(300, 2, nolog);
p_leg = p_fast;
p_leg.body_precomp = [];
T_det_star = 0.5 / p_fast.scales.T_0;
opts_base = odeset('AbsTol', 1e-4, 'RelTol', 1e-4, 'Stats', 'off');

t0 = tic;
sol_leg = ode15s(@(t, Z) ode_coupled_rhs(t, Z, p_leg), [0, T_det_star], Z0, opts_base);
wall_leg = toc(t0);

opts_fast = odeset(opts_base, 'Jacobian', @(t, Z) ode_coupled_jacobian(t, Z, p_fast));
t0 = tic;
sol_fast = ode15s(@(t, Z) ode_coupled_rhs(t, Z, p_fast), [0, T_det_star], Z0, opts_fast);
wall_fast = toc(t0);

pos_diff_det = norm(sol_fast.y(1:3, end) - sol_leg.y(1:3, end)) / ...
    max(norm(sol_leg.y(1:3, end) - Z0(1:3)), eps);
e2e_rows{1} = table("det_pili1_0p5s", p_fast.N_body, p_fast.M_pili, 0.5, ...
    wall_leg, wall_fast, wall_leg / wall_fast, pos_diff_det, ...
    'VariableNames', {'case_id', 'N_body_actual', 'M_pili', 'T_obs_s', ...
    'wall_legacy_s', 'wall_v10_s', 'speedup', 'end_pos_rel_diff'});
fprintf('  det_pili1_0p5s  | legacy %.2f s -> v10 %.2f s (%.2fx) | end pos rel diff = %.2e\n', ...
    wall_leg, wall_fast, wall_leg / wall_fast, pos_diff_det);

% B2: Brownian motion (N_body=300, 1 pilus at site 2, T_obs=0.2 s, dt_B=0.05 s)
dt_B_s = 0.05;
n_steps_b = 4;
dt_B_star = dt_B_s / p_fast.scales.T_0;

% Legacy: dense + default FD Jacobian at every Brownian segment
rng_leg = RandStream('mt19937ar', 'Seed', 12345);
Z_l = Z0;
t0 = tic;
for s = 1:n_steps_b
    p_leg.brown_force_dimless = sample_brownian_body_force(Z_l, p_leg, dt_B_s, rng_leg);
    sol_s = ode15s(@(t, Z) ode_coupled_rhs(t, Z, p_leg), [0, dt_B_star], Z_l, opts_base);
    Z_l = sol_s.y(:, end);
    Z_l(4:7) = Z_l(4:7) / norm(Z_l(4:7));
    Z_l = enforce_clamped_geometry(Z_l, p_leg);
end
wall_b_leg = toc(t0);

% v10.0.0: Block Schur + 3-colour frozen-hydro Jacobian + LaggedJacobian
rng_v10 = RandStream('mt19937ar', 'Seed', 12345);
jc = LaggedJacobian(50);
Z_f = Z0;
t0 = tic;
for s = 1:n_steps_b
    p_fast.brown_force_dimless = sample_brownian_body_force(Z_f, p_fast, dt_B_s, rng_v10);
    rhs_f = @(t, Z) ode_coupled_rhs(t, Z, p_fast);
    jac_f = @(t, Z) ode_coupled_jacobian(t, Z, p_fast);
    jc.tick();
    jc.beginSegment(rhs_f, jac_f);
    opts_b = odeset(opts_base, 'Jacobian', @(t, Z) jc.evaluate(t, Z));
    sol_s = ode15s(rhs_f, [0, dt_B_star], Z_f, opts_b);
    Z_f = sol_s.y(:, end);
    Z_f(4:7) = Z_f(4:7) / norm(Z_f(4:7));
    Z_f = enforce_clamped_geometry(Z_f, p_fast);
end
wall_b_v10 = toc(t0);

pos_diff_b = norm(Z_f(1:3) - Z_l(1:3)) / max(norm(Z_l(1:3) - Z0(1:3)), eps);
e2e_rows{2} = table("brown_pili1_0p2s", p_fast.N_body, p_fast.M_pili, n_steps_b * dt_B_s, ...
    wall_b_leg, wall_b_v10, wall_b_leg / wall_b_v10, pos_diff_b, ...
    'VariableNames', {'case_id', 'N_body_actual', 'M_pili', 'T_obs_s', ...
    'wall_legacy_s', 'wall_v10_s', 'speedup', 'end_pos_rel_diff'});
fprintf('  brown_pili1_0p2s| legacy %.2f s -> v10 %.2f s (%.2fx) | end pos rel diff = %.2e\n', ...
    wall_b_leg, wall_b_v10, wall_b_leg / wall_b_v10, pos_diff_b);

T_e2e = vertcat(e2e_rows{:});

op_csv  = fullfile(outdir, 'performance_benchmark_operators.csv');
e2e_csv = fullfile(outdir, 'performance_benchmark_e2e.csv');
writetable(T_ops, op_csv);
writetable(T_e2e, e2e_csv);

passed = all(T_ops.rhs_rel_diff < 1e-8) && all(T_e2e.end_pos_rel_diff < 1e-3);
elapsed_s = toc(t_all);

summary = struct( ...
    'operators',      T_ops, ...
    'e2e',            T_e2e, ...
    'passed',         passed, ...
    'matlab_version', version, ...
    'computer',       computer, ...
    'elapsed_s',      elapsed_s);

fprintf('----------------------------------------------------------------\n');
fprintf('Benchmark %s in %.2f s.\n  Tables: %s\n          %s\n', ...
    pass_word(passed), elapsed_s, op_csv, e2e_csv);
end


function [params, Z0] = build_case_params(n_body, sites, nolog)
if isempty(sites)
    Pili_Matrix = zeros(8, 0);
    n_seg = zeros(1, 0);
    rc = {};
    mt = {};
else
    M = numel(sites);
    pili_cfg = repmat(struct('type', 'linear', 'length_um', 3, 'site', 1), 1, M);
    for k = 1:M
        pili_cfg(k).site = sites(k);
    end
    [Pili_Matrix, ext] = parse_pili_config(pili_cfg);
    n_seg = repmat(30, 1, M);
    rc = {ext.root_coupling};
    mt = {ext.motor_torque};
end
params = init_coupled_params( ...
    'geometry_mode',       'parametric', ...
    'minor_axis',          0.4e-6, ...
    'major_axis',          1.25e-6, ...
    'N_body',              n_body, ...
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
    'N_segments',          n_seg, ...
    'Pili_Matrix',         Pili_Matrix, ...
    'root_coupling',       rc, ...
    'motor_torque',        mt, ...
    'logf',                nolog, ...
    'logf_file',           nolog);
Z0 = init_coupled_state(params, nolog);
[Z0, params] = adaptive_basis_check(Z0, params);
end


function word = pass_word(flag)
if flag
    word = 'PASSED';
else
    word = 'FAILED';
end
end
