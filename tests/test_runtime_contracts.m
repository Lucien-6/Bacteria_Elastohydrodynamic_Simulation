function tests = test_runtime_contracts
%TEST_RUNTIME_CONTRACTS  Runtime regression tests for repaired contracts
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.root = root;
addpath(root);
addpath(genpath(fullfile(root, 'src')));
end


function testStrictShearOrthogonality(testCase)
cfg = load_sim_config(fullfile(testCase.TestData.root, ...
    'configs', 'config_example.json'));
cfg.background_flow = struct( ...
    'type', 'shear_z', ...
    'U_max', 1e-6, ...
    'direction', [1, 0, 0], ...
    'shear_axis', [0, 0, 1], ...
    'shear_height_m', 5e-6);
validate_simulation_params(cfg);

cfg.background_flow.direction = [1, 0, 1];
verifyError(testCase, @() validate_simulation_params(cfg), ...
    'validate_params:InvalidParams');
end


function testUniformFlowIgnoresInvalidShearAxis(testCase)
cfg = load_sim_config(fullfile(testCase.TestData.root, ...
    'configs', 'config_example.json'));
cfg.background_flow = struct( ...
    'type', 'uniform', ...
    'U_max', 1e-6, ...
    'direction', [1, 0, 0], ...
    'shear_axis', [NaN, 0, 1]);
validate_simulation_params(cfg);

params = init_coupled_params( ...
    'background_flow', cfg.background_flow, ...
    'Pili_Matrix', zeros(8, 0), ...
    'N_segments', zeros(1, 0));

verifyEqual(testCase, params.bg_flow.type, 'uniform');
verifyEqual(testCase, params.bg_flow.shear_axis, [0; 0; 1]);
end


function testStrictShearRejectsInvalidVectors(testCase)
cfg = load_sim_config(fullfile(testCase.TestData.root, ...
    'configs', 'config_example.json'));
cfg.background_flow = struct( ...
    'type', 'shear_z', 'U_max', 1e-6, ...
    'direction', [1, 0, 0], 'shear_axis', [0, 0, 1], ...
    'shear_height_m', 5e-6);

bad = {
    rmfield(cfg.background_flow, 'direction')
    setfield(cfg.background_flow, 'direction', [0, 0, 0]) %#ok<SFLD>
    setfield(cfg.background_flow, 'direction', [NaN, 0, 0]) %#ok<SFLD>
    setfield(cfg.background_flow, 'shear_axis', [Inf, 0, 1]) %#ok<SFLD>
    setfield(cfg.background_flow, 'shear_axis', [1, 0, 0]) %#ok<SFLD>
    setfield(cfg.background_flow, 'shear_axis', [1, 1, 0]) %#ok<SFLD>
    };
for i = 1:numel(bad)
    cfg_i = cfg;
    cfg_i.background_flow = bad{i};
    verifyError(testCase, @() validate_simulation_params(cfg_i), ...
        'validate_params:InvalidParams');
end

cfg.background_flow.direction = [1, 1e-14, 0];
cfg.background_flow.shear_axis = [0, 1, 0];
validate_simulation_params(cfg);
end


function testSimulationResultClassification(testCase)
complete = struct('completed', true, 'termination_reason', 'completed');
partial = struct('completed', false, 'termination_reason', 'timeout');

[status, tf] = simulation_result_status(complete);
verifyEqual(testCase, status, "OK");
verifyTrue(testCase, tf);

[status, tf] = simulation_result_status(partial);
verifyEqual(testCase, status, "TIMEOUT");
verifyFalse(testCase, tf);
end


function testBrownianRemainderStep(testCase)
cfg = minimal_runtime_cfg(testCase.TestData.root);
cfg.case_name = unique_test_case('brownian_remainder');
cfg.enable_brown = true;
cfg.Temper = 293.15;
cfg.brown_rng_seed = 24680;
cfg.T_obs = 0.02;
cfg.dt = 0.01;
cfg.dt_B = 0.015;
cleanup = onCleanup(@() remove_test_output(testCase.TestData.root, cfg.case_name)); %#ok<NASGU>

original_rng = rng;
rng_cleanup = onCleanup(@() rng(original_rng)); %#ok<NASGU>
rng(13579, 'twister');
global_rng_before = rng;
[results, params] = run_simulation(cfg);
global_rng_after = rng;

verifyTrue(testCase, results.completed);
verifyEqual(testCase, results.termination_reason, 'completed');
verifyEqual(testCase, results.brownian_step_sizes_s, [0.015, 0.005], ...
    'AbsTol', 1e-12);
verifyEqual(testCase, sum(results.brownian_step_sizes_s), cfg.T_obs, ...
    'AbsTol', 1e-12);
verifyTrue(testCase, isfield(results, 'diagnostics'));
verifyGreaterThan(testCase, results.diagnostics.linear_solve_count, 0);
verifyLessThanOrEqual(testCase, ...
    results.diagnostics.max_backward_residual, ...
    results.diagnostics.backward_residual_error_threshold);
verifyGreaterThan(testCase, results.diagnostics.min_rcond, 0);

Z_final = [results.X_cm(:,end); results.q(:,end)];
[factor, diffusion, mobility] = ...
    compute_effective_body_brownian(Z_final, params);
verifyLessThan(testCase, norm(factor*factor' - 2*diffusion, 'fro') / ...
    norm(2*diffusion, 'fro'), 1e-12);
verifyGreaterThan(testCase, min(eig(0.5*(mobility + mobility'))), 0);
verifyEqual(testCase, global_rng_after, global_rng_before);
end


function testTimeoutReturnsPartialResults(testCase)
cfg = minimal_runtime_cfg(testCase.TestData.root);
cfg.case_name = unique_test_case('timeout');
cfg.T_obs = 1.0;
cfg.dt = 0.1;
cfg.tlim = 1e-6;
cfg.generate_animation = true;
cfg.compute_flowfield = true;
cleanup = onCleanup(@() remove_test_output(testCase.TestData.root, cfg.case_name)); %#ok<NASGU>

[results, ~, output_dir] = run_simulation(cfg);

verifyFalse(testCase, results.completed);
verifyEqual(testCase, results.termination_reason, 'timeout');
verifyLessThan(testCase, results.achieved_time_s, results.target_time_s);
verifyTrue(testCase, isfile(fullfile(output_dir, 'simulation_results.mat')));
verifyTrue(testCase, isfile(fullfile(output_dir, 'timeseries.csv')));
verifyEqual(testCase, width(results.X_cm), numel(results.t));
verifyEqual(testCase, width(results.q), numel(results.t));
verifyTrue(testCase, all(isfinite(results.t)));
verifyTrue(testCase, all(isfinite(results.X_cm), 'all'));
verifyTrue(testCase, all(isfinite(results.q), 'all'));
if numel(results.t) > 1
    verifyGreaterThan(testCase, min(diff(results.t)), 0);
end

log_text = fileread(fullfile(output_dir, 'simulation.log'));
verifyTrue(testCase, contains(log_text, 'SIMULATION PARTIAL: INTEGRATION TIMEOUT'));
verifyTrue(testCase, contains(log_text, '[ANIM] Skipped (integration timeout)'));
verifyTrue(testCase, contains(log_text, '[FLOW] Skipped (integration timeout)'));
verifyTrue(testCase, contains(log_text, 'warn 1.0e-10'));
end


function testHpcTimeoutReturnsPartialResults(testCase)
cfg = minimal_runtime_cfg(testCase.TestData.root);
cfg.case_name = unique_test_case('timeout_hpc');
cfg.T_obs = 1.0;
cfg.dt = 0.1;
cfg.tlim = 1e-6;
cleanup = onCleanup(@() remove_test_output(testCase.TestData.root, cfg.case_name)); %#ok<NASGU>

[results, ~, output_dir] = run_simulation_HPC_Linux(cfg);

verifyFalse(testCase, results.completed);
verifyEqual(testCase, results.termination_reason, 'timeout');
verifyTrue(testCase, isfile(fullfile(output_dir, 'simulation_results.mat')));
verifyTrue(testCase, isfile(fullfile(output_dir, 'timeseries.csv')));
verifyEqual(testCase, results.segment_count, 1);
end


function testDualEntryDeterministicAgreement(testCase)
cfg = minimal_runtime_cfg(testCase.TestData.root);
cfg.T_obs = 1e-4;
cfg.dt = 1e-4;
cfg.tlim = 120;
cfg.case_name = unique_test_case('dual_windows');
cleanup_main = onCleanup(@() remove_test_output( ...
    testCase.TestData.root, cfg.case_name)); %#ok<NASGU>
main_results = run_simulation(cfg);

cfg.case_name = unique_test_case('dual_hpc');
cleanup_hpc = onCleanup(@() remove_test_output( ...
    testCase.TestData.root, cfg.case_name)); %#ok<NASGU>
hpc_results = run_simulation_HPC_Linux(cfg);

verifyTrue(testCase, main_results.completed);
verifyTrue(testCase, hpc_results.completed);
verifyEqual(testCase, hpc_results.t_dimless, main_results.t_dimless, ...
    'AbsTol', 1e-14);
verifyEqual(testCase, hpc_results.X_cm, main_results.X_cm, ...
    'AbsTol', 1e-12);
verifyEqual(testCase, hpc_results.q, main_results.q, 'AbsTol', 1e-12);
end


function cfg = minimal_runtime_cfg(root)
cfg = load_sim_config(fullfile(root, 'configs', 'config_example.json'));
cfg.geometry_mode = 'parametric';
cfg.N_body = 100;
cfg.Pili_Matrix = zeros(8, 0);
cfg.pili_ext = struct('intrinsic_curvature', {}, 'active_moment', {}, ...
    'root_coupling', {}, 'motor_torque', {});
cfg.N_segments = zeros(1, 0);
cfg.preview_initial = false;
cfg.generate_animation = false;
cfg.compute_flowfield = false;
cfg.auto_confirm = true;
cfg.silent = true;
cfg.enable_brown = false;
cfg.tlim = 60;
cfg.AbsTol = 1e-7;
cfg.RelTol = 1e-5;
end


function case_name = unique_test_case(prefix)
[~, token] = fileparts(tempname);
case_name = fullfile('__tests__', [prefix, '_', token]);
end


function remove_test_output(root, case_name)
folder = fullfile(root, 'output', case_name);
if isfolder(folder)
    rmdir(folder, 's');
end
test_root = fullfile(root, 'output', '__tests__');
if isfolder(test_root)
    entries = dir(test_root);
    if numel(entries) <= 2
        rmdir(test_root);
    end
end
end
