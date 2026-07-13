function tests = test_smoke_baselines
%TEST_SMOKE_BASELINES  Five-family end-to-end baseline acceptance cases
%
%  This suite complements the focused invariant tests with short complete
%  runs for the standard, motor-driven, custom-geometry, and sphere cases.
%  The seeded Brownian family is covered by test_runtime_contracts.
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


function testSphereResistanceBaseline(testCase)
N = 120;
X = sphere_particle_sample_cached(N);
distance = sqrt(sum((reshape(X, N, 1, 3) - reshape(X, 1, N, 3)).^2, 3));
distance(1:N+1:end) = Inf;
h_mean = mean(min(distance, [], 2));

params = struct('N_body', N, 'M_pili', 0, 'N_segments', zeros(1,0), ...
    'eps_body', h_mean / 2.6667, 'epsilon', 1e-3);
[mobility, twist] = build_global_resistance(X, params, struct());
velocity = repmat([1; 0; 0], N, 1);
force = mobility \ velocity;
drag = sum(force(1:3:end));

verifyEmpty(testCase, twist);
verifyGreaterThan(testCase, drag, 0);
% build_global_resistance includes the 1/(8*pi) Stokeslet prefactor, so
% the unit-sphere drag in the project's force scaling is 6*pi/(8*pi)=3/4.
verifyEqual(testCase, drag, 0.75, 'RelTol', 0.05);
end


function testStandardPiliBaseline(testCase)
cfg = load_sim_config(fullfile(testCase.TestData.root, ...
    'configs', 'config_example.json'));
cfg = common_smoke_cfg(cfg, 'standard');
cfg.N_body = 100;
cfg.seg_per_um = 2;
cleanup = onCleanup(@() remove_smoke_output(testCase.TestData.root, cfg.case_name)); %#ok<NASGU>

results = run_simulation(cfg);

verifyTrue(testCase, results.completed);
verifyTrue(testCase, all(isfinite(results.X_cm), 'all'));
verifyTrue(testCase, all(isfinite(results.q), 'all'));
end


function testMotorBaseline(testCase)
cfg = load_sim_config(fullfile(testCase.TestData.root, ...
    'configs', 'config_motor_flagella_example.json'));
cfg = common_smoke_cfg(cfg, 'motor');
cfg.N_body = 100;
cfg.seg_per_um = 2;
cleanup = onCleanup(@() remove_smoke_output(testCase.TestData.root, cfg.case_name)); %#ok<NASGU>

results = run_simulation(cfg);

verifyTrue(testCase, results.completed);
verifyTrue(testCase, all(isfinite(results.X_cm), 'all'));
verifyTrue(testCase, all(isfinite(results.q), 'all'));
end


function testCustomGeometryBaseline(testCase)
cfg = load_sim_config(fullfile(testCase.TestData.root, ...
    'configs', 'config_custom_example.json'));
cfg = common_smoke_cfg(cfg, 'custom');
cfg.custom_geometry.body_points_file = fullfile(testCase.TestData.root, ...
    cfg.custom_geometry.body_points_file);
for k = 1:numel(cfg.custom_geometry.pili)
    cfg.custom_geometry.pili(k).points_file = ...
        fullfile(testCase.TestData.root, ...
        cfg.custom_geometry.pili(k).points_file);
end
cleanup = onCleanup(@() remove_smoke_output(testCase.TestData.root, cfg.case_name)); %#ok<NASGU>

results = run_simulation(cfg);

verifyTrue(testCase, results.completed);
verifyTrue(testCase, all(isfinite(results.X_cm), 'all'));
verifyGreaterThan(testCase, numel(results.X_pili), 0);
end


function cfg = common_smoke_cfg(cfg, label)
[~, token] = fileparts(tempname);
cfg.case_name = fullfile('__tests__', ['smoke_', label, '_', token]);
cfg.T_obs = 1e-4;
cfg.dt = 1e-4;
cfg.tlim = 120;
cfg.enable_brown = false;
cfg.preview_initial = false;
cfg.generate_animation = false;
cfg.compute_flowfield = false;
cfg.auto_confirm = true;
cfg.silent = true;
cfg.AbsTol = 1e-7;
cfg.RelTol = 1e-5;
end


function ext = empty_pili_ext()
ext = struct('intrinsic_curvature', {}, 'active_moment', {}, ...
    'root_coupling', {}, 'motor_torque', {});
end


function remove_smoke_output(root, case_name)
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
