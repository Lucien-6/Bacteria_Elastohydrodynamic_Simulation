function tests = test_hydrodynamic_invariants
%TEST_HYDRODYNAMIC_INVARIANTS  Hydrodynamic reciprocity regression tests
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
src_dir = fullfile(project_root, 'src');
addpath(genpath(src_dir));
testCase.TestData.src_dir = src_dir;
end


function teardownOnce(testCase)
rmpath(genpath(testCase.TestData.src_dir));
end


function testWeightedSelfAdjointUnequalSegments(testCase)
A = [2.0, -1.0, 0.5; 3.0, 4.0, -2.0; 1.0, 0.25, 5.0];
lengths = [0.4; 1.1; 2.3];
A_corr = weighted_self_adjoint(A, lengths);
weighted_operator = diag(lengths) * A_corr;

testCase.verifyLessThan(norm(weighted_operator - weighted_operator', 'fro'), ...
    1e-13 * norm(weighted_operator, 'fro'));
end


function testWeightedSelfAdjointUniformLimit(testCase)
A = [1.0, 2.0, -3.0; 4.0, 5.0, 6.0; 7.0, -8.0, 9.0];
lengths = 1.7 * ones(3, 1);
expected = 0.5 * (A + A');

testCase.verifyEqual(weighted_self_adjoint(A, lengths), expected, ...
    'AbsTol', 2e-15);
end


function testGlobalTwistUsesActualSegmentLengths(testCase)
X_all = [0.0, 1.0, 0.0; ...
         0.0, 0.0, 0.0; ...
         1.0, 0.0, 0.0; ...
         3.0, 0.0, 0.0];
params.N_body = 1;
params.M_pili = 1;
params.N_segments = 2;
params.eps_body = 0.2;
params.epsilon = 0.1;
geom.d3 = {[1.0, 0.0, 0.0; 1.0, 0.0, 0.0]};

[~, M_twist] = build_global_resistance(X_all, params, geom);
lengths = [1.0; 2.0];
weighted_operator = diag(lengths) * M_twist;

testCase.verifyTrue(all(isfinite(M_twist(:))));
testCase.verifyLessThan(norm(weighted_operator - weighted_operator', 'fro'), ...
    1e-11 * max(1, norm(weighted_operator, 'fro')));
end


function testWeightedCrossVirtualWork(testCase)
rng(20260713, 'twister');
M_rotlet = randn(9, 6);
D_twist = randn(6, 2);
segment_lengths = [0.65; 1.4];

quadrature_sets = {[1.0; 0.5; 1.5], [0.25; 2.0; 0.75]};
for k = 1:numel(quadrature_sets)
    quad_weights = quadrature_sets{k};
    [C_forward, C_reverse, W_f_diag, W_m_diag] = ...
        assemble_twist_coupling(M_rotlet, D_twist, ...
                                quad_weights, segment_lengths);

    expected_reverse = diag(W_m_diag) \ ...
        (C_forward' * diag(W_f_diag));
    testCase.verifyEqual(C_forward, 0.5 * M_rotlet * D_twist, ...
        'AbsTol', 5e-15);
    testCase.verifyEqual(C_reverse, expected_reverse, 'AbsTol', 5e-15);

    velocity = randn(9, 1);
    torque_density = randn(2, 1);
    forward_work = velocity' * ...
        (W_f_diag .* (C_forward * torque_density));
    reverse_work = torque_density' * ...
        (W_m_diag .* (C_reverse * velocity));
    testCase.verifyEqual(forward_work, reverse_work, ...
        'AbsTol', 2e-13 * max(1, abs(forward_work)));
end
end
