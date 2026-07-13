function tests = test_numerical_guards
%TEST_NUMERICAL_GUARDS  Regression tests for centralized solve diagnostics
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


function testHealthyLinearSolveResidual(testCase)
A = [4.0, 1.0; 2.0, 3.0];
b = [1.0; -2.0];
x = A \ b;
backward_residual = check_linear_solve(A, x, b, 'healthy test');

testCase.verifyLessThanOrEqual(backward_residual, 10 * eps);
end


function testElevatedResidualWarns(testCase)
A = eye(2);
b = ones(2, 1);
x = [1.0; 1.0 + 1e-8];

numerical_diagnostics_state('reset');
testCase.verifyWarning(@() check_linear_solve(A, x, b, 'warning test'), ...
    'numerical_diagnostics:ElevatedBackwardResidual');
snapshot = numerical_diagnostics_state('snapshot');
testCase.verifyEqual(snapshot.rcond_check_count, 1);
end


function testBadMatrixResidualErrors(testCase)
A = diag([1.0, 1e-12]);
b = [1.0; 1.0];
x_bad = [1.0; 0.0];

testCase.verifyError(@() check_linear_solve(A, x_bad, b, 'bad test'), ...
    'numerical_diagnostics:LargeBackwardResidual');
end


function testNonFiniteValueErrors(testCase)
A = eye(2);
b = ones(2, 1);
x = [1.0; NaN];

testCase.verifyError(@() check_linear_solve(A, x, b, 'nonfinite test'), ...
    'numerical_diagnostics:NonFiniteValue');
end


function testSingularConditionErrors(testCase)
A = diag([1.0, 1e-16]);
testCase.verifyError(@() check_matrix_condition(A, 'singular test'), ...
    'numerical_diagnostics:SingularMatrix');
end


function testNegativeDissipationErrors(testCase)
Mtt = diag([1.0, -0.1]);
testCase.verifyError(@() check_weighted_mobility_health( ...
    Mtt, [], [], [], ones(2,1), [], 'negative dissipation test'), ...
    'numerical_diagnostics:NegativeDissipation');
end


function testWeightedReciprocityErrors(testCase)
Mtt = eye(2);
Cf = [1.0; 0.0];
Cr = [0.0, 0.0];
Mtw = 1.0;
testCase.verifyError(@() check_weighted_mobility_health( ...
    Mtt, Cf, Cr, Mtw, ones(2,1), 1.0, 'reciprocity test'), ...
    'numerical_diagnostics:ReciprocityFailure');
end


function testMultipleMotorsConvergeToTrueResidual(testCase)
a = [0.1; -0.2];
B = [0.08, 0.01; -0.02, 0.06];
tau_motor = [1.0; 0.8];
om_max = [2.0; 1.5];

tau = solve_motor_self_consistent(a, B, tau_motor, om_max);
omega = a + B * tau;
tau_expected = tau_motor .* (1 - abs(omega) ./ om_max);
tau_expected(abs(omega) >= om_max) = 0;

true_residual = norm(tau - tau_expected, inf) / ...
    (1 + norm(tau, inf));
testCase.verifyLessThanOrEqual(true_residual, 1e-12);
end


function testMultipleMotorNonconvergenceErrors(testCase)
a = zeros(2, 1);
B = 10 * eye(2);
tau_motor = ones(2, 1);
om_max = ones(2, 1);

testCase.verifyError( ...
    @() solve_motor_self_consistent(a, B, tau_motor, om_max), ...
    'solve_motor_self_consistent:NoConvergence');
end
