function tests = test_geometry_objectivity
%TEST_GEOMETRY_OBJECTIVITY  Objective director initialization tests
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

tests = functiontests(localfunctions);

end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(project_root, 'src')));
testCase.TestData.project_root = project_root;
end


function testCurvedFrozenPilusIsObjective(testCase)
[params, reference_points] = curvedReferenceParams();
Z0 = init_coupled_state(params, @(varargin) []);
geom = reconstruct_geometry(Z0, params);

R_body = quat_to_rotmat(params.q_0_init);
expected_points = (R_body * reference_points')';
verifyEqual(testCase, geom.X_pili{1}, expected_points, 'AbsTol', 2e-11);

ref = params.directors_ref{1};
verifyEqual(testCase, geom.d1{1}, (R_body * ref.d1')', 'AbsTol', 2e-11);
verifyEqual(testCase, geom.d2{1}, (R_body * ref.d2')', 'AbsTol', 2e-11);
verifyEqual(testCase, geom.d3{1}, (R_body * ref.d3')', 'AbsTol', 2e-11);

[kappa1, kappa2, kappa3] = compute_pili_elastic_moments( ...
    geom.d1_basis{1}, geom.d2_basis{1}, geom.d3_basis{1}, ...
    params.N_segments(1), params.sigma, params.kb0{1}, ...
    params.root_coupling{1});
verifyLessThan(testCase, max(abs([kappa1; kappa2; kappa3])), 5e-10);
end


function testRootProjectionLeavesInteriorSegmentsUnchanged(testCase)
[params, ~] = curvedReferenceParams();
Z0 = init_coupled_state(params, @(varargin) []);

axis = [0.3; -0.5; 0.8];
axis = axis / norm(axis);
angle = 0.43;
Z0(4:7) = [cos(angle / 2); axis * sin(angle / 2)];
Z_projected = enforce_clamped_geometry(Z0, params, true);

N = params.N_segments(1);
interior = true(size(Z0));
interior(1:7) = false;
interior([8, 8 + N, 8 + 2 * N]) = false;
verifyEqual(testCase, Z_projected(interior), Z0(interior), 'AbsTol', 0);
end


function testIdentityStraightPilusDoesNotRegress(testCase)
N = 8;
root = [0.1, -0.2, 0.3];
segment_length = 0.07;
points = root + (0:N)' * segment_length .* [1, 0, 0];
theta = (pi / 2) * ones(N, 1);
phi = zeros(N, 1);
psi = zeros(N, 1);
[d1, d2, d3] = pili_directors(theta, phi, psi);

params = baseParams(points, d1, d2, d3, [1; 0; 0; 0]);
Z0 = init_coupled_state(params, @(varargin) []);
geom = reconstruct_geometry(Z0, params);

verifyEqual(testCase, geom.X_pili{1}, points, 'AbsTol', 2e-12);
verifyEqual(testCase, geom.d1_basis{1}, d1, 'AbsTol', 2e-12);
verifyEqual(testCase, geom.d2_basis{1}, d2, 'AbsTol', 2e-12);
verifyEqual(testCase, geom.d3_basis{1}, d3, 'AbsTol', 2e-12);
end


function [params, points] = curvedReferenceParams()
s = linspace(0, 1, 27)'.^1.3;
raw = [0.15 + 0.9 * s, ...
       -0.1 + 0.24 * sin(1.7 * pi * s), ...
       0.2 + 0.18 * (1 - cos(1.3 * pi * s))];
N = 18;
points = canonical_pili_centerline(raw, N);
tangents = diff(points, 1, 1);
tangents = tangents ./ sqrt(sum(tangents.^2, 2));
theta = acos(max(-1, min(1, tangents(:, 3))));
phi = atan2(tangents(:, 2), tangents(:, 1));
psi = zeros(N, 1);
[d1, d2, d3] = pili_directors(theta, phi, psi);

axis = [1; 2; -0.7];
axis = axis / norm(axis);
angle = 0.91;
q = [cos(angle / 2); axis * sin(angle / 2)];
params = baseParams(points, d1, d2, d3, q);
end


function params = baseParams(points, d1, d2, d3, q)
N = size(points, 1) - 1;
segment_lengths = sqrt(sum(diff(points, 1, 1).^2, 2));

params.q_0_init = q;
params.M_pili = 1;
params.N_segments = N;
params.rot = {eye(3)};
params.pili_x0 = {points};
params.r_local = {points(1, :)'};
params.ell = sum(segment_lengths);
params.body_local = zeros(3, 1);
params.directors_ref = {struct('d1', d1, 'd2', d2, 'd3', d3)};
params.D_local = {[d1(1, :)', d2(1, :)', d3(1, :)']};
params.root_coupling = {'clamped'};
params.sigma = 0.2;
params.kb0 = {referenceCurvature(d1, d2, d3, N)};
end


function kb0 = referenceCurvature(d1, d2, d3, N)
[d1_node, d2_node, d3_node] = compute_node_directors(d1, d2, d3, N);
d1_ds = [d1_node(2, :) - d1_node(1, :); ...
         d1_node(2:end, :) - d1_node(1:end-1, :)] * N;
d2_ds = [d2_node(2, :) - d2_node(1, :); ...
         d2_node(2:end, :) - d2_node(1:end-1, :)] * N;
d3_ds = [d3_node(2, :) - d3_node(1, :); ...
         d3_node(2:end, :) - d3_node(1:end-1, :)] * N;

kb0 = zeros(N + 1, 3);
kb0(:, 1) = dot(d3_node, d2_ds, 2);
kb0(:, 2) = dot(d1_node, d3_ds, 2);
kb0(:, 3) = dot(d2_node, d1_ds, 2);
end
