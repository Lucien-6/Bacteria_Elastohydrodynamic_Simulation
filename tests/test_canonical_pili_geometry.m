function tests = test_canonical_pili_geometry
%TEST_CANONICAL_PILI_GEOMETRY  Equal-chord canonical centreline tests
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

tests = functiontests(localfunctions);

end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.project_root = project_root;
addpath(genpath(fullfile(project_root, 'src')));
end


function testCircularArcHasEqualChords(testCase)
theta_base = linspace(0, 1, 19)'.^1.6;
theta = 0.85 * pi * theta_base;
radius = 2.3e-6;
raw = [radius * cos(theta), radius * sin(theta), 0.2e-6 * theta];

[canonical, chord, parameters] = canonical_pili_centerline(raw, 48);
lengths = sqrt(sum(diff(canonical, 1, 1).^2, 2));

verifyEqual(testCase, canonical(1, :), raw(1, :), 'AbsTol', 1e-15);
verifyEqual(testCase, canonical(end, :), raw(end, :), 'AbsTol', 1e-15);
verifyGreaterThan(testCase, min(diff(parameters)), 0);
verifyEqual(testCase, mean(lengths), chord, 'RelTol', 1e-12);
verifyLessThan(testCase, relativeSpread(lengths), 1e-10);
end


function testHelixHasEqualChords(testCase)
t_base = linspace(0, 1, 31)'.^1.35;
t = 3.4 * pi * t_base;
raw = [0.7e-6 * cos(t), 0.7e-6 * sin(t), 0.18e-6 * t];

[canonical, ~, parameters] = canonical_pili_centerline(raw, 72);
lengths = sqrt(sum(diff(canonical, 1, 1).^2, 2));

verifyGreaterThan(testCase, min(diff(parameters)), 0);
verifyLessThan(testCase, relativeSpread(lengths), 1e-10);
verifyEqual(testCase, canonical([1, end], :), raw([1, end], :), ...
    'AbsTol', 1e-15);
end


function testCustomBuilderUsesCanonicalMetadata(testCase)
t = 1.6 * pi * linspace(0, 1, 16)'.^1.4;
raw = [0.55e-6 * cos(t), 0.55e-6 * sin(t), 0.24e-6 * t];
points_file = [tempname, '.csv'];
writematrix(raw, points_file);
cleanup = onCleanup(@() deleteIfPresent(points_file));

root_direction = raw(2, :) - raw(1, :);
root_direction = root_direction(:) / norm(root_direction);
site = struct('position', zeros(3, 1), ...
    'direction', root_direction, 'label', 'root');
cfg = struct('points_file', points_file, 'site_index', 1);

Pili = custom_pili_building(cfg, site, 8, 1e-6, 0, ...
    3.279e8, 3e-9, @(varargin) []);
lengths = sqrt(sum(diff(Pili.x, 1, 2).^2, 1));
root_tangent = Pili.x(:, 2) - Pili.x(:, 1);
root_tangent = root_tangent / norm(root_tangent);

verifyLessThan(testCase, relativeSpread(lengths(:)), 1e-10);
verifyEqual(testCase, Pili.ell, mean(lengths), 'RelTol', 1e-12);
verifyEqual(testCase, Pili.L, sum(lengths), 'RelTol', 1e-12);
verifyEqual(testCase, Pili.root_dir, root_tangent, 'AbsTol', 1e-12);
end


function testDegenerateCurvesAreRejected(testCase)
verifyError(testCase, ...
    @() canonical_pili_centerline(zeros(4, 3), 8), ...
    'canonical_pili_centerline:DegenerateCurve');

points_with_duplicate = [0, 0, 0; 1, 0, 0; 1, 0, 0; 2, 0, 0];
verifyError(testCase, ...
    @() canonical_pili_centerline(points_with_duplicate, 8), ...
    'canonical_pili_centerline:DegenerateInputSegment');
end


function value = relativeSpread(values)
value = (max(values) - min(values)) / mean(values);
end


function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
