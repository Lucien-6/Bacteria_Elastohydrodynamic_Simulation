function tests = test_uniform_stiffness_contract
%TEST_UNIFORM_STIFFNESS_CONTRACT  Global pilus stiffness contract tests
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


function testPiliBuildingRejectsNinthRow(testCase)
matrix9 = zeros(9, 1);
verifyError(testCase, ...
    @() pili_building(matrix9, 10, 1.25e-6, 0.4e-6, 1e-8, 3e8, 3e-9), ...
    'pili_building:InvalidPiliMatrixRows');
end


function testInitParamsRejectsNinthRowBeforeGeometry(testCase)
verifyError(testCase, ...
    @() init_coupled_params('Pili_Matrix', zeros(9, 1)), ...
    'init_coupled_params:InvalidPiliMatrixRows');
end


function testLoadConfigRejectsRawNinthRow(testCase)
json_path = [tempname, '.json'];
cleanup = onCleanup(@() deleteIfPresent(json_path));
fid = fopen(json_path, 'w', 'n', 'UTF-8');
assert(fid ~= -1, 'Could not create temporary JSON file.');
fprintf(fid, ['{"geometry_mode":"parametric","pili":', ...
    '[[1],[3],[0],[0],[0],[0],[0],[1],[1.0e-20]]}']);
fclose(fid);

verifyError(testCase, @() load_sim_config(json_path), ...
    'load_sim_config:InvalidPiliMatrixRows');
end


function testNamedParametricStiffnessIsRejected(testCase)
pili = struct('type', 'linear', 'length_um', 2, 'site', 1, ...
    'bending_stiffness', 1e-20);
verifyError(testCase, @() parse_pili_config(pili), ...
    'parse_pili_config:UnsupportedBendingStiffness');
end


function testParserRejectsNinthMatrixRow(testCase)
verifyError(testCase, @() parse_pili_config(zeros(9, 1)), ...
    'parse_pili_config:InvalidPiliMatrixRows');
end


function testCustomStiffnessIsValidationError(testCase)
custom_pilus = struct('points_file', 'missing-pilus.csv', ...
    'site_index', 1, 'bending_stiffness', 1e-20);
site = struct('position', [0; 0; 0], ...
    'direction', [1; 0; 0], 'label', 'root');
cfg = struct('body_points_file', 'missing-body.csv', ...
    'attachment_sites', site, 'pili', custom_pilus);

[errors, ~] = validate_custom_geometry(cfg, {}, {});
message = strjoin(errors, newline);
verifyTrue(testCase, contains(message, 'bending_stiffness is not supported'));
verifyTrue(testCase, contains(message, 'global material E and r_pili'));
end


function testAllPiliUseGlobalEI(testCase)
matrix = [1, 1; ...
          2, 3; ...
          0, 0; ...
          0, 0; ...
          NaN, NaN; ...
          NaN, NaN; ...
          NaN, NaN; ...
          1, 6];
E = 3.279e8;
r_pili = 3e-9;
expected_EI = E * pi * r_pili^4 / 4;

Pili = pili_building(matrix, 4, 1.25e-6, 0.4e-6, 1e-8, E, r_pili);
verifyEqual(testCase, [Pili.A], repmat(expected_EI, 1, 2), ...
    'RelTol', 10 * eps);
end


function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
