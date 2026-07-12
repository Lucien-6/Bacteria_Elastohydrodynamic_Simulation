%% verify_migration.m
%  Migration Verification Script
%
%  Description:
%    Verifies that the directory restructuring is complete and correct.
%    Checks:
%      1. All source files are in their expected locations
%      2. Root directory only contains main entry points
%      3. All subdirectories have the expected files
%      4. No orphaned files remain in root
%
%  Usage:
%    verify_migration
%
%  Working Directory:
%    Run this script with the project root as the current working directory.
%    All checked paths are resolved relative to pwd; the script does not
%    change the current working directory.
%
%  Side Effects:
%    Executes close all force, clear, and clc before verification.  This
%    closes all figures, clears the base workspace, and clears the Command
%    Window.  It then prints the verification report to standard output.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-05-15
%  Updated: 2026-07-12

close all force;
clear; clc;

fprintf('========================================\n');
fprintf('  Migration Verification Script v1.0\n');
fprintf('  Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:ss'));
fprintf('========================================\n\n');

%% Define expected file locations
expected_structure = struct();

% Root directory: only main entry points
expected_structure.root = {
    'main_coupled_simulation.m'
    'run_simulation.m'
    'run_simulation_HPC_Linux.m'
    };

% Core simulation engine
expected_structure.core = {
    'build_coupled_lin_sys.m'
    'build_global_resistance.m'
    'ode_coupled_rhs.m'
    };

% Geometry building
expected_structure.geometry = {
    'capsule_body_building.m'
    'capsule_sphere_cap_points.m'
    'capsule_sphere_mesh_cache.m'
    'custom_body_building.m'
    'custom_pili_building.m'
    'pili_building.m'
    'pili_type_building.m'
    'read_geometry_csv_mat.m'
    'reconstruct_geometry.m'
    'snap_sphere_cap_count.m'
    'sphere_particle_sample_cached.m'
    'validate_custom_geometry.m'
    };

% Physics computations
expected_structure.physics = {
    'compute_curvature_binormal.m'
    'compute_effective_body_brownian.m'
    'compute_node_directors.m'
    'compute_pili_elastic_moments.m'
    'coupled_active_moment.m'
    'coupled_background_flow.m'
    };

% Mathematical utilities
expected_structure.math = {
    'adaptive_basis_check.m'
    'enforce_clamped_geometry.m'
    'pili_directors.m'
    'pili_euler_angles.m'
    'pili_spatial_coords.m'
    'quat_kinematics.m'
    'quat_multiply.m'
    'quat_to_rotmat.m'
    'rotation_between_vectors.m'
    'skew_matrix.m'
    };

% Configuration management
expected_structure.config = {
    'init_coupled_params.m'
    'init_coupled_state.m'
    'load_sim_config.m'
    'normalise_json_struct_array.m'
    'parse_initial_orientation.m'
    'parse_pili_config.m'
    'save_sim_config.m'
    'validate_simulation_params.m'
    };

% Visualization and export
expected_structure.visualization = {
    'compute_plot_flowfield.m'
    'export_figure_png_pdf_html.m'
    'render_progress_bar.m'
    };

% General utilities
expected_structure.utils = {
    'brown_rng_seed_from_clock.m'
    'odetpbar.m'
    'pack_state_rates.m'
    'project_version.m'
    'unpack_state.m'
    'verify_migration.m'
    };

%% Check each location
all_pass = true;

% Check root directory
fprintf('[1/8] Checking root directory ...\n');
root_files = dir('*.m');
root_names = {root_files.name}';
for i = 1:length(expected_structure.root)
    if ~ismember(expected_structure.root{i}, root_names)
        fprintf('  MISSING: %s\n', expected_structure.root{i});
        all_pass = false;
    end
end
% Check for unexpected files in root
for i = 1:length(root_names)
    if ~ismember(root_names{i}, expected_structure.root)
        fprintf('  UNEXPECTED: %s (should be in src/)\n', root_names{i});
        all_pass = false;
    end
end
if all_pass
    fprintf('  PASS - Root directory contains only expected files\n');
end

% Check src/core
fprintf('[2/8] Checking src/core/ ...\n');
if ~verify_subdir('src/core', expected_structure.core)
    all_pass = false;
end

% Check src/geometry
fprintf('[3/8] Checking src/geometry/ ...\n');
if ~verify_subdir('src/geometry', expected_structure.geometry)
    all_pass = false;
end

% Check src/physics
fprintf('[4/8] Checking src/physics/ ...\n');
if ~verify_subdir('src/physics', expected_structure.physics)
    all_pass = false;
end

% Check src/math
fprintf('[5/8] Checking src/math/ ...\n');
if ~verify_subdir('src/math', expected_structure.math)
    all_pass = false;
end

% Check src/config
fprintf('[6/8] Checking src/config/ ...\n');
if ~verify_subdir('src/config', expected_structure.config)
    all_pass = false;
end

% Check src/visualization
fprintf('[7/8] Checking src/visualization/ ...\n');
if ~verify_subdir('src/visualization', expected_structure.visualization)
    all_pass = false;
end

% Check src/utils
fprintf('[8/8] Checking src/utils/ ...\n');
if ~verify_subdir('src/utils', expected_structure.utils)
    all_pass = false;
end

%% Summary
fprintf('\n========================================\n');
if all_pass
    fprintf('  ALL CHECKS PASSED\n');
    fprintf('  Directory restructuring is complete!\n');
else
    fprintf('  SOME CHECKS FAILED\n');
    fprintf('  Please review the messages above.\n');
end
fprintf('========================================\n');

%% Helper function
function pass = verify_subdir(dir_path, expected_files)
    pass = true;
    if ~isfolder(dir_path)
        fprintf('  MISSING DIRECTORY: %s\n', dir_path);
        pass = false;
        return;
    end
    
    actual_files = dir(fullfile(dir_path, '*.m'));
    actual_names = {actual_files.name}';
    
    % Check expected files exist
    for i = 1:length(expected_files)
        if ~ismember(expected_files{i}, actual_names)
            fprintf('  MISSING: %s/%s\n', dir_path, expected_files{i});
            pass = false;
        end
    end
    
    % Check for unexpected files
    for i = 1:length(actual_names)
        if ~ismember(actual_names{i}, expected_files)
            fprintf('  UNEXPECTED: %s/%s\n', dir_path, actual_names{i});
            pass = false;
        end
    end
    
    if pass
        fprintf('  PASS - All %d files accounted for\n', length(expected_files));
    end
end
