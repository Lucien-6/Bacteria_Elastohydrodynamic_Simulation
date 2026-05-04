%% test_motor_reaction_torque.m
%  Diagnostic test: motor reaction-torque application audit (v7.3.0+)
%
%  Description:
%    In free_torsion mode the constant motor exerts a torque tau along
%    the pilus root tangent d3 (lab frame).  Historical revisions of
%    ode_coupled_rhs.m attempted to add a body-side reaction term to
%    R_global(4:6):
%        v7.0.x : R_global(4:6) -= tau_eff * J_euler(3,:)'
%        v7.2.x : R_global(4:6) -= tau_eff * d3_lab_seg1
%
%    Both choices were physically incorrect.  Rows 4-6 of R_global are
%    the FULL SYSTEM angular-momentum balance (B_global rows 4-6 sum
%    r x f over body + every pilus node, see build_coupled_lin_sys.m).
%    Motor torque between body and pilus is INTERNAL: by Newton's third
%    law the action +tau on the pilus and the reaction -tau on the
%    body cancel in the system sum, and nothing should be added to
%    R_global(4:6).  Adding a non-zero term injects a fictitious
%    external moment that over-prescribes the body's rotation.
%
%    Production code in v7.3.0+ therefore adds NEITHER candidate
%    vector; the body's reaction motion propagates implicitly through
%    the global mobility solve via the pilus moment-balance row that
%    holds the +tau term and the root kinematic constraints.
%
%    Purpose of this script:
%      1. Document the geometric discrepancy between the two legacy
%         candidates (J_euler(3,:)' vs d3_lab_seg1) at a non-trivial
%         configuration — the discrepancy is non-zero whenever the
%         pilus tangent leaves the body XY-plane.
%      2. Confirm that under zero gravity, single free_torsion pilus,
%         and constant motor torque, the body angular velocity
%         observed at t = 0 is consistent with the new "no reaction
%         in row 4-6" rule (Omega is driven purely by hydrodynamic
%         coupling through the pilus moment balance).
%
%    The procedure:
%      1. Build a single linear free_torsion pilus deflected so that its
%         root tangent has a non-trivial 3D direction.
%      2. Disable gravity by setting rho_bac = rho_fluid.
%      3. Initialise params + Z0 directly (no run_simulation).
%      4. Reconstruct geometry to obtain d3_lab(seg 1) at t = 0.
%      5. Reproduce the J_euler(3,:) row exactly as build_coupled_lin_sys.m
%         exports it (params field geom.J_euler3{k}).
%      6. Run a single ode_coupled_rhs eval and inspect the body angular
%         velocity Omega.
%      7. Print the two legacy candidate vectors, their separation,
%         the body Omega, and a sanity hint on the magnitude.
%
%    NOTE: this test does NOT assert that Omega matches any particular
%    closed-form expression — anisotropic body drag and pilus
%    hydrodynamic feedback couple non-trivially.  It is a forensic
%    print-only probe; treat the output as a record, not a pass/fail.
%
%    See also CHANGELOG.md (v7.3.0) and ode_coupled_rhs.m §"NOTE (v7.3.0)".
%
%  Usage:
%    Run from the project root or from inside tests/.  No CLI arguments.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-05-04
%  Updated: 2026-05-04 — Reframe around v7.3.0 fix: production no
%                         longer applies any reaction term to rows 4-6.

close all force;
clear; clc;

%% ---- Locate project root and add it to the path ----
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);    % parent of tests/
if ~isfile(fullfile(projectRoot, 'run_simulation.m'))
    % Fallback: assume current pwd is project root
    projectRoot = pwd;
    if ~isfile(fullfile(projectRoot, 'run_simulation.m'))
        error('test_motor_reaction_torque:WrongDirectory', ...
            ['Cannot locate run_simulation.m. Place this script in tests/ ', ...
             'under Coupled_Body_Pili, or cd to the project root before running.']);
    end
end
addpath(projectRoot);

fprintf('=========================================================\n');
fprintf(' Diagnostic: motor reaction torque (v7.3.0 audit)\n');
fprintf(' Production: NO reaction term added to R_global(4:6).\n');
fprintf(' This script only inspects legacy candidate vectors.\n');
fprintf('=========================================================\n\n');

%% ---- Construct a free_torsion + motor configuration ----
%
%  The pilus is linear, attached at site 1 (capsule -X pole).  Deflection
%  angles tilt the root tangent away from the principal axes so that
%  d3_lab and J_euler(3,:)' are unambiguously distinct.
pili_config(1).type         = 'linear';
pili_config(1).length_um    = 3;
pili_config(1).deflection_z = pi/3;     % rotate about lab Z (tilt)
pili_config(1).deflection_x = pi/4;     % rotate about lab X (twist about tangent)
pili_config(1).site         = 1;
pili_config(1).intrinsic_curvature = 'frozen';
pili_config(1).active_moment       = struct('type', 'off');
pili_config(1).root_coupling       = 'free_torsion';
pili_config(1).motor_torque        = struct( ...
    'type', 'constant', 'torque_pN_nm', 2000, 'max_omega_Hz', 300);

[Pili_Matrix, pili_ext] = parse_pili_config(pili_config);

% Root_coupling / motor_torque cell arrays for init_coupled_params
root_coupling_cfg = {pili_ext.root_coupling};
motor_torque_cfg  = {pili_ext.motor_torque};
intrinsic_cfg     = {pili_ext.intrinsic_curvature};
active_cfg        = {pili_ext.active_moment};

%% ---- Initialise parameters and state (no logging) ----
silent = @(varargin) [];   % no-op logger to keep stdout clean
params = init_coupled_params( ...
    'geometry_mode', 'parametric', ...
    'N_body',        100, ...                    % small for fast initialisation
    'minor_axis',    0.4e-6, ...
    'major_axis',    1.25e-6, ...
    'alpha_0',       0, ...
    'mu',            1e-3, ...
    'E',             3.279e8, ...
    'r_pili',        3e-9, ...
    'sigma',         0.2, ...
    'rho_bac',       998.232, ...                % == rho_fluid -> zero gravity
    'rho_fluid',     998.232, ...
    'Pili_Matrix',   Pili_Matrix, ...
    'seg_per_um',    10, ...
    'intrinsic_curvature', intrinsic_cfg, ...
    'active_moment',       active_cfg, ...
    'root_coupling',       root_coupling_cfg, ...
    'motor_torque',        motor_torque_cfg, ...
    'logf',                silent, ...
    'logf_file',           silent);

% Sanity check: gravity should be zero
fprintf('[Setup] Eg (elastogravitational number) = %.3e (should be 0)\n', params.Eg);
if abs(params.Eg) > 1e-12
    warning('test_motor_reaction_torque:GravityNotZero', ...
        'Eg is not exactly zero (%.3e); gravity will contaminate the diagnostic.', ...
        params.Eg);
end

Z0 = init_coupled_state(params, silent);

%% ---- Reconstruct geometry at t = 0 ----
geom = reconstruct_geometry(Z0, params);

% Pilus 1: lab-frame d3 of segment 1 (Newton's third-law direction)
d3_lab_seg1 = geom.d3{1}(1, :)';                  % [3 x 1]
d3_lab_seg1 = d3_lab_seg1 / norm(d3_lab_seg1);    % defensive renormalise

% Reproduce J_euler(3,:) exactly as build_coupled_lin_sys.m does
rot_k    = params.rot{1};
theta1   = geom.pili_angles{1}.theta(1);
phi1     = geom.pili_angles{1}.phi(1);
sT = sin(theta1);  cT = cos(theta1);
sP = sin(phi1);    cP = cos(phi1);
J_euler_basis_row3 = [cP/sT, sP/sT, 0];           % 1 x 3
J_euler_row3       = J_euler_basis_row3 * rot_k;  % 1 x 3 (lab frame)
J_euler3_col       = J_euler_row3.';              % [3 x 1]

% Normalised candidate reaction-torque directions (sign omitted for clarity)
v_correct = d3_lab_seg1;                              % already unit
v_v7_path = J_euler3_col / max(norm(J_euler3_col), 1e-30);

% Discrepancy metrics
diff_vec = v_correct - v_v7_path;
diff_norm = norm(diff_vec);
cos_angle = max(-1, min(1, dot(v_correct, v_v7_path)));
angle_deg = acos(cos_angle) * 180 / pi;

%% ---- Report ----
fprintf('\n[Geometry] First-segment Euler angles:\n');
fprintf('   theta_1 = %.6f rad (%.2f deg)\n', theta1, theta1*180/pi);
fprintf('   phi_1   = %.6f rad (%.2f deg)\n', phi1,   phi1  *180/pi);

fprintf('\n[Legacy A] d3_lab(seg 1)         = [%+.6f, %+.6f, %+.6f]   (v7.2.x path)\n', ...
    v_correct(1), v_correct(2), v_correct(3));
fprintf('[Legacy B] J_euler(3,:)'' (norm.) = [%+.6f, %+.6f, %+.6f]   (v7.0.x path)\n', ...
    v_v7_path(1), v_v7_path(2), v_v7_path(3));

fprintf('\n[Discrepancy]  ||A - B||_2  = %.6e\n', diff_norm);
fprintf('               angle(A, B) = %.4f deg\n', angle_deg);
fprintf('               z-component: A_z = %+.6f,  B_z = %+.6f  (B_z is structurally 0)\n', ...
    v_correct(3), v_v7_path(3));
fprintf('\n[v7.3.0 rule] Neither A nor B is added to R_global(4:6).\n');
fprintf('              Body reaction is propagated through hydrodynamic coupling.\n');

%% ---- Single ode_coupled_rhs evaluation (zero gravity) ----
%  Show the body angular acceleration that the current implementation
%  produces and compare its direction with -A and -B.  Anisotropic body
%  drag prevents Omega from being exactly anti-parallel to either, but
%  the closer alignment indicates which reaction-torque rule is in use.
dZ = ode_coupled_rhs(0, Z0, params);
Omega_body = dZ(4:6);                              % [3 x 1] body angular velocity (lab)
Omega_norm = norm(Omega_body);

if Omega_norm < 1e-30
    fprintf('\n[ode_coupled_rhs] |Omega_body| = %.3e (no body rotation observed).\n', Omega_norm);
else
    Omega_hat = Omega_body / Omega_norm;
    cos_align_A = -dot(Omega_hat, v_correct);     % >0 if Omega || -A
    cos_align_B = -dot(Omega_hat, v_v7_path);     % >0 if Omega || -B
    fprintf('\n[ode_coupled_rhs] |Omega_body| = %.3e (dimensionless rate)\n', Omega_norm);
    fprintf('                  Omega_hat   = [%+.6f, %+.6f, %+.6f]\n', ...
        Omega_hat(1), Omega_hat(2), Omega_hat(3));
    fprintf('                  cos(<Omega, -A>)  = %+.6f   (legacy v7.2 candidate)\n', cos_align_A);
    fprintf('                  cos(<Omega, -B>)  = %+.6f   (legacy v7.0 candidate)\n',  cos_align_B);
    fprintf('                  Note: in v7.3.0+ the body Omega is NOT expected to be\n');
    fprintf('                  parallel to either -A or -B; it is determined by the\n');
    fprintf('                  full hydrodynamic coupling and the anisotropic body drag.\n');
end

fprintf('\n=========================================================\n');
fprintf(' Diagnostic complete.\n');
fprintf('=========================================================\n');
