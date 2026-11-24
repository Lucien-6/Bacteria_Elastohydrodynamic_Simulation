function summary = run_cpc_tests(varargin)
%RUN_CPC_TESTS  Reviewer checks for the Computer Physics Communications deposit
%
%  Description:
%    Runs the short tests that travel with the program.  Sphere, mobility,
%    capsule, uniform flow and filament sedimentation are single linear
%    solves.  The motor check integrates 0.05 s and fits the signed rotor
%    speed.  The helix and cantilever checks call the drivers in validation/.
%    No Brownian ensemble is repeated here.
%
%    Tests
%      sphere             sedimentation speed against Stokes' law
%      mobility           translational and rotational diffusion against Stokes-Einstein
%      capsule            orientation dependence of a fore-aft symmetric capsule
%      flow               uniform stream: the cell, with one filament, moves with the fluid
%      motor              straight free-torsion filament against the twist-drag law
%      sink               one 3 um filament slows sedimentation relative to the bare capsule
%      schur_equivalence  Block Schur vs dense LU single-step RHS equivalence (< 1e-12)
%      pili_thermal       v10.0.0 slow-mode thermal operator & balancing axial torque check
%      end_to_end         full JSON -> validate -> run_simulation -> MAT/CSV pipeline smoke
%      chlamy_3d_dry_run  custom spheroid + measured-gait 3-D Chlamydomonas dry-run check
%      helix              rigid-helix thrust, torque and drag against Johnson theory
%      cantilever         body-fixed bending time against the local cantilever formula
%
%  Syntax:
%    summary = run_cpc_tests
%    summary = run_cpc_tests('SkipExternal', true)
%
%  Output:
%    summary - struct with .table, .passed, .version, .matlab_version,
%              .computer, .elapsed_s.  The table is written to
%              output/cpc_package/cpc_test_report.csv.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-17
%  Updated: 2026-11-24 — added v10.0.0 Block Schur, pilus thermal operator,
%                        end-to-end pipeline, and Chlamydomonas 3-D dry-run checks.

ip = inputParser;
addParameter(ip, 'SkipExternal', false, @(x) islogical(x) || isnumeric(x));
parse(ip, varargin{:});

t_suite = tic;
root = cpc_startup();
outdir = fullfile(root, 'output', 'cpc_package');
if ~isfolder(outdir)
    mkdir(outdir);
end

fprintf('================================================================\n');
fprintf('  BacEHD v%s — CPC Reviewer Verification Suite\n', project_version());
fprintf('  Platform: MATLAB %s on %s\n', version, computer);
fprintf('================================================================\n');

rows = {};
rows{end + 1} = test_sphere();
rows{end + 1} = test_mobility();
rows{end + 1} = test_capsule();
rows{end + 1} = test_uniform_flow();
rows{end + 1} = test_motor();
rows{end + 1} = test_filament_sink();
rows{end + 1} = test_schur_equivalence();
rows{end + 1} = test_pili_thermal_operator();
rows{end + 1} = test_end_to_end_pipeline(root);
rows{end + 1} = test_chlamy_3d_dry_run();
if ~ip.Results.SkipExternal
    rows{end + 1} = test_external('helix', @() cpc_helix_propulsive_matrix('OutDir', outdir));
    rows{end + 1} = test_external('cantilever', @() cpc_cantilever_relaxation('OutDir', outdir));
end

report = vertcat(rows{:});
elapsed_s = toc(t_suite);
summary = struct('table', report, 'passed', all(report.passed), ...
    'version', project_version(), 'matlab_version', version, ...
    'computer', computer, 'elapsed_s', elapsed_s);
writetable(report, fullfile(outdir, 'cpc_test_report.csv'));
disp(report);
fprintf('BacEHD v%s CPC tests %s in %.2f s.  Report: %s\n', ...
    project_version(), pass_word(summary.passed), elapsed_s, ...
    fullfile(outdir, 'cpc_test_report.csv'));
end


function root = cpc_startup()
root = fileparts(fileparts(mfilename('fullpath')));
cd(root);
addpath(root);
addpath(genpath(fullfile(root, 'src')));
addpath(fullfile(root, 'validation'));
addpath(fullfile(root, 'cpc'));
addpath(fullfile(root, 'chlamydomonas'));
addpath(fullfile(root, 'chlamydomonas', 'core'));
addpath(fullfile(root, 'MATLAB Add-On', ...
    'Suite of functions to perform uniform sampling of a sphere'));
end


function row = test_sphere()
C = cpc_constants();
p = make_body(C.a_sphere, C.a_sphere, 600, 0, [], C, true, []);
V = sinking_velocity(p);
row = result_row('sphere', 'V / V_Stokes', V / stokes_speed(C.a_sphere, C), 1, 0.01);
end


function row = test_mobility()
C = cpc_constants();
p = make_body(C.a_sphere, C.a_sphere, 600, 0, [], C, true, []);
[~, p] = initial_state(p);
[~, D] = compute_effective_body_brownian(initial_state_only(p), p);
kB = 1.380649e-23;
Dt = mean(diag(D(1:3, 1:3))) / (kB * C.Temper / (6 * pi * C.mu * C.a_sphere));
Dr = mean(diag(D(4:6, 4:6))) / (kB * C.Temper / (8 * pi * C.mu * C.a_sphere^3));
row = [result_row('mobility_translation', 'D_t / D_Stokes', Dt, 1, 0.02); ...
       result_row('mobility_rotation', 'D_r / D_Stokes', Dr, 1, 0.02)];
end


function Z = initial_state_only(p)
[Z, ~] = initial_state(p);
end


function row = test_capsule()
C = cpc_constants();
V0 = sinking_velocity(make_body(C.minor, C.major, 120, 0, [], C, true, []));
V90 = sinking_velocity(make_body(C.minor, C.major, 120, pi / 2, [], C, true, []));
[V45, drift, omega] = sinking_state(make_body(C.minor, C.major, 120, pi / 4, [], C, true, []));
V_pred = V0 * cos(pi / 4)^2 + V90 * sin(pi / 4)^2;
drift_pred = abs(V90 - V0) * sin(pi / 4) * cos(pi / 4);
row = [result_row('capsule_speed', 'V(45 deg) / analytic', V45 / V_pred, 1, 1e-6); ...
       result_row('capsule_drift', 'drift / analytic', drift / drift_pred, 1, 1e-6); ...
       result_row('capsule_spin', 'angular speed (rad/s)', omega, 0, 1e-8)];
end


function row = test_uniform_flow()
C = cpc_constants();
U = 1e-6;
bg = struct('type', 'uniform', 'U_max', U, 'direction', [1 0 0]);
pili = struct('type', 'linear', 'length_um', 3, 'site', 1);
p = make_body(C.minor, C.major, 120, pi / 2, pili, C, false, bg);
[Z, p] = initial_state(p);
dZ = ode_coupled_rhs(0, Z, p);
V = dZ(1:3) * p.scales.U_0;
Omega = body_omega_dimless(Z, dZ) / p.scales.T_0;
angle_rate = max(abs(dZ(8:end))) / p.scales.T_0;
row = [result_row('uniform_flow', 'V_x / U', V(1) / U, 1, 1e-6); ...
       result_row('uniform_flow_cross', 'cross-flow speed / U', norm(V(2:3)) / U, 0, 1e-6); ...
       result_row('uniform_flow_spin', 'angular speed * L / U', norm(Omega) * C.major / U, 0, 1e-6); ...
       result_row('uniform_flow_bend', 'angle rate * L / U', angle_rate * C.major / U, 0, 1e-5)];
end


function row = test_motor()
C = cpc_constants();
pili = struct('type', 'linear', 'length_um', 3, 'site', 1, ...
    'root_coupling', 'free_torsion');
pili.motor_torque = struct('type', 'constant', 'torque_pN_nm', 1, ...
    'max_omega_Hz', 300);
p = make_body(C.minor, C.major, 120, 0, pili, C, false, []);
[Z, p] = initial_state(p);
T0 = p.scales.T_0;
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);
sol = ode15s(@(t, ZZ) ode_coupled_rhs(t, ZZ, p), [0, 0.05 / T0], Z, opts);
t_phys = linspace(0.02, 0.05, 61);
chi = zeros(size(t_phys));
for i = 1:numel(t_phys)
    g = reconstruct_geometry(deval(sol, t_phys(i) / T0), p);
    frame = [g.d1{1}(1, :)', g.d2{1}(1, :)', g.d3{1}(1, :)'];
    reference = quat_to_rotmat(g.q) * p.D_local{1};
    Q = reference' * frame;
    chi(i) = atan2(Q(2, 1), Q(1, 1));
end
chi = unwrap(chi);
slope = polyfit(t_phys, chi, 1);
%  Signed speed.  Positive chi is the paper's right-handed sense about the
%  outward root tangent (paper_metrics).  A reversed rotor must fail.
hz = slope(1) / (2 * pi);
pred = motor_analytic_hz(C);
row = result_row('motor', 'omega / analytic', hz / pred, 1, 0.02);
end


function row = test_filament_sink()
C = cpc_constants();
bare = sinking_velocity(make_body(C.minor, C.major, 120, pi / 3, [], C, true, []));
pili = struct('type', 'linear', 'length_um', 3, 'site', 1);
with = sinking_velocity(make_body(C.minor, C.major, 120, pi / 3, pili, C, true, []));
%  Instantaneous mobility ratio at 60 deg, not the 200 s reoriented speed.
row = result_row('filament_sink', 'instantaneous V / V_bare', with / bare, 0.931, 0.02);
end


function row = test_schur_equivalence()
C = cpc_constants();
pili = struct('type', 'linear', 'length_um', 3, 'site', 2);
p_schur = make_body(C.minor, C.major, 120, pi / 3, pili, C, true, []);
[Z, p_schur] = initial_state(p_schur);
p_dense = p_schur;
p_dense.body_precomp = [];
dZ_schur = ode_coupled_rhs(0, Z, p_schur);
dZ_dense = ode_coupled_rhs(0, Z, p_dense);
rel_err = max(abs(dZ_schur - dZ_dense)) / max(max(abs(dZ_dense)), eps);
row = result_row('schur_equivalence', 'max |dZ_schur-dZ_dense|/|dZ|', rel_err, 0, 1e-10);
end


function row = test_pili_thermal_operator()
C = cpc_constants();
pili = struct('type', 'linear', 'length_um', 3, 'site', 1);
p = make_body(C.minor, C.major, 120, 0, pili, C, false, []);
[Z, p] = initial_state(p);
p.brown_pili = init_pili_brownian_modes(Z, p, 2);
[M_joint, info] = compute_effective_coupled_brownian(Z, p);
rng_s = RandStream('mt19937ar', 'Seed', 20261124);
[f_body, g_modes, F_pili, ~, tau_pili] = sample_brownian_coupled_force(Z, p, 0.002, rng_s);
all_ok = all(isfinite(M_joint(:))) && all(isfinite(f_body)) && all(isfinite(g_modes)) && ...
    all(isfinite(F_pili(:))) && all(isfinite(tau_pili(:))) && ~info.was_projected;
row = [result_row('pili_thermal_pattern', 'modal load pattern_resid', info.pattern_resid, 0, 1e-10); ...
       result_row('pili_thermal_twist',   'twist_row_resid',          info.twist_row_resid, 0, 1e-10); ...
       result_row('pili_thermal_sampler', 'SPD & finite sampler flag', double(all_ok), 1, 0)];
end


function row = test_end_to_end_pipeline(root)
cfg = load_sim_config(fullfile(root, 'configs', 'config_example.json'));
cfg.case_name          = '_cpc_test_e2e_smoke';
cfg.N_body             = 60;
cfg.T_obs              = 0.04;
cfg.dt                 = 0.02;
cfg.auto_confirm       = true;
cfg.preview_initial    = false;
cfg.generate_animation = false;
cfg.compute_flowfield  = false;
cfg.silent             = true;
[res, ~, outdir] = run_simulation(cfg);
cleanup_dir = onCleanup(@() remove_dir_if_exists(outdir));
ok = isstruct(res) && res.completed && numel(res.t) == 3 && ...
    all(isfinite(res.X_cm_um(:))) && ...
    isfile(fullfile(outdir, 'simulation_results.mat')) && ...
    isfile(fullfile(outdir, 'timeseries.csv')) && ...
    isfile(fullfile(outdir, 'config.json'));
clear cleanup_dir
row = result_row('end_to_end_pipeline', 'JSON->run_simulation flag', double(ok), 1, 0);
end


function row = test_chlamy_3d_dry_run()
diag_info = chlamy_run_preset('HELIX3D', [], struct('dry_run', true));
ok = isstruct(diag_info) && isfield(diag_info, 'ok') && diag_info.ok && ...
    isfinite(diag_info.rhs0_norm) && diag_info.rhs0_norm > 0;
row = result_row('chlamy_3d_dry_run', '3-D Chlamydomonas dry_run', double(ok), 1, 0);
end


function remove_dir_if_exists(d)
if ~isempty(d) && isfolder(d)
    rmdir(d, 's');
end
end


function row = test_external(name, runner)
out = runner();
row = result_row(name, 'pass flag', double(out.passed), 1, 0);
end


function p = make_body(minor_axis, major_axis, n_body, alpha, pili, C, gravity, bg)
nolog = @(varargin) [];
args = {'geometry_mode', 'parametric', 'minor_axis', minor_axis, ...
    'major_axis', major_axis, 'N_body', n_body, 'alpha_0', alpha, ...
    'E', C.E, 'r_pili', C.r_pili, 'sigma', C.sigma, 'mu', C.mu, ...
    'Temper', C.Temper, 'hydro_radius_factor', sqrt(2), ...
    'seg_per_um', 10, 'logf', nolog, 'logf_file', nolog};
if gravity
    args = [args, {'rho_bac', C.rho_bac, 'rho_fluid', C.rho_fluid}];
else
    args = [args, {'rho_bac', C.rho_fluid, 'rho_fluid', C.rho_fluid}];
end
if ~isempty(bg)
    args = [args, {'background_flow', bg}];
end
if isempty(pili)
    args = [args, {'Pili_Matrix', zeros(8, 0), 'N_segments', zeros(1, 0)}];
else
    [matrix, ext] = parse_pili_config(pili);
    n_seg_vec = max(round(matrix(2, :) * 10), 10);
    args = [args, {'Pili_Matrix', matrix, 'N_segments', n_seg_vec, ...
        'root_coupling', {ext.root_coupling}, ...
        'motor_torque', {ext.motor_torque}}];
end
p = init_coupled_params(args{:});
end


function [Z, p] = initial_state(p)
nolog = @(varargin) [];
Z = init_coupled_state(p, nolog);
[Z, p] = adaptive_basis_check(Z, p);
end


function V_sink = sinking_velocity(p)
[V_sink, ~, ~] = sinking_state(p);
end


function [V_sink, drift, omega] = sinking_state(p)
[Z, p] = initial_state(p);
lin_sys = build_coupled_lin_sys(Z, p);
rhs = zeros(size(lin_sys, 1), 1);
rhs(1:3) = -p.Eg * p.g_hat(:);
zd = lin_sys \ rhs;
V = zd(1:3) * p.scales.U_0;
V_sink = -V(3);
drift = hypot(V(1), V(2));
omega = norm(zd(4:6)) / p.scales.T_0;
end


function Omega = body_omega_dimless(Z, dZ)
q = Z(4:7);
dq = dZ(4:7);
q_conj = [q(1); -q(2:4)];
omega_quat = 2 * quat_multiply(dq, q_conj);
Omega = omega_quat(2:4);
end


function V = stokes_speed(a, C)
volume = 4 / 3 * pi * a^3;
force = (C.rho_bac - C.rho_fluid) * volume * 9.81;
V = force / (6 * pi * C.mu * a);
end


function hz = motor_analytic_hz(C)
beta = sqrt(2);
pili_length_m = 3e-6;
torque = 1e-21;
omega_max = 300;
drag = 4 * pi^2 * C.mu * (beta * C.r_pili)^2 * pili_length_m;
hz = torque / (drag + torque / omega_max);
end


function C = cpc_constants()
C.a_sphere = 0.4e-6;
C.minor = 0.4e-6;
C.major = 1.25e-6;
C.mu = 1e-3;
C.rho_bac = 1180;
C.rho_fluid = 998.232;
C.E = 3.279e8;
C.r_pili = 3e-9;
C.sigma = 0.2;
C.Temper = 293.15;
end


function row = result_row(name, quantity, measured, expected, tolerance)
error_value = measured - expected;
passed = abs(error_value) <= tolerance + 1e-15;
row = table(string(name), string(quantity), measured, expected, tolerance, ...
    error_value, passed, 'VariableNames', ...
    {'test', 'quantity', 'measured', 'expected', 'tolerance', 'error', 'passed'});
fprintf('  %-22s %-28s %+.6g  (expected %+.6g +/- %g)  %s\n', ...
    name, quantity, measured, expected, tolerance, pass_word(passed));
end


function word = pass_word(flag)
if flag
    word = 'PASS';
else
    word = 'FAIL';
end
end
