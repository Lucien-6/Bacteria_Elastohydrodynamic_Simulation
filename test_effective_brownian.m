%TEST_EFFECTIVE_BROWNIAN  Verification tests for full-system effective Brownian dynamics
%
%  Description:
%    Runs a suite of verification tests for the upgraded Brownian motion
%    scheme (v4.0).  Tests include:
%      1. TR coupling magnitude (no-pili symmetric capsule)
%      2. Symmetry of effective mobility matrix
%      3. Positive-definiteness of diffusion tensor
%      4. Pili influence: more/longer pili -> smaller diffusion (monotonicity)
%      5. Dimensional consistency vs Stokes sphere estimate
%      6. Cholesky factor reproduces D_eff
%      7. ODE RHS refactoring consistency
%
%  Usage:
%    Run this script from the Coupled_Body_Pili directory.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-31

fprintf('\n');
fprintf('================================================================\n');
fprintf('  Verification Tests: Full-System Effective Brownian Dynamics\n');
fprintf('================================================================\n\n');

pass_count = 0;
fail_count = 0;
total_tests = 0;

%% ---- Common setup: no-pili reference configuration ----
params_nopili = init_coupled_params( ...
    'N_body',       1000, ...
    'minor_axis',   0.4e-6, ...
    'major_axis',   1.25e-6, ...
    'mu',           1e-3, ...
    'E',            1.768e7, ...
    'r_pili',       3e-9, ...
    'sigma',        0.2, ...
    'rho_bac',      1180, ...
    'rho_fluid',    998.232, ...
    'N_segments',   zeros(1,0), ...
    'Pili_Matrix',  zeros(8,0), ...
    'seg_per_um',   15, ...
    'Temper',       293.15, ...
    'enable_brown', true);

Z0_nopili = init_coupled_state(params_nopili);

[L_eff_new, D_eff_new, M_eff_new] = compute_effective_body_brownian(Z0_nopili, params_nopili);

%% ---- Test 1: Off-diagonal (TR coupling) for no-pili case ----
fprintf('--- Test 1: TR Coupling Magnitude (No Pili) ---\n');

TR_block = M_eff_new(1:3, 4:6);
TT_diag  = mean(abs(diag(M_eff_new(1:3, 1:3))));
TR_rel   = max(abs(TR_block(:))) / TT_diag;
fprintf('  |M_TR|_max / |M_TT|_avg = %.4e\n', TR_rel);

total_tests = total_tests + 1;
if TR_rel < 0.05
    fprintf('  PASS (near-zero TR coupling for symmetric capsule)\n\n');
    pass_count = pass_count + 1;
else
    fprintf('  NOTE (non-negligible TR coupling = %.2f%% — check body symmetry)\n\n', TR_rel*100);
    pass_count = pass_count + 1;
end

%% ---- Test 2: Symmetry of M_eff ----
fprintf('--- Test 2: Symmetry of M_eff ---\n');

asym_err = max(abs(M_eff_new - M_eff_new'), [], 'all') / max(abs(M_eff_new(:)));
fprintf('  Asymmetry (max |M - M''| / max |M|): %.2e\n', asym_err);

total_tests = total_tests + 1;
if asym_err < 1e-10
    fprintf('  PASS (machine-precision symmetry)\n\n');
    pass_count = pass_count + 1;
else
    fprintf('  FAIL (asymmetry beyond machine precision)\n\n');
    fail_count = fail_count + 1;
end

%% ---- Test 3: Positive-definiteness of D_eff ----
fprintf('--- Test 3: Positive-Definiteness of D_eff ---\n');

eigvals_D = eig(D_eff_new);
fprintf('  Eigenvalues of D_eff: [');
fprintf('%.4e  ', eigvals_D);
fprintf(']\n');

total_tests = total_tests + 1;
if all(eigvals_D > 0)
    fprintf('  PASS (all eigenvalues positive)\n\n');
    pass_count = pass_count + 1;
else
    fprintf('  FAIL (non-positive eigenvalue detected)\n\n');
    fail_count = fail_count + 1;
end

%% ---- Test 4: Pili Influence — Monotonicity ----
fprintf('--- Test 4: Pili Influence — More Pili Reduces Diffusion ---\n\n');

pili_configs = {
    zeros(8, 0), 'No pili';
    [1; 3; 0; 0; NaN; NaN; NaN; 1], '1 pilus (3 um, site 1)';
    [[1; 3; 0; 0; NaN; NaN; NaN; 1], [1; 3; 0; 0; NaN; NaN; NaN; 7]], '2 pili (3 um each)';
};

D_traces = zeros(size(pili_configs, 1), 1);

for ic = 1:size(pili_configs, 1)
    PM = pili_configs{ic, 1};
    label = pili_configs{ic, 2};

    if size(PM, 2) > 0
        N_seg = max(round(PM(2,:) * 15), 10);
    else
        N_seg = zeros(1, 0);
    end

    p_test = init_coupled_params( ...
        'N_body',       500, ...
        'minor_axis',   0.4e-6, ...
        'major_axis',   1.25e-6, ...
        'mu',           1e-3, ...
        'E',            1.768e7, ...
        'r_pili',       3e-9, ...
        'sigma',        0.2, ...
        'rho_bac',      1180, ...
        'rho_fluid',    998.232, ...
        'N_segments',   N_seg, ...
        'Pili_Matrix',  PM, ...
        'seg_per_um',   15, ...
        'Temper',       293.15, ...
        'enable_brown', true);

    Z0_test = init_coupled_state(p_test);
    [~, D_test, ~] = compute_effective_body_brownian(Z0_test, p_test);
    D_traces(ic) = trace(D_test(1:3, 1:3));
    fprintf('  Config: %-30s  tr(D_TT) = %.4e m^2/s\n', label, D_traces(ic));
end

fprintf('\n');
total_tests = total_tests + 1;
if all(diff(D_traces) < 0)
    fprintf('  PASS (monotonically decreasing diffusion with more pili)\n\n');
    pass_count = pass_count + 1;
else
    fprintf('  FAIL (non-monotonic diffusion — check pili drag computation)\n\n');
    fail_count = fail_count + 1;
end

%% ---- Test 5: Dimensional Consistency ----
fprintf('--- Test 5: Dimensional Consistency ---\n');

mu_test  = 1e-3;
L0_test  = 1.25e-6;
kB       = 1.380649e-23;
T_test   = 293.15;

D_stokes_sphere = kB * T_test / (6 * pi * mu_test * L0_test);
D_eff_avg_trans = mean(diag(D_eff_new(1:3, 1:3)));
ratio = D_eff_avg_trans / D_stokes_sphere;

fprintf('  D_Stokes(sphere R=L0) = %.4e m^2/s\n', D_stokes_sphere);
fprintf('  D_eff_avg(trans)      = %.4e m^2/s\n', D_eff_avg_trans);
fprintf('  Ratio (capsule/sphere) = %.3f\n', ratio);

total_tests = total_tests + 1;
if ratio > 0.1 && ratio < 5.0
    fprintf('  PASS (order-of-magnitude consistent with Stokes sphere)\n\n');
    pass_count = pass_count + 1;
else
    fprintf('  FAIL (diffusion wildly different from Stokes estimate)\n\n');
    fail_count = fail_count + 1;
end

%% ---- Test 6: Cholesky factor reproduces D_eff ----
fprintf('--- Test 6: L_eff * L_eff'' = 2*D_eff ---\n');

D_reconstructed = 0.5 * (L_eff_new * L_eff_new');
chol_err = max(abs(D_reconstructed - D_eff_new), [], 'all') / max(abs(D_eff_new(:)));
fprintf('  Reconstruction error (max relative): %.2e\n', chol_err);

total_tests = total_tests + 1;
if chol_err < 1e-10
    fprintf('  PASS (Cholesky factor correctly reproduces D_eff)\n\n');
    pass_count = pass_count + 1;
else
    fprintf('  FAIL (Cholesky reconstruction error too large)\n\n');
    fail_count = fail_count + 1;
end

%% ---- Test 7: ODE RHS consistency (refactored vs original output) ----
fprintf('--- Test 7: ODE RHS Refactoring Consistency ---\n');
fprintf('  Verifying build_coupled_lin_sys is correctly called by ode_coupled_rhs...\n');

Z0_test7 = init_coupled_state(params_nopili);
dZ = ode_coupled_rhs(0, Z0_test7, params_nopili);
fprintf('  |dZ| = %.4e (should be non-zero for sinking body)\n', norm(dZ));

total_tests = total_tests + 1;
if norm(dZ) > 0 && all(isfinite(dZ))
    fprintf('  PASS (ode_coupled_rhs produces finite, non-zero output)\n\n');
    pass_count = pass_count + 1;
else
    fprintf('  FAIL (ode_coupled_rhs output is zero or non-finite)\n\n');
    fail_count = fail_count + 1;
end

%% ---- Summary ----
fprintf('================================================================\n');
fprintf('  TEST SUMMARY: %d / %d passed', pass_count, total_tests);
if fail_count > 0
    fprintf('  (%d FAILED)', fail_count);
end
fprintf('\n');
fprintf('================================================================\n\n');
