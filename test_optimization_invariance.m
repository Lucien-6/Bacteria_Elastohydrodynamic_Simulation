function test_optimization_invariance()
%TEST_OPTIMIZATION_INVARIANCE  Numerical equivalence test for Stage 2/3 optimisations
%
%  Description:
%    Validates that refactored versions of build_global_resistance and
%    build_coupled_lin_sys produce numerically identical outputs to the
%    frozen reference copies in _legacy/.  Two states are tested:
%      (1) The initial state Z0 from init_coupled_state.
%      (2) A randomly perturbed state Z0 + 1e-3 * randn(size(Z0)).
%
%    All matrix comparisons use the relative Frobenius-norm criterion
%        ||diff||_F / max(||ref||_F, eps_floor)  <  1e-14
%
%  Syntax:
%    test_optimization_invariance()
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-21

root = fileparts(mfilename('fullpath'));
if isempty(root); root = pwd; end
addpath(fullfile(root, '_legacy'));

TOL      = 1e-14;        % Relative tolerance
EPS_FLR  = 1e-300;       % Floor for denominator
rng(20260421, 'twister');  % Reproducible perturbation

%% ---------- Build config (mirrors profile_short_run) ----------
cfg = local_build_cfg();

%% ---------- Build params and initial state ----------
params = init_coupled_params( ...
    'N_body',       cfg.N_body, ...
    'minor_axis',   cfg.minor_axis, ...
    'major_axis',   cfg.major_axis, ...
    'alpha_0',      cfg.alpha_0, ...
    'mu',           cfg.mu, ...
    'E',            cfg.E, ...
    'r_pili',       cfg.r_pili, ...
    'sigma',        cfg.sigma_p, ...
    'rho_bac',      cfg.rho_bac, ...
    'rho_fluid',    cfg.rho_fluid, ...
    'Pili_Matrix',  cfg.Pili_Matrix, ...
    'seg_per_um',   cfg.seg_per_um, ...
    'delta',        cfg.delta, ...
    'Temper',       293.15, ...
    'enable_brown', false, ...
    'fixman_midpoint', false, ...
    'intrinsic_curvature', cfg.intrinsic_curvature_cfg, ...
    'active_moment',       cfg.active_moment_cfg, ...
    'background_flow',     struct('type', 'off'));

Z0 = init_coupled_state(params);

%% ---------- Run test at Z0 and at a perturbed Z ----------
fprintf('\n================== test_optimization_invariance ==================\n');
fprintf('TOL (relative Frobenius) : %.0e\n', TOL);
fprintf('N_body = %d, M_pili = %d, N_segs = [%s], state DOF = %d\n\n', ...
    params.N_body, params.M_pili, num2str(params.N_segments(:)'), numel(Z0));

all_pass = true;

% --- Case 1: initial state ---
fprintf('---- Case 1: initial state Z0 ----\n');
case1 = run_comparison(Z0, params, TOL, EPS_FLR);
all_pass = all_pass && case1;

% --- Case 2: perturbed state ---
Zp = Z0 + 1e-3 * randn(size(Z0));
fprintf('\n---- Case 2: perturbed state Z0 + 1e-3 * randn ----\n');
case2 = run_comparison(Zp, params, TOL, EPS_FLR);
all_pass = all_pass && case2;

fprintf('\n=================================================================\n');
if all_pass
    fprintf('RESULT: ALL CHECKS PASSED (below %.0e)\n', TOL);
else
    fprintf('RESULT: FAILED — see relative errors above\n');
end
fprintf('=================================================================\n\n');

end

%% ------------------------------------------------------------------------
function pass = run_comparison(Z, params, TOL, EPS_FLR)

% Geometry (shared across both paths)
geom = reconstruct_geometry(Z, params);

%% Compare build_global_resistance
[M_ref, Atw_ref, Mr_ref] = build_global_resistance_v0(geom.X_all, params, geom);
[M_new, Atw_new, Mr_new] = build_global_resistance(geom.X_all, params, geom);

e_M   = rel_err(M_ref,  M_new,  EPS_FLR);
e_Mr  = rel_err(Mr_ref, Mr_new, EPS_FLR);
e_Atw = 0;
for k = 1:numel(Atw_ref)
    e_Atw = max(e_Atw, rel_err(Atw_ref{k}, Atw_new{k}, EPS_FLR));
end

report('M_weighted      ', e_M,   TOL);
report('M_rotlet        ', e_Mr,  TOL);
report('A_twist_cells   ', e_Atw, TOL);

%% Compare build_coupled_lin_sys
[ls_ref, geom_ref, hydro_ref] = build_coupled_lin_sys_v0(Z, params);
[ls_new, geom_new, hydro_new] = build_coupled_lin_sys(Z, params);

e_ls  = rel_err(ls_ref,  ls_new,  EPS_FLR);
e_Bg  = rel_err(hydro_ref.B_global,   hydro_new.B_global,   EPS_FLR);
e_Mw  = rel_err(hydro_ref.M_weighted, hydro_new.M_weighted, EPS_FLR);
if ~isempty(hydro_ref.A_twist_total)
    e_At = rel_err(hydro_ref.A_twist_total, hydro_new.A_twist_total, EPS_FLR);
else
    e_At = 0;
end
if ~isempty(hydro_ref.M_rotlet)
    e_Mrh = rel_err(hydro_ref.M_rotlet, hydro_new.M_rotlet, EPS_FLR);
else
    e_Mrh = 0;
end

report('lin_sys         ', e_ls,  TOL);
report('hydro.B_global  ', e_Bg,  TOL);
report('hydro.M_weighted', e_Mw,  TOL);
report('hydro.A_twist_tot', e_At, TOL);
report('hydro.M_rotlet  ', e_Mrh, TOL);

% Geometry should be identical (deterministic from Z)
e_X = rel_err(geom_ref.X_all, geom_new.X_all, EPS_FLR);
report('geom.X_all      ', e_X, TOL);

pass = max([e_M, e_Mr, e_Atw, e_ls, e_Bg, e_Mw, e_At, e_Mrh, e_X]) < TOL;

end

%% ------------------------------------------------------------------------
function e = rel_err(A, B, eps_flr)
%REL_ERR Relative Frobenius error.
d = norm(A(:) - B(:));
n = max(norm(A(:)), eps_flr);
e = d / n;
end

function report(name, err, tol)
if err < tol
    mark = 'PASS';
else
    mark = 'FAIL';
end
fprintf('  [%s] %s : rel_err = %.3e (tol %.0e)\n', mark, name, err, tol);
end

%% ------------------------------------------------------------------------
function cfg = local_build_cfg()
% Matches profile_short_run.m defaults.

cfg.N_body      = 300;
cfg.minor_axis  = 0.4e-6;
cfg.major_axis  = 1.25e-6;
cfg.alpha_0     = 0;

pili_config(1).type = 'linear';
pili_config(1).length_um = 3;
pili_config(1).deflection_z = pi / 2;
pili_config(1).deflection_x = pi / 2;
pili_config(1).site = 3;
pili_config(1).intrinsic_curvature = 'frozen';
pili_config(1).active_moment = struct('type', 'off');

pili_config(2).type = 'linear';
pili_config(2).length_um = 3;
pili_config(2).deflection_z = pi / 2;
pili_config(2).deflection_x = pi / 2;
pili_config(2).site = 5;
pili_config(2).intrinsic_curvature = 'frozen';
pili_config(2).active_moment = struct('type', 'off');

[cfg.Pili_Matrix, pili_ext] = parse_pili_config(pili_config);

% Mirror run_simulation.m lines 82-84 (convert struct array to cell arrays)
cfg.intrinsic_curvature_cfg = {pili_ext.intrinsic_curvature};
cfg.active_moment_cfg       = {pili_ext.active_moment};

cfg.seg_per_um = 10;
cfg.E          = 3.279e8;
cfg.r_pili     = 3e-9;
cfg.sigma_p    = 0.2;
cfg.mu         = 1e-3;
cfg.rho_bac    = 998.232;
cfg.rho_fluid  = 998.232;
cfg.delta      = pi / 20;

end
