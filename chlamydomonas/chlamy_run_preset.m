function s = chlamy_run_preset(tag, n_threads, opts)
%CHLAMY_RUN_PRESET  Authoritative preset runner for Chlamydomonas swimming cases
%
%  Description:
%    Runs one tracked Chlamydomonas reinhardtii (CC-124 / CC-125, cell
%    1p00 A at 57.75 Hz) simulation preset using the measured-gait dual-k
%    active moment (chlamydomonas/data/fit_1p00A.mat) and optional 3-D
%    out-of-plane beating (1w + 2w proximal envelope + 17 deg chiral
%    beat-plane tilt), then runs chlamy_swim_analysis and saves
%    swim_analysis.mat inside the case output folder.
%
%  Presets (TAG, case-insensitive):
%    3-D beating (tilt=[+17,+17] deg, r_perp=0.25, r_perp2=0.05,
%                 s_w=0.20, delta=delta2=+pi/2, perp_sign=[+1,-1], tol=1e-4):
%      'SYM3D'        C2-symmetric 3-D breaststroke (straight + ~1.17 Hz spin),
%                     3 beats
%      'HELIX3D'      5% cis/trans amp asymmetry (left-handed superhelix,
%                     ~1.09 Hz spin, R~1.5 um, P~-70 um), 3 beats
%      'HELIX3D_20'   Same as HELIX3D, 20 beats (~0.35 s) with 3-D MP4 animation
%      'HELIX3D_2S'   Realistic long 3-D helical run: T_obs = 2.0 s (~116 beats,
%                     ~2.2 full helical turns), dt = 1 ms, tol = 1e-4, tlim = 120 h,
%                     with 3-D MP4 animation -> output/Chlamy-3D-Helix-2s/
%
%    Planar reference / control presets:
%      'PLANAR_SYM50' (alias 'SYM50')  Coplanar symmetric breaststroke, 50 beats
%      'PLANAR_A75F75'(alias 'A75F75') Coplanar flagellum-2 amp x0.75 & freq x0.75, 1 s
%
%  Syntax:
%    s = chlamy_run_preset(tag)
%    s = chlamy_run_preset(tag, n_threads)
%    s = chlamy_run_preset(tag, n_threads, opts)
%
%  Inputs:
%    tag       - Preset name listed above
%    n_threads - (Optional) BLAS thread cap for this MATLAB session
%    opts      - (Optional) struct:
%                  .dry_run  true: validate geometry, drive handles and initial
%                            ODE RHS without integrating (default false)
%
%  Output:
%    s - Struct from chlamy_swim_analysis (or dry-run diagnostic struct)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

arguments
    tag       (1,:) char
    n_threads       = []
    opts      (1,1) struct = struct()
end

if ~isempty(n_threads)
    maxNumCompThreads(n_threads);
end

chlamy_dir = fileparts(mfilename('fullpath'));
core_dir   = fullfile(chlamy_dir, 'core');
proj_root  = fileparts(chlamy_dir);
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));
addpath(chlamy_dir);
addpath(core_dir);

orig_dir = pwd;
cleanup_cd = onCleanup(@() cd(orig_dir));
w_prev = warning('off', 'all');
cleanup_w = onCleanup(@() warning(w_prev));
cd(proj_root);

data_file = fullfile(chlamy_dir, 'data', 'fit_1p00A.mat');
assert(isfile(data_file), 'chlamy_run_preset:MissingData', ...
    'Moment table not found: %s', data_file);
S_fit = load(data_file, 'fit');
tab_base = S_fit.fit;

[ov, n_skip, descr] = build_preset_overrides(tag, tab_base);
cfg_file = 'configs/config_chlamy_biflagellate_gz.json';

if isfield(opts, 'dry_run') && logical(opts.dry_run)
    s = run_dry_check(cfg_file, ov, tag, descr);
    return
end

fprintf('[CHLAMY] Preset %-12s | %s\n', upper(tag), descr);
[res, ~, out_dir] = chlamy_run_case(cfg_file, ov);
s = chlamy_swim_analysis(res, struct('f_Hz', tab_base.f_Hz, 'n_skip', n_skip));
s.tag = upper(tag);
s.description = descr;
save(fullfile(out_dir, 'swim_analysis.mat'), 's');
end


function [ov, n_skip, descr] = build_preset_overrides(tag, tab_base)
%BUILD_PRESET_OVERRIDES  Map a preset tag to chlamy_run_case overrides
key = upper(strtrim(tag));
ov = struct('moment_table', tab_base);

% 3-D out-of-plane moment table (1w + 2w proximal envelope, C2 chirality)
tab_3d = chlamy_make_perp_envelope(tab_base, struct( ...
    'r_perp', 0.25, 'r_perp2', 0.05, 's_w', 0.20, ...
    'delta_perp', pi/2, 'delta_perp2', pi/2, 'quiet', true));

switch key
    case {'SYM3D', 'OOP3D_SYM'}
        ov.case_name      = 'Chlamy-3D-Sym';
        ov.moment_table   = tab_3d;
        ov.tilt_deg       = [17, 17];
        ov.perp_sign      = [1, -1];
        ov.drive_amp      = [1, 1];
        ov.drive_amp_perp = [1, 1];
        ov.T_obs          = 3 / tab_base.f_Hz;
        ov.rel_tol        = 1e-4;
        ov.abs_tol        = 1e-4;
        n_skip            = 1;
        descr             = '3-D symmetric breaststroke (3 beats, tol=1e-4)';

    case {'HELIX3D', 'OOP3D_HELIX'}
        ov.case_name      = 'Chlamy-3D-Helix';
        ov.moment_table   = tab_3d;
        ov.tilt_deg       = [17, 17];
        ov.perp_sign      = [1, -1];
        ov.drive_amp      = [1, 0.95];
        ov.drive_amp_perp = [1, 1];
        ov.T_obs          = 3 / tab_base.f_Hz;
        ov.rel_tol        = 1e-4;
        ov.abs_tol        = 1e-4;
        n_skip            = 1;
        descr             = '3-D left-handed superhelix, 5% asym (3 beats, tol=1e-4)';

    case {'HELIX3D_20', 'OOP3D_HELIX20'}
        ov.case_name      = 'Chlamy-3D-Helix-20b';
        ov.moment_table   = tab_3d;
        ov.tilt_deg       = [17, 17];
        ov.perp_sign      = [1, -1];
        ov.drive_amp      = [1, 0.95];
        ov.drive_amp_perp = [1, 1];
        ov.T_obs          = 20 / tab_base.f_Hz;
        ov.dt             = 0.001;
        ov.rel_tol        = 1e-4;
        ov.abs_tol        = 1e-4;
        ov.animation      = true;
        ov.tlim           = 120 * 3600;
        n_skip            = 2;
        descr             = '3-D helical swimming, 5% asym (20 beats, animation, tol=1e-4)';

    case {'HELIX3D_2S', 'OOP3D_HELIX2S'}
        ov.case_name      = 'Chlamy-3D-Helix-2s';
        ov.moment_table   = tab_3d;
        ov.tilt_deg       = [17, 17];
        ov.perp_sign      = [1, -1];
        ov.drive_amp      = [1, 0.95];
        ov.drive_amp_perp = [1, 1];
        ov.T_obs          = 2.0;
        ov.dt             = 0.001;
        ov.rel_tol        = 1e-4;
        ov.abs_tol        = 1e-4;
        ov.animation      = true;
        ov.tlim           = 120 * 3600;
        n_skip            = 2;
        descr             = '3-D realistic helical swimming (2.0 s ~ 116 beats ~ 2.2 turns, animation, tol=1e-4)';

    case {'PLANAR_SYM50', 'SYM50'}
        ov.case_name      = 'Chlamy-SB-SYM50';
        ov.T_obs          = 50 / tab_base.f_Hz;
        ov.animation      = true;
        n_skip            = 2;
        descr             = 'Coplanar symmetric reference (50 beats, animation)';

    case {'PLANAR_A75F75', 'A75F75'}
        ov.case_name        = 'Chlamy-SB-A75F75';
        ov.drive_amp        = [1, 0.75];
        ov.drive_freq_ratio = [1, 0.75];
        ov.T_obs            = 1.0;
        ov.dt               = 0.001;
        ov.tlim             = 120 * 3600;
        n_skip              = 2;
        descr               = 'Coplanar flagellum-2 amp x0.75 & freq x0.75 (1.0 s)';

    otherwise
        error('chlamy_run_preset:UnknownTag', ...
            ['Unknown preset ''%s''. Valid tags: SYM3D, HELIX3D, ', ...
             'HELIX3D_20, HELIX3D_2S, PLANAR_SYM50, PLANAR_A75F75.'], tag);
end
end


function diag = run_dry_check(cfg_file, ov, tag, descr)
%RUN_DRY_CHECK  Verify config, 3-D geometry, moment handles and RHS without ODE
cfg = load_sim_config(cfg_file);
cfg.case_name = ov.case_name;
cfg.T_obs     = ov.T_obs;
if isfield(ov, 'dt'),      cfg.dt     = ov.dt;      end
if isfield(ov, 'rel_tol'), cfg.RelTol = ov.rel_tol; end
if isfield(ov, 'abs_tol'), cfg.AbsTol = ov.abs_tol; end
if isfield(ov, 'tlim'),    cfg.tlim   = ov.tlim;    end
cfg.NT = round(cfg.T_obs / cfg.dt) + 1;

n_pili = numel(cfg.pili_ext);
tilt_rad = zeros(1, n_pili);
if isfield(ov, 'tilt_deg')
    td = deg2rad(ov.tilt_deg(:).');
    for k = 1:n_pili
        tilt_rad(k) = td(min(k, numel(td)));
        cfg.custom_geometry.pili(k).deflection_z = tilt_rad(k);
    end
end

tab = ov.moment_table;
p0 = chlamy_build_params(cfg, @(varargin) []);
for k = 1:n_pili
    t0 = p0.d_local{k}(:);
    g  = tilt_rad(k);
    n0 = [0; 0; 1]*cos(g) + cross(t0, [0; 0; 1])*sin(g) + t0*dot(t0, [0; 0; 1])*(1 - cos(g));
    amp_k = 1;
    if isfield(ov, 'drive_amp'), amp_k = ov.drive_amp(min(k, numel(ov.drive_amp))); end
    ps_k = 1;
    if isfield(ov, 'perp_sign'), ps_k = ov.perp_sign(min(k, numel(ov.perp_sign))); end
    dopt = struct('normal_local', p0.D_local{k}.' * n0, 'amp', amp_k, 'amp_perp', ps_k);
    cfg.pili_ext(k).active_moment = struct('type', 'custom', ...
        'omega_Hz', tab.f_Hz, 'fun', chlamy_table_moment_fun(tab, k, dopt));
end
[params, Z0] = chlamy_build_params(cfg, @(varargin) []);
dZ0 = ode_coupled_rhs(0, Z0, params);
assert(all(isfinite(dZ0)), 'chlamy_run_preset:DryRunNonFinite', 'Initial RHS is non-finite.');

diag = struct('tag', upper(tag), 'case_name', cfg.case_name, 'description', descr, ...
    'T_obs_s', cfg.T_obs, 'dt_s', cfg.dt, 'NT', cfg.NT, ...
    'n_beats_approx', cfg.T_obs * tab.f_Hz, 'RelTol', cfg.RelTol, ...
    'N_seg', params.N_segments, 'rhs0_norm', norm(dZ0), 'ok', true);
fprintf(['[DRY-RUN PASS] %-10s -> output/%s | T_obs=%.2f s (~%.1f beats, NT=%d) | ', ...
    'tol=%.0e | ||dZ(0)||=%.3e\n'], diag.tag, diag.case_name, diag.T_obs_s, ...
    diag.n_beats_approx, diag.NT, diag.RelTol, diag.rhs0_norm);
end
