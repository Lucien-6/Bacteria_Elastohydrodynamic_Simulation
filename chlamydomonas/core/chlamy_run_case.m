function [results, params, output_dir] = chlamy_run_case(config_file, ov)
%CHLAMY_RUN_CASE  Run one biflagellate Chlamydomonas simulation case
%
%  Description:
%    Loads the tracked Chlamydomonas JSON configuration, injects the
%    measured-gait tabulated active moment (defaulting to
%    chlamydomonas/data/fit_1p00A.mat when ov.moment_table is omitted),
%    applies optional 3-D beat-plane tilt and per-flagellum symmetry-breaking
%    overrides, and invokes run_simulation in headless mode.
%
%  Syntax:
%    [results, params, output_dir] = chlamy_run_case()
%    [results, params, output_dir] = chlamy_run_case(config_file)
%    [results, params, output_dir] = chlamy_run_case(config_file, ov)
%
%  Inputs:
%    config_file - (Optional) JSON path relative to the project root
%                  (default 'configs/config_chlamy_biflagellate_gz.json')
%    ov          - (Optional) struct with any of the fields:
%                    .case_name         output folder name
%                    .T_obs             observation time [s]
%                    .dt                output stride [s]
%                    .seg_per_um        segments per um
%                    .hydro_radius_factor  eps_pili / r_pili (default sqrt(2))
%                    .rel_tol, .abs_tol ODE tolerances
%                    .tlim              wall-clock limit of one ODE segment [s]
%                    .animation         true/false (default false)
%                    .preview           true/false (default false)
%                    .silent            true/false (default false)
%                    .moment_table      table from chlamy_fit_moment /
%                                       chlamy_make_perp_envelope (default:
%                                       chlamydomonas/data/fit_1p00A.mat)
%                    .tilt_deg          [g1 g2] beat-plane tilt [deg]:
%                                       right-handed rotation about each
%                                       flagellum's outward root tangent
%                                       ([+17, +17] = C2-chiral Wan tilt)
%                    .drive_amp         per-flagellum scale for q_par^osc (1)
%                    .drive_phase       per-flagellum beat phase shift [rad] (0)
%                    .drive_freq_ratio  per-flagellum frequency ratio (1)
%                    .drive_amp_perp    per-flagellum scale for q_perp (1)
%                    .perp_sign         per-flagellum chirality sign for q_perp
%                                       (default [+1, -1], applied once here)
%
%  Outputs:
%    results, params, output_dir - as returned by run_simulation
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-06
%  Updated: 2026-11-24  remove obsolete traveling-wave/40-deg overrides,
%                       restore cwd via onCleanup, load tracked data/fit_1p00A

arguments
    config_file (1,:) char = 'configs/config_chlamy_biflagellate_gz.json'
    ov          (1,1) struct = struct()
end

core_dir  = fileparts(mfilename('fullpath'));
chlamy_dir = fileparts(core_dir);
proj_root = fileparts(chlamy_dir);
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));
addpath(core_dir);

orig_dir = pwd;
cleanup_cd = onCleanup(@() cd(orig_dir));
cd(proj_root);   % run_simulation writes to <pwd>/output/<case_name>

cfg_path = config_file;
if ~isfile(cfg_path), cfg_path = fullfile(proj_root, config_file); end
cfg = load_sim_config(cfg_path);

if isfield(ov, 'case_name'),           cfg.case_name           = ov.case_name;           end
if isfield(ov, 'T_obs'),               cfg.T_obs               = ov.T_obs;               end
if isfield(ov, 'tlim'),                cfg.tlim                = ov.tlim;                end
if isfield(ov, 'dt'),                  cfg.dt                  = ov.dt;                  end
if isfield(ov, 'seg_per_um'),          cfg.seg_per_um          = ov.seg_per_um;          end
if isfield(ov, 'hydro_radius_factor'), cfg.hydro_radius_factor = ov.hydro_radius_factor; end
if isfield(ov, 'rel_tol'),             cfg.RelTol              = ov.rel_tol;             end
if isfield(ov, 'abs_tol'),             cfg.AbsTol              = ov.abs_tol;             end
cfg.NT = round(cfg.T_obs / cfg.dt) + 1;

n_pili = numel(cfg.pili_ext);
tilt_rad = zeros(1, n_pili);
if isfield(ov, 'tilt_deg') && ~isempty(ov.tilt_deg)
    % Beat-plane tilt: rigid rotation of each flagellum about its outward root
    % tangent (custom_pili_building 'deflection_z'); attachment and root tangent
    % remain unchanged.  Same sign [+beta, +beta] gives C2 chiral tilt.
    td = deg2rad(ov.tilt_deg(:).');
    for k = 1:n_pili
        tilt_rad(k) = pick(td, k);
        if iscell(cfg.custom_geometry.pili)
            cfg.custom_geometry.pili{k}.deflection_z = tilt_rad(k);
        else
            cfg.custom_geometry.pili(k).deflection_z = tilt_rad(k);
        end
    end
end

% Load moment table (explicit override or tracked chlamydomonas/data/fit_1p00A.mat)
if isfield(ov, 'moment_table') && ~isempty(ov.moment_table)
    tab = ov.moment_table;
else
    default_mat = fullfile(chlamy_dir, 'data', 'fit_1p00A.mat');
    assert(isfile(default_mat), 'chlamy_run_case:MissingTable', ...
        'Default moment table not found: %s', default_mat);
    S_tab = load(default_mat, 'fit');
    tab = S_tab.fit;
end

p0 = chlamy_build_params(cfg, @(varargin) []);   % scales & root frames only
assert(abs(p0.scales.T_0 - tab.T_0) < 1e-12 * tab.T_0, 'chlamy_run_case:MomentTable', ...
    'Time scale of the moment table differs from this configuration.');

for k = 1:n_pili
    fr = pick(getf(ov, 'drive_freq_ratio', 1), k);
    t0 = p0.d_local{k}(:);
    g  = tilt_rad(k);
    n0 = [0; 0; 1] * cos(g) + cross(t0, [0; 0; 1]) * sin(g) + ...
         t0 * dot(t0, [0; 0; 1]) * (1 - cos(g));
    amp_perp = pick(getf(ov, 'drive_amp_perp', 1), k) * ...
               pick(getf(ov, 'perp_sign', [1, -1]), k);
    dopt = struct('normal_local', p0.D_local{k}.' * n0, ...
        'amp',        pick(getf(ov, 'drive_amp', 1), k), ...
        'amp_perp',   amp_perp, ...
        'phase',      pick(getf(ov, 'drive_phase', 0), k), ...
        'freq_ratio', fr);
    newam = struct('type', 'custom', ...
        'omega_Hz', tab.f_Hz * fr, ...
        'Q_peak_Nm', tab.Q_peak_Nm(min(k, numel(tab.Q_peak_Nm))), ...
        'fun', chlamy_table_moment_fun(tab, k, dopt));
    cfg.pili_ext(k).active_moment = newam;
end

cfg.auto_confirm       = true;
cfg.silent             = isfield(ov, 'silent') && logical(ov.silent);
cfg.preview_initial    = isfield(ov, 'preview') && logical(ov.preview);
cfg.compute_flowfield  = false;
cfg.generate_animation = isfield(ov, 'animation') && logical(ov.animation);

[results, params, output_dir] = run_simulation(cfg);
end

function v = pick(x, k)
v = x(min(k, numel(x)));
end

function v = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
