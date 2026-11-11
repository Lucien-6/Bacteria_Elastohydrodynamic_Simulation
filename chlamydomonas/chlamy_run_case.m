function [results, params, output_dir] = chlamy_run_case(config_file, ov)
%CHLAMY_RUN_CASE  Run one Chlamydomonas case with optional parameter overrides
%
%  Description:
%    Loads a JSON configuration, applies scalar overrides, forces headless
%    execution (auto_confirm = true) and calls run_simulation.  Used by the
%    calibration and convergence drivers so every run starts from the same
%    tracked JSON file.
%
%  Syntax:
%    [results, params, output_dir] = chlamy_run_case(config_file)
%    [results, params, output_dir] = chlamy_run_case(config_file, ov)
%
%  Inputs:
%    config_file - JSON path relative to the project root
%    ov          - (Optional) struct with any of the fields
%                    .case_name         output folder name
%                    .T_obs             observation time [s]
%                    .dt                output stride [s]
%                    .seg_per_um        segments per um
%                    .hydro_radius_factor  eps_pili / r_pili of the pilus
%                                       regularisation (default sqrt(2))
%                    .root_half_angle_deg  root half-angle [deg] between each
%                                       flagellum root tangent and the cell
%                                       axis; sets attachment_sites(1:2).direction
%                                       to (cos b, +-sin b, 0).  The solver
%                                       aligns the (shape-preserving) rest
%                                       arc to it by a rigid rotation about
%                                       Z, so the beat plane stays XY.  Keep
%                                       |b - 40 deg| < 45 deg (the CSV files
%                                       are generated for 40 deg).
%                    .amplitude_Nm      active-moment amplitude (all flagella)
%                    .wavenumber        active-moment k_w (all flagella)
%                    .omega_Hz          beat frequency (all flagella)
%                    .animation         true/false (default false)
%                    .preview           true/false (default false; saves the
%                                       initial-pose figure, ~40 MB PDF)
%                    .rel_tol, .abs_tol ODE tolerances
%                    .moment            struct(.Q_Nm, .lambda_rel, .s_taper,
%                                       .standing, .root_env [proximal ramp,
%                                       see chlamy_moment_fun]); Q_Nm and omega_Hz may be
%                                       vectors (one entry per flagellum):
%                                       use the custom moment of
%                                       chlamy_moment_fun instead of
%                                       traveling_wave (peak moment Q_Nm [N m],
%                                       wavelength lambda_rel*L, taper from
%                                       s_taper); phase/frequency are taken
%                                       from the JSON
%                    .moment_table      table from chlamy_inverse_moment: tabulated
%                                       q(s,t) of the measured gait for every
%                                       flagellum (chlamy_table_moment_fun)
%                    .tilt_deg          [g1 g2] beat-plane tilt [deg]: rotation of each
%                                       flagellum about its own root tangent
%                                       (right-handed about the outward tangent);
%                                       the drive moment acts about the tilted
%                                       plane normal (used with .moment_table)
%                    .drive_amp, .drive_phase, .drive_freq_ratio
%                                       per-flagellum amplitude factor, beat
%                                       phase [rad] and frequency ratio of the
%                                       tabulated drive (.moment_table; see
%                                       chlamy_table_moment_fun)
%                    .pili_index        restrict amplitude/wavenumber/omega
%                                       overrides to these flagella
%
%  Outputs:
%    results, params, output_dir - as returned by run_simulation
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-06

proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));
cd(proj_root);   % run_simulation writes to <pwd>/output/<case_name>

if nargin < 2, ov = struct(); end
cfg = load_sim_config(fullfile(proj_root, config_file));

if isfield(ov, 'case_name'),  cfg.case_name  = ov.case_name;  end
if isfield(ov, 'T_obs'),      cfg.T_obs      = ov.T_obs;      end
if isfield(ov, 'dt'),         cfg.dt         = ov.dt;         end
if isfield(ov, 'seg_per_um'), cfg.seg_per_um = ov.seg_per_um; end
if isfield(ov, 'hydro_radius_factor'), cfg.hydro_radius_factor = ov.hydro_radius_factor; end
if isfield(ov, 'root_half_angle_deg')
    b_rh = deg2rad(ov.root_half_angle_deg);
    cfg.custom_geometry.attachment_sites(1).direction = [cos(b_rh);  sin(b_rh); 0];
    cfg.custom_geometry.attachment_sites(2).direction = [cos(b_rh); -sin(b_rh); 0];
end
tilt_rad = zeros(1, 2);
if isfield(ov, 'tilt_deg')
    % beat-plane tilt: rigid rotation of each flagellum about its own root
    % tangent (right-handed about the outward tangent), custom_pili_building
    % 'deflection_z'; the root tangent and the attachment are unchanged
    tilt_rad(1:numel(ov.tilt_deg)) = deg2rad(ov.tilt_deg(:).');
    for k = 1:numel(cfg.custom_geometry.pili)
        if iscell(cfg.custom_geometry.pili)
            cfg.custom_geometry.pili{k}.deflection_z = tilt_rad(k);
        else
            cfg.custom_geometry.pili(k).deflection_z = tilt_rad(k);
        end
    end
end
if isfield(ov, 'rel_tol'),    cfg.RelTol     = ov.rel_tol;    end
if isfield(ov, 'abs_tol'),    cfg.AbsTol     = ov.abs_tol;    end
cfg.NT = round(cfg.T_obs / cfg.dt) + 1;

idx = 1:numel(cfg.pili_ext);
if isfield(ov, 'pili_index'), idx = ov.pili_index; end
for k = idx
    am = cfg.pili_ext(k).active_moment;
    if isstruct(am) && strcmp(am.type, 'traveling_wave')
        if isfield(ov, 'amplitude_Nm'), am.amplitude_Nm = ov.amplitude_Nm; end
        if isfield(ov, 'wavenumber'),   am.wavenumber   = ov.wavenumber;   end
        if isfield(ov, 'omega_Hz') && ~isfield(ov, 'moment')
            am.omega_Hz = ov.omega_Hz;
        end
        cfg.pili_ext(k).active_moment = am;
    end
end

% Optional: replace the built-in traveling_wave by the custom moment of
% chlamy_moment_fun (uniform wave of wavelength lambda_rel*L with a tip taper)
if isfield(ov, 'moment')
    mo = ov.moment;
    p0 = chlamy_build_params(cfg, @(varargin) []);   % scales only
    EI = p0.scales.EI;  T_0 = p0.scales.T_0;  L_0 = p0.scales.L_0;
    for k = idx
        am = cfg.pili_ext(k).active_moment;
        if ~(isstruct(am) && strcmp(am.type, 'traveling_wave'))
            error('chlamy_run_case:Moment', ...
                'Flagellum %d must have a traveling_wave placeholder in the JSON.', k);
        end
        f_Hz = am.omega_Hz;
        if isfield(ov, 'omega_Hz')      % scalar, or one value per flagellum
            f_Hz = ov.omega_Hz(min(k, numel(ov.omega_Hz)));
        end
        phi = 0;
        if isfield(am, 'phase_rad'), phi = am.phase_rad; end
        Lk = p0.ell(k) * L_0;
        newam = struct('type', 'custom', 'phase_rad', phi, 'omega_Hz', f_Hz, ...
            'Q_Nm', mo.Q_Nm(min(k, numel(mo.Q_Nm))), ...
            'lambda_rel', getf(mo, 'lambda_rel', 0.70), ...
            's_taper', getf(mo, 's_taper', 0.6), ...
            'standing', getf(mo, 'standing', false));
        newam.fun = chlamy_moment_fun(newam.Q_Nm, f_Hz, newam.lambda_rel, phi, ...
            newam.s_taper, EI, Lk, T_0, newam.standing, getf(mo, 'root_env', []));
        cfg.pili_ext(k).active_moment = newam;
    end
end

% Optional: tabulated active moment q(s,t) from the measured gait
% (chlamy_inverse_moment); one table per flagellum
if isfield(ov, 'moment_table')
    tab = ov.moment_table;
    p0 = chlamy_build_params(cfg, @(varargin) []);   % scales only
    assert(abs(p0.scales.T_0 - tab.T_0) < 1e-12 * tab.T_0, 'chlamy_run_case:MomentTable', ...
        'Time scale of the moment table differs from this configuration.');
    for k = idx
        am = cfg.pili_ext(k).active_moment;
        if ~(isstruct(am) && strcmp(am.type, 'traveling_wave'))
            error('chlamy_run_case:Moment', ...
                'Flagellum %d must have a traveling_wave placeholder in the JSON.', k);
        end
        newam = struct('type', 'custom', 'phase_rad', 0, 'omega_Hz', tab.f_Hz, ...
            'Q_Nm', tab.Q_peak_Nm(k), 'lambda_rel', 0, 's_taper', 0, 'standing', false);
        % beat-plane normal: planar-gait normal z, rotated about the root tangent by the tilt,
        % expressed in the (body-fixed) clamped root frame
        t0 = p0.d_local{k}(:);
        g  = tilt_rad(k);
        n0 = [0; 0; 1] * cos(g) + cross(t0, [0; 0; 1]) * sin(g) + t0 * dot(t0, [0; 0; 1]) * (1 - cos(g));
        dopt = struct('normal_local', p0.D_local{k}.' * n0, ...
            'amp',        pick(getf(ov, 'drive_amp', 1), k), ...
            'phase',      pick(getf(ov, 'drive_phase', 0), k), ...
            'freq_ratio', pick(getf(ov, 'drive_freq_ratio', 1), k));
        newam.fun = chlamy_table_moment_fun(tab, k, dopt);
        cfg.pili_ext(k).active_moment = newam;
    end
end

cfg.auto_confirm       = true;
cfg.preview_initial    = isfield(ov, 'preview') && ov.preview;   % the vector PDF is ~40 MB; geometry is audited by chlamy_static_check
cfg.compute_flowfield  = false;
cfg.generate_animation = isfield(ov, 'animation') && ov.animation;

[results, params, output_dir] = run_simulation(cfg);
end

function v = pick(x, k)
v = x(min(k, numel(x)));
end

function v = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
