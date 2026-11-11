function [params, Z0, cfg] = chlamy_build_params(cfg_or_path, logf_h)
%CHLAMY_BUILD_PARAMS  Build params and initial state exactly as run_simulation does
%
%  Description:
%    Loads a JSON configuration (or takes a cfg struct), validates it and
%    runs init_coupled_params / init_coupled_state with the same argument
%    list as run_simulation.m, without integrating the ODE.  Used by the
%    static checks of the Chlamydomonas validation plan.
%
%  Syntax:
%    [params, Z0, cfg] = chlamy_build_params('configs/config_chlamy_single_flagellum.json')
%    [params, Z0, cfg] = chlamy_build_params(cfg, logf_h)
%
%  Inputs:
%    cfg_or_path - JSON path (relative to the project root) or a cfg struct
%    logf_h      - (Optional) logger handle, default fprintf
%
%  Outputs:
%    params - Parameter structure from init_coupled_params
%    Z0     - Initial state vector from init_coupled_state
%    cfg    - The cfg struct actually used
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-06

proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));

if nargin < 2 || isempty(logf_h)
    logf_h = @(varargin) fprintf(varargin{:});
end

if ischar(cfg_or_path) || isstring(cfg_or_path)
    cfg = load_sim_config(fullfile(proj_root, char(cfg_or_path)));
else
    cfg = cfg_or_path;
end

validate_simulation_params(cfg, logf_h);

hrf = sqrt(2);
if isfield(cfg, 'hydro_radius_factor') && ~isempty(cfg.hydro_radius_factor)
    hrf = cfg.hydro_radius_factor;
end
eps_body_cfg = [];
if isfield(cfg, 'eps_body') && ~isempty(cfg.eps_body)
    eps_body_cfg = cfg.eps_body;
end
init_orient = [];
if isfield(cfg, 'initial_orientation') && ~isempty(cfg.initial_orientation)
    init_orient = cfg.initial_orientation;
end

ic  = {cfg.pili_ext.intrinsic_curvature};
am  = {cfg.pili_ext.active_moment};
rc  = {cfg.pili_ext.root_coupling};
mt  = {cfg.pili_ext.motor_torque};

params = init_coupled_params( ...
    'geometry_mode',       cfg.geometry_mode, ...
    'custom_geometry',     cfg.custom_geometry, ...
    'N_body',              NaN, ...
    'minor_axis',          NaN, ...
    'major_axis',          NaN, ...
    'alpha_0',             cfg.alpha_0, ...
    'initial_orientation', init_orient, ...
    'mu',                  cfg.mu, ...
    'E',                   cfg.E, ...
    'r_pili',              cfg.r_pili, ...
    'hydro_radius_factor', hrf, ...
    'sigma',               cfg.sigma_p, ...
    'rho_bac',             cfg.rho_bac, ...
    'rho_fluid',           cfg.rho_fluid, ...
    'N_segments',          zeros(1, 0), ...
    'Pili_Matrix',         cfg.Pili_Matrix, ...
    'seg_per_um',          cfg.seg_per_um, ...
    'delta',               cfg.delta, ...
    'Temper',              cfg.Temper, ...
    'enable_brown',        cfg.enable_brown, ...
    'intrinsic_curvature', ic, ...
    'active_moment',       am, ...
    'root_coupling',       rc, ...
    'motor_torque',        mt, ...
    'background_flow',     cfg.background_flow, ...
    'eps_body',            eps_body_cfg, ...
    'logf',                logf_h);

Z0 = init_coupled_state(params, logf_h);
end
