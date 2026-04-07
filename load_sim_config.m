function cfg = load_sim_config(json_path)
%LOAD_SIM_CONFIG  Load simulation configuration from a JSON file
%
%  Description:
%    Reads a JSON configuration file and returns a flat struct with all
%    user-configurable parameters.  Missing fields are filled with default
%    values identical to the inline defaults in main_coupled_simulation.m.
%    Pili defined as a struct array in the JSON are automatically converted
%    to a [8 x M] Pili_Matrix via parse_pili_config.
%
%  Syntax:
%    cfg = load_sim_config(json_path)
%
%  Input:
%    json_path - Path to a JSON configuration file
%
%  Output:
%    cfg - Flat struct with fields matching workspace variables in main:
%      case_name, minor_axis, major_axis, N_body, alpha_0,
%      Pili_Matrix, seg_per_um, E, r_pili, sigma_p,
%      mu, rho_bac, rho_fluid, T_obs, AbsTol, RelTol,
%      delta, tlim, dt, enable_brown, Temper, dt_B, fixman_midpoint,
%      preview_initial, generate_animation, animation_fps, animation_skip,
%      auto_confirm, flowfield_grid_spacing, flowfield_grid_range,
%      flowfield_quiver_skip, flowfield_chunk_size
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17
%  Updated: 2026-03-17 — v3.2: Fixed dead code in get_json (islogical guard)
%  Updated: 2026-04-03 — v5.0: background_flow, per-pili intrinsic_curvature / active_moment

%% Read and decode JSON
fid = fopen(json_path, 'r', 'n', 'UTF-8');
if fid == -1
    error('load_sim_config:FileNotFound', 'Cannot open config file: %s', json_path);
end
raw = fread(fid, '*char')';
fclose(fid);
j = jsondecode(raw);

%% Map nested JSON fields to flat cfg with defaults
cfg.case_name  = get_json(j, 'case_name', '');

% Body
jb = get_json(j, 'body', struct());
cfg.minor_axis = get_json(jb, 'minor_axis', 0.4e-6);
cfg.major_axis = get_json(jb, 'major_axis', 1.25e-6);
cfg.N_body     = get_json(jb, 'N_body',     1000);
cfg.alpha_0    = get_json(jb, 'alpha_0',    0);
cfg.eps_body   = get_json(jb, 'eps_body',   []);

% Pili — accept struct array (named) or raw matrix
jp = get_json(j, 'pili', []);
if isstruct(jp) && ~isempty(jp)
    [cfg.Pili_Matrix, cfg.pili_ext] = parse_pili_config(jp);
elseif isnumeric(jp) && ~isempty(jp)
    cfg.Pili_Matrix = jp;
    cfg.pili_ext    = struct('intrinsic_curvature', {}, 'active_moment', {});
else
    cfg.Pili_Matrix = zeros(8, 0);
    cfg.pili_ext    = struct('intrinsic_curvature', {}, 'active_moment', {});
end

% Material
jm = get_json(j, 'material', struct());
cfg.E       = get_json(jm, 'E',       1.768e7);
cfg.r_pili  = get_json(jm, 'r_pili',  3e-9);
cfg.sigma_p = get_json(jm, 'sigma_p', 0.2);

% Fluid
jf = get_json(j, 'fluid', struct());
cfg.mu        = get_json(jf, 'mu',        1e-3);
cfg.rho_bac   = get_json(jf, 'rho_bac',   1180.0);
cfg.rho_fluid = get_json(jf, 'rho_fluid', 998.232);

% Simulation
js = get_json(j, 'simulation', struct());
cfg.T_obs      = get_json(js, 'T_obs',      30.0);
cfg.AbsTol     = get_json(js, 'AbsTol',     1e-4);
cfg.RelTol     = get_json(js, 'RelTol',     1e-4);
delta_frac     = get_json(js, 'delta_frac', 20);
cfg.delta      = pi / delta_frac;
cfg.tlim       = get_json(js, 'tlim',       120);
cfg.dt         = get_json(js, 'dt',         0.1);
cfg.seg_per_um = get_json(js, 'seg_per_um', 15);

% Brownian
jbr = get_json(j, 'brownian', struct());
cfg.enable_brown     = logical(get_json(jbr, 'enable', false));
cfg.Temper           = get_json(jbr, 'temperature', 293.15);
dt_B_ratio           = get_json(jbr, 'dt_B_ratio', 0.5);
cfg.dt_B             = cfg.dt * dt_B_ratio;
cfg.fixman_midpoint  = logical(get_json(jbr, 'fixman_midpoint', false));

% Output / visualisation
jo = get_json(j, 'output', struct());
cfg.preview_initial    = logical(get_json(jo, 'preview_initial',    true));
cfg.generate_animation = logical(get_json(jo, 'generate_animation', true));
cfg.animation_fps      = get_json(jo, 'animation_fps',  30);
cfg.animation_skip     = get_json(jo, 'animation_skip', 1);
cfg.auto_confirm       = logical(get_json(jo, 'auto_confirm', false));

% Flow field
jff = get_json(j, 'flowfield', struct());
cfg.flowfield_grid_spacing = get_json(jff, 'grid_spacing_um', 0.05);
cfg.flowfield_grid_range   = get_json(jff, 'grid_range_um',   10.0);
cfg.flowfield_quiver_skip  = get_json(jff, 'quiver_skip',     5);
cfg.flowfield_chunk_size   = get_json(jff, 'chunk_size',      2000);

% Background flow
jbg = get_json(j, 'background_flow', struct());
cfg.background_flow.type           = get_json(jbg, 'type', 'off');
cfg.background_flow.U_max          = get_json(jbg, 'U_max', 0);
cfg.background_flow.direction      = get_json(jbg, 'direction', [1, 0, 0]);
cfg.background_flow.shear_axis     = get_json(jbg, 'shear_axis', [0, 0, 1]);
cfg.background_flow.shear_height_m = get_json(jbg, 'shear_height_m', 0);

%% Derive N_segments
if size(cfg.Pili_Matrix, 2) > 0
    pili_lengths_um  = cfg.Pili_Matrix(2, :);
    cfg.N_segments   = max(round(pili_lengths_um * cfg.seg_per_um), 10);
else
    cfg.N_segments   = zeros(1, 0);
end
cfg.NT = round(cfg.T_obs / cfg.dt) + 1;

end


function v = get_json(s, fname, default)
%GET_JSON  Read field from decoded JSON struct with fallback default
    if isstruct(s) && isfield(s, fname)
        v = s.(fname);
        if islogical(default) && isnumeric(v)
            v = logical(v);
        end
    else
        v = default;
    end
end
