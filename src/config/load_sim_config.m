function cfg = load_sim_config(json_path)
%LOAD_SIM_CONFIG  Load simulation configuration from a JSON file
%
%  Description:
%    Reads a UTF-8 JSON configuration file and maps its nested schema to the
%    flat configuration struct consumed by run_simulation.  Missing JSON
%    fields use the standard sinking-workflow fallbacks defined in this
%    function.  configs/config_example.json is an illustrative configuration
%    and intentionally overrides some fallbacks (for example case_name and
%    T_obs); it is not a canonical default-value table.  The fallbacks are
%    also distinct from the intentional motor-diagnostic inline demo in
%    main_coupled_simulation.m.  Named parametric pili are converted to a
%    fixed [8 x M] Pili_Matrix and a parallel pili_ext struct array by
%    parse_pili_config.
%
%  Syntax:
%    cfg = load_sim_config(json_path)
%
%  Input:
%    json_path - Path to a UTF-8 JSON configuration file
%
%  Output:
%    cfg - Flat struct containing:
%      Identity and geometry:
%        case_name, geometry_mode, minor_axis [m], major_axis [m], N_body,
%        alpha_0 [rad], eps_body [L_0], initial_orientation, Pili_Matrix,
%        pili_ext, plus custom_geometry in custom mode.  Pili_Matrix row 2
%        is length [um], rows 3-4 are angles [rad], and other rows are
%        dimensionless indices or shape parameters.  The matrix must have
%        exactly eight rows; per-pilus stiffness rows are rejected.
%      Material and fluid:
%        E [Pa], r_pili [m], sigma_p, mu [Pa*s], rho_bac [kg/m^3],
%        rho_fluid [kg/m^3]
%      Time integration and Brownian motion:
%        T_obs [s], AbsTol, RelTol, delta [rad], tlim [s], dt [s],
%        seg_per_um [1/um], enable_brown, Temper [K], dt_B [s],
%        brown_rng_seed
%      Output and flow-field controls:
%        preview_initial, generate_animation, animation_fps,
%        animation_skip, auto_confirm, compute_flowfield,
%        flowfield_grid_spacing [um], flowfield_grid_range [um],
%        flowfield_quiver_skip, flowfield_chunk_size
%      Background flow:
%        background_flow.type, .U_max [m/s], .direction, .shear_axis,
%        .shear_height_m [m]
%      Derived values:
%        N_segments and NT
%
%    JSON-to-cfg mapping:
%      case_name, geometry_mode -> cfg fields with the same names
%      custom_geometry -> custom_geometry and custom-mode geometry fields
%      custom_geometry.initial_orientation -> initial_orientation
%      body.* -> minor_axis, major_axis, N_body, alpha_0, eps_body
%      pili -> Pili_Matrix and pili_ext
%      material.* -> E, r_pili, sigma_p
%      fluid.* -> mu, rho_bac, rho_fluid
%      simulation.{T_obs,AbsTol,RelTol,tlim,dt,seg_per_um} -> same names
%      simulation.delta_frac -> delta = pi/delta_frac
%      brownian.enable -> enable_brown
%      brownian.temperature -> Temper
%      brownian.dt_B_ratio -> dt_B = dt*dt_B_ratio
%      brownian.rng_seed -> brown_rng_seed
%      output.* -> preview_initial, generate_animation, animation_fps,
%        animation_skip, auto_confirm, compute_flowfield
%      flowfield.grid_spacing_um -> flowfield_grid_spacing
%      flowfield.grid_range_um -> flowfield_grid_range
%      flowfield.quiver_skip -> flowfield_quiver_skip
%      flowfield.chunk_size -> flowfield_chunk_size
%      background_flow.* -> fields of cfg.background_flow
%
%  Errors:
%    Throws load_sim_config:FileNotFound when json_path cannot be opened.
%    JSON decoding, struct normalisation, and pili parsing errors propagate
%    from jsondecode, normalise_json_struct_array, and parse_pili_config.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17
%  Updated: 2026-07-13

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

% Geometry mode: 'parametric' (default) or 'custom'
cfg.geometry_mode = get_json(j, 'geometry_mode', 'parametric');

switch cfg.geometry_mode
    case 'custom'
        % --- Custom geometry mode ---
        jcg = get_json(j, 'custom_geometry', struct());
        cfg.custom_geometry = jcg;

        % Set parametric placeholders to NaN/empty (enforce no-mixing)
        cfg.minor_axis   = NaN;
        cfg.major_axis   = NaN;
        cfg.N_body       = NaN;
        cfg.Pili_Matrix  = zeros(8, 0);

        % Parse per-pilus extended fields from custom_geometry.pili (M1)
        cfg.pili_ext = parse_custom_pili_ext(get_json(jcg, 'pili', []));

        % Default alpha_0 and eps_body for custom mode
        cfg.alpha_0  = get_json(jcg, 'alpha_0', 0);
        cfg.eps_body = get_json(jcg, 'eps_body', []);

        % Initial orientation (C2) — overrides alpha_0 when set
        cfg.initial_orientation = get_json(jcg, 'initial_orientation', []);

    otherwise
        % --- Parametric mode (default) ---
        jb = get_json(j, 'body', struct());
        cfg.minor_axis = get_json(jb, 'minor_axis', 0.4e-6);
        cfg.major_axis = get_json(jb, 'major_axis', 1.25e-6);
        cfg.N_body     = get_json(jb, 'N_body',     300);
        cfg.alpha_0    = get_json(jb, 'alpha_0',    0);
        cfg.eps_body   = get_json(jb, 'eps_body',   []);

        % initial_orientation is a custom-mode-only feature; leave empty
        cfg.initial_orientation = [];

        % Pili — accept struct array (named), cell-array-of-struct (heterogeneous
        % JSON), or raw matrix.  jsondecode emits a cell array when the JSON
        % objects in `pili` have different field sets (e.g. one entry has
        % `motor_torque` and another does not, as in
        % config_motor_flagella_example.json).  normalise_json_struct_array
        % collapses both shapes into a single struct array using the union of
        % fields.  Raw [8 x M] numeric matrices are still accepted for
        % backward compatibility.
        jp = get_json(j, 'pili', []);
        if isnumeric(jp) && ~isempty(jp)
            if size(jp, 1) ~= 8
                error('load_sim_config:InvalidPiliMatrixRows', ...
                    ['Raw numeric pili matrices must have exactly 8 rows; ', ...
                     'got %d. Pilus stiffness is derived from the global ', ...
                     'material E and r_pili values.'], size(jp, 1));
            end
            cfg.Pili_Matrix = jp;
            cfg.pili_ext    = struct('intrinsic_curvature', {}, 'active_moment', {}, ...
                                     'root_coupling', {}, 'motor_torque', {});
        elseif isempty(jp)
            cfg.Pili_Matrix = zeros(8, 0);
            cfg.pili_ext    = struct('intrinsic_curvature', {}, 'active_moment', {}, ...
                                     'root_coupling', {}, 'motor_torque', {});
        else
            jp_norm = normalise_json_struct_array(jp);
            if isempty(jp_norm)
                cfg.Pili_Matrix = zeros(8, 0);
                cfg.pili_ext    = struct('intrinsic_curvature', {}, 'active_moment', {}, ...
                                         'root_coupling', {}, 'motor_torque', {});
            else
                [cfg.Pili_Matrix, cfg.pili_ext] = parse_pili_config(jp_norm);
            end
        end
end

% Material
jm = get_json(j, 'material', struct());
cfg.E       = get_json(jm, 'E',       3.279e8);
cfg.r_pili  = get_json(jm, 'r_pili',  3e-9);
cfg.sigma_p = get_json(jm, 'sigma_p', 0.2);
% Pilus hydrodynamic regularisation length / physical radius (default sqrt(2))
cfg.hydro_radius_factor = get_json(jm, 'hydro_radius_factor', sqrt(2));

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
cfg.tlim       = get_json(js, 'tlim',       3600);
cfg.dt         = get_json(js, 'dt',         0.1);
cfg.seg_per_um = get_json(js, 'seg_per_um', 10);

% Brownian
jbr = get_json(j, 'brownian', struct());
cfg.enable_brown     = logical(get_json(jbr, 'enable', false));
cfg.Temper           = get_json(jbr, 'temperature', 293.15);
dt_B_ratio           = get_json(jbr, 'dt_B_ratio', 0.5);
cfg.dt_B             = cfg.dt * dt_B_ratio;

% Optional: positive scalar = fixed twister seed; omit / NaN / nonpositive = time-based default
if isfield(jbr, 'rng_seed')
    v = jbr.rng_seed;
    if isempty(v) || ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v <= 0
        cfg.brown_rng_seed = NaN;
    else
        cfg.brown_rng_seed = double(v);
    end
else
    cfg.brown_rng_seed = NaN;
end

% Output / visualisation
jo = get_json(j, 'output', struct());
cfg.preview_initial    = logical(get_json(jo, 'preview_initial',    true));
cfg.generate_animation = logical(get_json(jo, 'generate_animation', true));
cfg.animation_fps      = get_json(jo, 'animation_fps',  30);
cfg.animation_skip     = get_json(jo, 'animation_skip', 1);
cfg.auto_confirm       = logical(get_json(jo, 'auto_confirm', false));
cfg.compute_flowfield  = logical(get_json(jo, 'compute_flowfield', true));

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
if strcmp(cfg.geometry_mode, 'custom')
    % Custom mode: N_segments derived from actual pili lengths (computed later)
    cfg.N_segments = zeros(1, 0);
else
    if size(cfg.Pili_Matrix, 2) > 0
        pili_lengths_um  = cfg.Pili_Matrix(2, :);
        cfg.N_segments   = max(round(pili_lengths_um * cfg.seg_per_um), 10);
    else
        cfg.N_segments   = zeros(1, 0);
    end
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


function ext = parse_custom_pili_ext(pili_json)
%PARSE_CUSTOM_PILI_EXT  Build pili_ext struct array from custom-mode JSON.
%
%  Extracts per-pilus optional fields (intrinsic_curvature, active_moment,
%  root_coupling, motor_torque) from custom_geometry.pili and returns a
%  struct array with the same schema as parse_pili_config (parametric mode).
%  Missing fields are filled with the same defaults as parse_pili_config:
%  intrinsic_curvature = 'frozen', active_moment = struct('type','off'),
%  root_coupling = 'clamped', motor_torque = struct('type','off').
%
%  Inputs:
%    pili_json - JSON-decoded pili array (struct array / cell array / [])
%
%  Output:
%    ext - [1 x M] struct array with fields .intrinsic_curvature,
%          .active_moment, .root_coupling, .motor_torque
%          (M = 0 if pili_json is empty)

    if isempty(pili_json)
        ext = struct('intrinsic_curvature', {}, 'active_moment', {}, ...
                     'root_coupling', {}, 'motor_torque', {});
        return;
    end

    % Handle JSON heterogeneous arrays (cell-of-struct)
    pili = normalise_json_struct_array(pili_json);
    M    = numel(pili);

    if M == 0
        ext = struct('intrinsic_curvature', {}, 'active_moment', {}, ...
                     'root_coupling', {}, 'motor_torque', {});
        return;
    end

    ext(M) = struct('intrinsic_curvature', 'frozen', ...
                    'active_moment', struct('type', 'off'), ...
                    'root_coupling', 'clamped', ...
                    'motor_torque', struct('type', 'off'));
    for k = 1:M
        p = pili(k);

        if isfield(p, 'intrinsic_curvature') && ~isempty(p.intrinsic_curvature)
            ext(k).intrinsic_curvature = p.intrinsic_curvature;
        else
            ext(k).intrinsic_curvature = 'frozen';
        end

        if isfield(p, 'active_moment') && ~isempty(p.active_moment)
            ext(k).active_moment = p.active_moment;
        else
            ext(k).active_moment = struct('type', 'off');
        end

        if isfield(p, 'root_coupling') && ~isempty(p.root_coupling)
            ext(k).root_coupling = char(p.root_coupling);
        else
            ext(k).root_coupling = 'clamped';
        end

        if isfield(p, 'motor_torque') && ~isempty(p.motor_torque)
            ext(k).motor_torque = p.motor_torque;
        else
            ext(k).motor_torque = struct('type', 'off');
        end
    end
end
