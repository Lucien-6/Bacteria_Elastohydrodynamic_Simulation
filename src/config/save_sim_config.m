function save_sim_config(cfg, json_path)
%SAVE_SIM_CONFIG  Save current simulation configuration to a JSON file
%
%  Description:
%    Serialises the supported user-facing configuration fields to a
%    human-readable UTF-8 JSON file.  The resulting file is compatible with
%    load_sim_config, but it is not an exact snapshot of every cfg field or
%    all runtime state.
%
%    Round-trip limitations:
%      - Derived fields such as N_segments and NT are recomputed on load.
%      - delta is saved as the integer round(pi/delta) and reloads as
%        pi/delta_frac, so arbitrary delta values are quantised.
%      - Runtime function handles in active_moment and motor_torque entries
%        are removed because JSON cannot represent function handles.
%      - A Brownian RNG seed is omitted when Brownian motion is disabled or
%        the seed is absent/nonpositive, and therefore reloads as NaN.
%      - Streamline spacing and length are omitted when unset, and reload
%        as NaN.  The panel then uses a fraction of its shorter side.
%      - Floating-point values follow MATLAB jsonencode/jsondecode semantics.
%      - Runtime-only fields, including silent, progress_queue, progress_tag,
%        logger handles, and any other fields outside the schema written
%        below, are not preserved.
%
%  Syntax:
%    save_sim_config(cfg, json_path)
%
%  Inputs:
%    cfg       - Flat configuration struct (as used in main_coupled_simulation)
%    json_path - Destination UTF-8 JSON path; an existing file is replaced
%
%  Side Effects and Warnings:
%    Creates or truncates json_path.  If the destination cannot be opened,
%    emits save_sim_config:WriteError and returns without throwing.  JSON
%    encoding errors propagate; fwrite and fclose completion is not checked.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17
%  Updated: 2026-07-12
%  Updated: 2026-11-12 — v9.0.0: brownian.pili_thermal / pili_mode_orders
%  Updated: 2026-11-13 — optional flowfield stream_spacing_um / stream_length_um

%% Build nested structure matching config_example.json layout
j.case_name    = cfg.case_name;
j.geometry_mode = cfg.geometry_mode;

switch cfg.geometry_mode
    case 'custom'
        % --- Custom geometry mode ---
        j.custom_geometry = cfg.custom_geometry;
        if isfield(cfg, 'alpha_0')
            j.custom_geometry.alpha_0 = cfg.alpha_0;
        end
        if isfield(cfg, 'eps_body') && ~isempty(cfg.eps_body)
            j.custom_geometry.eps_body = cfg.eps_body;
        end
        if isfield(cfg, 'initial_orientation') && ~isempty(cfg.initial_orientation)
            j.custom_geometry.initial_orientation = cfg.initial_orientation;
        end

        % Re-inject per-pilus extended fields into custom_geometry.pili,
        % so the saved JSON round-trips back through load_sim_config.
        if isfield(cfg, 'pili_ext') && ~isempty(cfg.pili_ext) ...
                && isfield(j.custom_geometry, 'pili') ...
                && ~isempty(j.custom_geometry.pili)
            pili_raw = normalise_json_struct_array(j.custom_geometry.pili);
            M_ext = numel(cfg.pili_ext);
            for kk = 1:min(numel(pili_raw), M_ext)
                pe = cfg.pili_ext(kk);
                pili_raw(kk).intrinsic_curvature = pe.intrinsic_curvature;
                am_val = pe.active_moment;
                if isstruct(am_val) && isfield(am_val, 'fun')
                    am_val = rmfield(am_val, 'fun');
                end
                pili_raw(kk).active_moment = am_val;
                if isfield(pe, 'root_coupling')
                    pili_raw(kk).root_coupling = pe.root_coupling;
                end
                if isfield(pe, 'motor_torque')
                    mt_val = pe.motor_torque;
                    if isstruct(mt_val) && isfield(mt_val, 'fun')
                        mt_val = rmfield(mt_val, 'fun');
                    end
                    pili_raw(kk).motor_torque = mt_val;
                end
            end
            j.custom_geometry.pili = pili_raw;
        end
        j.pili = [];

    otherwise
        % --- Parametric mode ---
        j.body.minor_axis = cfg.minor_axis;
        j.body.major_axis = cfg.major_axis;
        j.body.N_body     = cfg.N_body;
        j.body.alpha_0    = cfg.alpha_0;
        if isfield(cfg, 'eps_body') && ~isempty(cfg.eps_body)
            j.body.eps_body = cfg.eps_body;
        end

        % Pili: convert [8 x M] matrix back to struct array (or cell array if heterogeneous)
        M = size(cfg.Pili_Matrix, 2);
        type_names = {'linear', 'parabolic', 'circular', 'sinusoidal', 'spiral'};
        pili_cells = cell(1, M);
        for k = 1:M
            p_entry = struct();
            type_idx = cfg.Pili_Matrix(1, k);
            if ~(isnumeric(type_idx) && isscalar(type_idx) && isfinite(type_idx) && ...
                    type_idx >= 1 && type_idx <= numel(type_names) && floor(type_idx) == type_idx)
                error('save_sim_config:InvalidPiliType', ...
                    'Pili %d: Pili_Matrix(1, %d) must be integer 1-%d, got %g.', ...
                    k, k, numel(type_names), type_idx);
            end
            p_entry.type         = type_names{type_idx};
            p_entry.length_um    = cfg.Pili_Matrix(2, k);
            p_entry.deflection_z = cfg.Pili_Matrix(3, k);
            p_entry.deflection_x = cfg.Pili_Matrix(4, k);
            sp = cfg.Pili_Matrix(5:7, k)';
            sp(isnan(sp)) = [];
            if isempty(sp)
                p_entry.shape_params = [];
            else
                p_entry.shape_params = sp;
            end
            p_entry.site        = cfg.Pili_Matrix(8, k);

            % Extended per-pili fields
            if isfield(cfg, 'pili_ext') && k <= numel(cfg.pili_ext)
                p_entry.intrinsic_curvature = cfg.pili_ext(k).intrinsic_curvature;
                am_val = cfg.pili_ext(k).active_moment;
                if isstruct(am_val) && isfield(am_val, 'fun')
                    am_val = rmfield(am_val, 'fun');
                end
                p_entry.active_moment = am_val;

                % Root coupling type
                if isfield(cfg.pili_ext(k), 'root_coupling')
                    p_entry.root_coupling = cfg.pili_ext(k).root_coupling;
                end
                % Motor torque configuration
                if isfield(cfg.pili_ext(k), 'motor_torque')
                    mt_val = cfg.pili_ext(k).motor_torque;
                    if isstruct(mt_val) && isfield(mt_val, 'fun')
                        mt_val = rmfield(mt_val, 'fun');
                    end
                    p_entry.motor_torque = mt_val;
                end
            end

            pili_cells{k} = p_entry;
        end
        if M == 0
            j.pili = [];
        else
            % Pack as struct array when all entries share identical field sets,
            % otherwise keep cell array so jsonencode serialises heterogeneous pili.
            f0 = fieldnames(pili_cells{1});
            same_fields = true;
            for k = 2:M
                if ~isequal(fieldnames(pili_cells{k}), f0)
                    same_fields = false;
                    break;
                end
            end
            if same_fields
                j.pili = [pili_cells{:}];
            else
                j.pili = pili_cells;
            end
        end
end

j.material.E       = cfg.E;
j.material.r_pili  = cfg.r_pili;
j.material.sigma_p = cfg.sigma_p;
if isfield(cfg, 'hydro_radius_factor') && ~isempty(cfg.hydro_radius_factor)
    j.material.hydro_radius_factor = cfg.hydro_radius_factor;
else
    j.material.hydro_radius_factor = sqrt(2);   % default used by init_coupled_params
end

j.fluid.mu        = cfg.mu;
j.fluid.rho_bac   = cfg.rho_bac;
j.fluid.rho_fluid = cfg.rho_fluid;

j.simulation.T_obs      = cfg.T_obs;
j.simulation.AbsTol     = cfg.AbsTol;
j.simulation.RelTol     = cfg.RelTol;
j.simulation.delta_frac = round(pi / cfg.delta);
j.simulation.tlim       = cfg.tlim;
j.simulation.dt         = cfg.dt;
j.simulation.seg_per_um = cfg.seg_per_um;

j.brownian.enable      = cfg.enable_brown;
j.brownian.temperature = cfg.Temper;
if cfg.dt > 0
    j.brownian.dt_B_ratio = cfg.dt_B / cfg.dt;
else
    j.brownian.dt_B_ratio = 0.5;
end
if logical(cfg.enable_brown) && isfield(cfg, 'brown_rng_seed') ...
        && isscalar(cfg.brown_rng_seed) && isnumeric(cfg.brown_rng_seed) ...
        && isfinite(cfg.brown_rng_seed) && cfg.brown_rng_seed > 0
    j.brownian.rng_seed = cfg.brown_rng_seed;
end
% Pilus thermal fluctuations (v9.0.0): written only when switched on, so the
% snapshot of every other run keeps its previous layout.
if isfield(cfg, 'enable_pili_brown') && ~isempty(cfg.enable_pili_brown) ...
        && logical(cfg.enable_pili_brown)
    j.brownian.pili_thermal = true;
    if isfield(cfg, 'pili_brown_orders') && ~isempty(cfg.pili_brown_orders)
        j.brownian.pili_mode_orders = cfg.pili_brown_orders;
    else
        j.brownian.pili_mode_orders = 2;
    end
    if isfield(cfg, 'pili_brown_tau_cut_s') && ~isempty(cfg.pili_brown_tau_cut_s) ...
            && isnumeric(cfg.pili_brown_tau_cut_s) && isscalar(cfg.pili_brown_tau_cut_s) ...
            && isfinite(cfg.pili_brown_tau_cut_s) && cfg.pili_brown_tau_cut_s > 0
        j.brownian.pili_tau_cut_s = double(cfg.pili_brown_tau_cut_s);
    end
    if isfield(cfg, 'pili_brown_protect_pairs') && ~isempty(cfg.pili_brown_protect_pairs)
        j.brownian.pili_protect_pairs = logical(cfg.pili_brown_protect_pairs);
    end
    if isfield(cfg, 'brown_pili_rng_seeds') && ~isempty(cfg.brown_pili_rng_seeds) ...
            && isnumeric(cfg.brown_pili_rng_seeds) ...
            && all(isfinite(cfg.brown_pili_rng_seeds(:))) ...
            && all(cfg.brown_pili_rng_seeds(:) > 0)
        j.brownian.pili_rng_seeds = round(double(cfg.brown_pili_rng_seeds(:)'));
    end
end

j.output.preview_initial    = cfg.preview_initial;
j.output.generate_animation = cfg.generate_animation;
j.output.animation_fps      = cfg.animation_fps;
j.output.animation_skip     = cfg.animation_skip;
j.output.auto_confirm       = cfg.auto_confirm;
if isfield(cfg, 'compute_flowfield')
    j.output.compute_flowfield = logical(cfg.compute_flowfield);
else
    j.output.compute_flowfield = true;   % Matches run_simulation backward-compatible default
end

% Background flow
if isfield(cfg, 'background_flow')
    j.background_flow = cfg.background_flow;
else
    j.background_flow.type = 'off';
end

j.flowfield.grid_spacing_um = cfg.flowfield_grid_spacing;
j.flowfield.grid_range_um   = cfg.flowfield_grid_range;
j.flowfield.quiver_skip     = cfg.flowfield_quiver_skip;
j.flowfield.chunk_size      = cfg.flowfield_chunk_size;
if isfield(cfg, 'flowfield_stream_spacing') && isnumeric(cfg.flowfield_stream_spacing) ...
        && isscalar(cfg.flowfield_stream_spacing) && isfinite(cfg.flowfield_stream_spacing)
    j.flowfield.stream_spacing_um = cfg.flowfield_stream_spacing;
end
if isfield(cfg, 'flowfield_stream_length') && isnumeric(cfg.flowfield_stream_length) ...
        && isscalar(cfg.flowfield_stream_length) && isfinite(cfg.flowfield_stream_length)
    j.flowfield.stream_length_um = cfg.flowfield_stream_length;
end

%% Write JSON
json_str = jsonencode(j, 'PrettyPrint', true);
fid = fopen(json_path, 'w', 'n', 'UTF-8');
if fid == -1
    warning('save_sim_config:WriteError', 'Cannot write config file: %s', json_path);
    return;
end
fprintf(fid, '%s', json_str);
fclose(fid);

end
