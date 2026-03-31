function save_sim_config(cfg, json_path)
%SAVE_SIM_CONFIG  Save current simulation configuration to a JSON file
%
%  Description:
%    Writes the full set of simulation parameters to a human-readable JSON
%    file for reproducibility.  The file can later be loaded with
%    load_sim_config to reproduce the exact same run.
%
%  Syntax:
%    save_sim_config(cfg, json_path)
%
%  Inputs:
%    cfg       - Flat configuration struct (as used in main_coupled_simulation)
%    json_path - Destination file path (e.g. 'output/MyCase/config.json')
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17

%% Build nested structure matching config_example.json layout
j.case_name = cfg.case_name;

j.body.minor_axis = cfg.minor_axis;
j.body.major_axis = cfg.major_axis;
j.body.N_body     = cfg.N_body;
j.body.alpha_0    = cfg.alpha_0;

% Pili: convert [8 x M] matrix back to struct array
M = size(cfg.Pili_Matrix, 2);
type_names = {'linear', 'parabolic', 'circular', 'sinusoidal', 'spiral'};
pili_arr = [];
for k = 1:M
    p_entry.type        = type_names{cfg.Pili_Matrix(1, k)};
    p_entry.length_um   = cfg.Pili_Matrix(2, k);
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
    if isempty(pili_arr)
        pili_arr = p_entry;
    else
        pili_arr(end+1) = p_entry; %#ok<AGROW>
    end
end
if isempty(pili_arr)
    j.pili = [];
else
    j.pili = pili_arr;
end

j.material.E       = cfg.E;
j.material.r_pili  = cfg.r_pili;
j.material.sigma_p = cfg.sigma_p;

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
if isfield(cfg, 'fixman_midpoint')
    j.brownian.fixman_midpoint = cfg.fixman_midpoint;
else
    j.brownian.fixman_midpoint = false;
end

j.output.preview_initial    = cfg.preview_initial;
j.output.generate_animation = cfg.generate_animation;
j.output.animation_fps      = cfg.animation_fps;
j.output.animation_skip     = cfg.animation_skip;
j.output.auto_confirm       = cfg.auto_confirm;

j.flowfield.grid_spacing_um = cfg.flowfield_grid_spacing;
j.flowfield.grid_range_um   = cfg.flowfield_grid_range;
j.flowfield.quiver_skip     = cfg.flowfield_quiver_skip;
j.flowfield.chunk_size      = cfg.flowfield_chunk_size;

%% Write JSON
json_str = jsonencode(j, 'PrettyPrint', true);
fid = fopen(json_path, 'w', 'n', 'UTF-8');
if fid == -1
    warning('save_sim_config:WriteError', 'Cannot write config file: %s', json_path);
    return;
end
fwrite(fid, json_str, 'char');
fclose(fid);

end
