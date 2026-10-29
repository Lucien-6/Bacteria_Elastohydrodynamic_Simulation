function compute_plot_flowfield(Z_final, params, output_dir, run_id, logf, flow_opts, t_final, silent)
%COMPUTE_PLOT_FLOWFIELD  Compute and plot flow field panels at final state
%
%  Description:
%    Computes the Stokes flow velocity field around the bacterium on the
%    XZ plane (Y = Y_cm) and YZ plane (X = X_cm) using force coefficients
%    extracted from the coupled body-pili system at the final simulation
%    state.  Produces a 2x2 panel figure:
%
%      (a) 3D view with XZ flow-field surface + 3D bacterium
%      (b) 3D view with YZ flow-field surface + 3D bacterium
%      (c) 2D XZ flow field (speed and normalised arrows)
%      (d) 2D YZ flow field (speed and normalised arrows)
%
%    This procedure has no return value.  It creates a visible figure,
%    writes flowfield_panels.png, flowfield_panels.pdf, and, when
%    supported, flowfield_panels.html in output_dir, then closes the
%    figure.  The colour scale is the 99th percentile of this frame.
%    It also reports progress and results through logf.
%
%    Evaluation uses the same hydrodynamic kernels as the simulation:
%      Body : RSM regularised Stokeslet (Cortez 2005)
%      Pili : RSS regularised Stokeslet segments + Rotlet (Cortez 2018)
%
%  Syntax:
%    compute_plot_flowfield(Z_final, params, output_dir, run_id, logf)
%    compute_plot_flowfield(Z_final, params, output_dir, run_id, logf, flow_opts)
%    compute_plot_flowfield(Z_final, params, output_dir, run_id, logf, flow_opts, t_final)
%    compute_plot_flowfield(Z_final, params, output_dir, run_id, logf, flow_opts, t_final, silent)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-06
%  Updated: 2026-10-28 — shared flow-field evaluation and panel layout

if nargin < 6 || isempty(flow_opts), flow_opts = struct(); end
if nargin < 7 || isempty(t_final), t_final = 0; end
if nargin < 8 || isempty(silent),  silent  = false; end

t_flow_start = tic;
hydro = flowfield_from_saved_state(Z_final, params, t_final, logf);
planes = flowfield_grid_velocity(hydro, params, flow_opts, logf, silent);
logf('[FLOW] Flow-field grids completed in %.1f s\n', toc(t_flow_start));

frame = pack_final_frame(hydro, planes, params, t_final, run_id);
logf('[FLOW] Generating 2x2 flow field panel...\n');

png_file  = fullfile(output_dir, 'flowfield_panels.png');
pdf_file  = fullfile(output_dir, 'flowfield_panels.pdf');
html_file = fullfile(output_dir, 'flowfield_panels.html');
draw_opts = struct('action', 'static', 'visible', 'on', ...
    'bundle_export', true, 'close_after', true, ...
    'png_path', png_file, 'pdf_path', pdf_file, 'html_path', html_file);
render_flowfield_panels(frame, draw_opts);

logf('[FLOW] Flow field panel saved: %s (600 dpi)  |  %s (vector)\n', png_file, pdf_file);
if isfile(html_file)
    logf('[FLOW] Interactive HTML: %s\n', html_file);
end

end


function frame = pack_final_frame(hydro, planes, params, t_final, run_id)
frame.planes = planes;
frame.geom = hydro.geom;
frame.params = params;
frame.L_to_um = hydro.L_to_um;
frame.cmax = frame_colour_limit(planes.samples);
frame.time_s = t_final * params.scales.T_0;
frame.run_id = run_id;
frame.scale_note = 'Frame colour scale';
end


function cmax = frame_colour_limit(samples)
samples = samples(isfinite(samples));
if isempty(samples)
    cmax = 1;
    return
end
cmax = prctile(samples, 99);
if ~(isfinite(cmax) && cmax > 0)
    cmax = max(samples);
end
if ~(isfinite(cmax) && cmax > 0)
    cmax = 1;
end
end
