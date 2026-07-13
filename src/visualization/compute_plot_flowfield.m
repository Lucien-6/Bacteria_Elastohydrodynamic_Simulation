function compute_plot_flowfield(Z_final, params, output_dir, run_id, logf, flow_opts, t_final, silent)
%COMPUTE_PLOT_FLOWFIELD  Compute and plot flow field panels at final state
%
%  Description:
%    Computes the Stokes flow velocity field around the bacterium on the
%    XZ plane (Y = Y_cm) and YZ plane (X = X_cm) using force coefficients
%    extracted from the coupled body-pili system at the final simulation
%    state.  Produces a 2x2 panel figure:
%
%      Top-left     : 3D view with XZ flow-field surface + 3D bacterium
%      Top-right    : 3D view with YZ flow-field surface + 3D bacterium
%      Bottom-left  : 2D XZ flow field (imagesc + quiver)
%      Bottom-right : 2D YZ flow field (imagesc + quiver)
%
%    This procedure has no return value.  It creates a visible figure,
%    writes flowfield_panels.png, flowfield_panels.pdf, and, when
%    supported, flowfield_panels.html in output_dir, then closes the
%    figure.  It also reports progress and results through logf.
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
%  Inputs:
%    Z_final    - [(7 + 3*sum(N_segments)) x 1] Final state vector from
%                 ODE integration
%    params     - Parameter structure (from init_coupled_params)
%    output_dir - Output directory path
%    run_id     - Run identifier string
%    logf       - Logging function handle (dual_log)
%    flow_opts  - (Optional) struct with fields:
%                   .grid_spacing_um  Grid spacing [um]       (default: 0.05)
%                   .grid_range_um    Grid extent [um]        (default: 10.0)
%                   .quiver_skip      Quiver subsample        (default: 5)
%                   .chunk_size       Eval-point chunk size   (default: 2000)
%    t_final    - (Optional) Dimensionless time associated with Z_final.
%                 Required for time-dependent active moments.
%    silent     - (Optional) logical, default false.  When true, suppresses
%                 only the chunked XZ/YZ render_progress_bar displays.
%                 It does not suppress logf messages, warnings, figure
%                 creation, file export, or any other side effect.
%
%  Outputs:
%    None.  Results are written to files as described above.
%
%  References:
%    [1] Cortez et al. (2005), Phys. Fluids 17, 031504.
%    [2] Cortez (2018), J. Comput. Phys. 375, 783-796.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-06
%  Updated: 2026-07-13

%% ---- Configuration ----
if nargin < 6 || isempty(flow_opts), flow_opts = struct(); end
if nargin < 7 || isempty(t_final), t_final = 0; end
if nargin < 8 || isempty(silent),  silent  = false; end
grid_spacing_um = get_flow_opt(flow_opts, 'grid_spacing_um', 0.05);
grid_range_um   = get_flow_opt(flow_opts, 'grid_range_um',   10.0);
quiver_skip     = get_flow_opt(flow_opts, 'quiver_skip',     5);
chunk_size      = get_flow_opt(flow_opts, 'chunk_size',      2000);

%% ---- Sanity check: params must be in sync with Z_final ----
%  Step 4 below reconstructs pili node velocities from dpsi_k/dphi_k/
%  dtheta_k using the per-pilus local basis params.rot{k}.  This basis
%  may be rotated mid-simulation by adaptive_basis_check.m when an
%  Euler-angle singularity is approached, in which case the Z vector
%  must be interpreted in the CURRENT params.rot{k} basis.  This
%  function is normally invoked on (Z_final, params) at simulation
%  end, when the two are guaranteed consistent.  We assert minimal
%  shape coherence here so that any future intermediate-time call
%  with stale params.rot triggers a clear error instead of silently
%  producing incorrect velocities.
if ~iscell(params.rot) || numel(params.rot) ~= params.M_pili
    error('compute_plot_flowfield:RotShapeMismatch', ...
        ['params.rot must be a cell array of length M_pili = %d, got %d. ', ...
         'This usually indicates params/Z desynchronisation; the function ', ...
         'expects params.rot{k} to be the local basis under which Z''s ', ...
         'pilus Euler angles are expressed.'], ...
        params.M_pili, numel(params.rot));
end
% Z layout: [X_cm(3); q(4); theta_1..N_k; phi_1..N_k; psi_1..N_k]_{k=1..M_pili}
expected_Z_len = 7 + 3*sum(params.N_segments);
if numel(Z_final) ~= expected_Z_len
    error('compute_plot_flowfield:ZShapeMismatch', ...
        'numel(Z_final) = %d does not match the expected layout 7 + 3*sum(N_segs) = %d.', ...
        numel(Z_final), expected_Z_len);
end

%% ---- Dimensional conversions ----
L_0     = params.scales.L_0;
U_0     = params.scales.U_0;
L_to_um = L_0 * 1e6;
um_to_L = 1 / L_to_um;

N_body  = params.N_body;
M_pili  = params.M_pili;
N_segs  = params.N_segments;

logf('[FLOW] Grid: %.1f x %.1f um, spacing = %.2f um, quiver skip = %d\n', ...
    grid_range_um, grid_range_um, grid_spacing_um, quiver_skip);

%% ======================================================================
%  Step 1: Reconstruct geometry from final state
%  ======================================================================
geom    = reconstruct_geometry(Z_final, params);
X_cm    = geom.X_cm;
X_all   = geom.X_all;
R_body  = geom.R_body;
N_total = size(X_all, 1);

%% ======================================================================
%  Step 2: Get state rates from ODE right-hand side
%  ======================================================================
logf('[FLOW] Computing ODE right-hand side for force extraction...\n');
t_flow_start = tic;
dZ = ode_coupled_rhs(t_final, Z_final, params);

%% ======================================================================
%  Step 3: Extract kinematic quantities from dZ
%  ======================================================================
V    = dZ(1:3);
dqdt = dZ(4:7);
q    = geom.q;

q_conj     = [q(1); -q(2:4)];
omega_quat = 2 * quat_multiply(dqdt, q_conj);
Omega      = omega_quat(2:4);

offset      = 7;
pili_dtheta = cell(M_pili, 1);
pili_dphi   = cell(M_pili, 1);
pili_dpsi   = cell(M_pili, 1);
for k = 1:M_pili
    N_k = N_segs(k);
    pili_dtheta{k} = dZ(offset + 1         : offset + N_k);
    pili_dphi{k}   = dZ(offset + N_k + 1   : offset + 2*N_k);
    pili_dpsi{k}   = dZ(offset + 2*N_k + 1 : offset + 3*N_k);
    offset = offset + 3 * N_k;
end

%% ======================================================================
%  Step 4: Compute translational velocity at all body / pili points
%  ======================================================================
v_trans = zeros(3 * N_total, 1);

for i = 1:N_body
    r_i = X_all(i,:)' - X_cm;
    v_i = V + cross(Omega, r_i);
    v_trans(3*(i-1)+1 : 3*i) = v_i;
end

pt_idx = N_body;
for k = 1:M_pili
    N_k   = N_segs(k);
    rot_k = params.rot{k};
    ds    = params.ell(k) / N_k;

    theta_k  = geom.pili_angles{k}.theta;
    phi_k    = geom.pili_angles{k}.phi;
    dtheta_k = pili_dtheta{k};
    dphi_k   = pili_dphi{k};

    v_prev = V + cross(Omega, X_all(pt_idx + 1, :)' - X_cm);
    v_trans(3*pt_idx + 1 : 3*pt_idx + 3) = v_prev;

    for j = 1:N_k
        sT = sin(theta_k(j));  cT = cos(theta_k(j));
        sP = sin(phi_k(j));    cP = cos(phi_k(j));

        de3_basis = [cT*cP*dtheta_k(j) - sT*sP*dphi_k(j);
                     cT*sP*dtheta_k(j) + sT*cP*dphi_k(j);
                     -sT*dtheta_k(j)];
        v_next = v_prev + ds * (rot_k' * de3_basis);
        v_trans(3*(pt_idx + j) + 1 : 3*(pt_idx + j) + 3) = v_next;
        v_prev = v_next;
    end
    pt_idx = pt_idx + N_k + 1;
end

%% ======================================================================
%  Step 5: Compute twist rates  (omega_twist_j = dpsi_j + cos(theta_j)*dphi_j)
%  ======================================================================
v_twist = zeros(sum(N_segs), 1);
tw_idx  = 0;
for k = 1:M_pili
    theta_k = geom.pili_angles{k}.theta;
    for j = 1:N_segs(k)
        tw_idx = tw_idx + 1;
        v_twist(tw_idx) = pili_dpsi{k}(j) + cos(theta_k(j)) * pili_dphi{k}(j);
    end
end

%% ======================================================================
%  Step 6: Build hydrodynamic matrices (re-use existing function)
%  ======================================================================
logf('[FLOW] Building hydrodynamic matrices...\n');
[M_weighted, M_twist, M_rotlet, D_twist] = build_global_resistance(X_all, params, geom);

%% ======================================================================
%  Step 7: Compute torque density and solve for force coefficients
%  ======================================================================
logf('[FLOW] Solving for force coefficients (%d DOFs)...\n', 3*N_total);
if ~isempty(M_twist)
    twist_segment_lengths = zeros(sum(N_segs), 1);
    seg_offset = 0;
    for k = 1:M_pili
        N_k = N_segs(k);
        X_k = geom.X_pili{k};
        twist_segment_lengths(seg_offset + (1:N_k)) = ...
            sqrt(sum(diff(X_k, 1, 1).^2, 2));
        seg_offset = seg_offset + N_k;
    end

    [C_twist, C_twist_reverse] = assemble_twist_coupling( ...
        M_rotlet, D_twist, params.quad_weights, twist_segment_lengths);
    Minv_v   = M_weighted \ v_trans;
    Minv_C   = M_weighted \ C_twist;
    check_linear_solve(M_weighted, [Minv_v, Minv_C], ...
        [v_trans, C_twist], 'flow-field translational solve');
    S_twist  = M_twist - C_twist_reverse * Minv_C;
    m_twist  = S_twist \ (v_twist - C_twist_reverse * Minv_v);
    check_linear_solve(S_twist, m_twist, ...
        v_twist - C_twist_reverse * Minv_v, ...
        'flow-field twist Schur solve');
    f_coeffs = Minv_v - Minv_C * m_twist;
    tau      = D_twist * m_twist;
else
    tau = zeros(0, 1);
    f_coeffs = M_weighted \ v_trans;
    check_linear_solve(M_weighted, f_coeffs, v_trans, ...
        'flow-field translational solve');
end

%% ======================================================================
%  Step 8: Precompute shared data structures for velocity evaluation
%  ======================================================================
X_cm_um     = X_cm * L_to_um;
half_range  = grid_range_um / 2;
um_s_factor = U_0 * 1e6;

% ---- Build body visualisation descriptor (dispatches parametric / custom) ----
body_vis = build_body_visualization(params, L_0);

% ---- Body kernel data ----
body_data.X  = X_all(1:N_body, :);
body_data.fx = f_coeffs(1:3:3*N_body);
body_data.fy = f_coeffs(2:3:3*N_body);
body_data.fz = f_coeffs(3:3:3*N_body);
body_data.e2 = params.eps_body^2;

% ---- Pili segment kernel data ----
seg_data = struct('active', false);
if M_pili > 0
    N_seg_all = sum(N_segs);
    seg_data.y0 = zeros(N_seg_all, 3);
    seg_data.y1 = zeros(N_seg_all, 3);
    seg_data.f0 = zeros(N_seg_all, 3);
    seg_data.f1 = zeros(N_seg_all, 3);

    s_count     = 0;
    pili_offset = N_body;
    for k = 1:M_pili
        N_k = N_segs(k);
        for s = 1:N_k
            s_count = s_count + 1;
            g0 = pili_offset + s;
            g1 = pili_offset + s + 1;
            seg_data.y0(s_count, :) = X_all(g0, :);
            seg_data.y1(s_count, :) = X_all(g1, :);
            seg_data.f0(s_count, :) = f_coeffs(3*(g0-1)+1 : 3*g0)';
            seg_data.f1(s_count, :) = f_coeffs(3*(g1-1)+1 : 3*g1)';
        end
        pili_offset = pili_offset + N_k + 1;
    end

    seg_data.v  = seg_data.y0 - seg_data.y1;
    seg_data.L  = sqrt(sum(seg_data.v.^2, 2));
    seg_data.L2 = seg_data.L.^2;
    seg_data.e2 = params.epsilon^2;
    seg_data.has_rotlet = ~isempty(tau);
    if seg_data.has_rotlet
        seg_data.tau_x = tau(1:3:end);
        seg_data.tau_y = tau(2:3:end);
        seg_data.tau_z = tau(3:3:end);
    end
    seg_data.active = true;
end

%% ======================================================================
%  Step 9: Evaluate XZ-plane flow field (Y = Y_cm)
%  ======================================================================
logf('[FLOW] Computing XZ-plane flow field...\n');

x_vec_xz = (X_cm_um(1) - half_range) : grid_spacing_um : (X_cm_um(1) + half_range);
z_vec_xz = (X_cm_um(3) - half_range) : grid_spacing_um : (X_cm_um(3) + half_range);
[XG_xz, ZG_xz] = meshgrid(x_vec_xz, z_vec_xz);
YG_xz = X_cm_um(2) * ones(size(XG_xz));

Nx_xz = length(x_vec_xz);
Nz_xz = length(z_vec_xz);

grid_xz_dl = [XG_xz(:)*um_to_L, YG_xz(:)*um_to_L, ZG_xz(:)*um_to_L];

xz_progress_tag = '[FLOW-XZ] ';
if silent
    xz_progress_tag = '';   % suppress chunked progress bar under cfg.silent
end
[Ux_xz, ~, Uz_xz, Umag_xz, mask_xz, dist_xz] = evaluate_plane_velocity( ...
    grid_xz_dl, [Nz_xz, Nx_xz], X_cm, R_body, body_vis, ...
    body_data, seg_data, um_s_factor, chunk_size, xz_progress_tag);

logf('[FLOW] XZ velocity range: [%.4e, %.4e] um/s\n', ...
    min(Umag_xz(:), [], 'omitnan'), max(Umag_xz(:), [], 'omitnan'));

%% ======================================================================
%  Step 10: Evaluate YZ-plane flow field (X = X_cm)
%  ======================================================================
logf('[FLOW] Computing YZ-plane flow field...\n');

y_vec_yz = (X_cm_um(2) - half_range) : grid_spacing_um : (X_cm_um(2) + half_range);
z_vec_yz = (X_cm_um(3) - half_range) : grid_spacing_um : (X_cm_um(3) + half_range);
[YG_yz, ZG_yz] = meshgrid(y_vec_yz, z_vec_yz);
XG_yz = X_cm_um(1) * ones(size(YG_yz));

Ny_yz = length(y_vec_yz);
Nz_yz = length(z_vec_yz);

grid_yz_dl = [XG_yz(:)*um_to_L, YG_yz(:)*um_to_L, ZG_yz(:)*um_to_L];

yz_progress_tag = '[FLOW-YZ] ';
if silent
    yz_progress_tag = '';   % suppress chunked progress bar under cfg.silent
end
[~, Uy_yz, Uz_yz, Umag_yz, mask_yz, dist_yz] = evaluate_plane_velocity( ...
    grid_yz_dl, [Nz_yz, Ny_yz], X_cm, R_body, body_vis, ...
    body_data, seg_data, um_s_factor, chunk_size, yz_progress_tag);

t_flow_elapsed = toc(t_flow_start);
logf('[FLOW] YZ velocity range: [%.4e, %.4e] um/s  (total compute: %.1f s)\n', ...
    min(Umag_yz(:), [], 'omitnan'), max(Umag_yz(:), [], 'omitnan'), t_flow_elapsed);

%% ======================================================================
%  Step 11: Build 3D body surface mesh for visualisation (mode-aware)
%  ======================================================================
body_mesh = build_body_mesh(body_vis, params, R_body, X_cm_um);

%% ======================================================================
%  Step 12: Extract body outlines for 2D subplots (mode-aware)
%  ======================================================================
[cap_xz_h, cap_xz_v] = extract_body_outline( ...
    dist_xz, [Nz_xz, Nx_xz], x_vec_xz, z_vec_xz, body_vis);
[cap_yz_h, cap_yz_v] = extract_body_outline( ...
    dist_yz, [Nz_yz, Ny_yz], y_vec_yz, z_vec_yz, body_vis);

%% ======================================================================
%  Step 13: Create 2x2 panel figure
%  ======================================================================
logf('[FLOW] Generating 2x2 flow field panel...\n');

fig_flow = figure('Position', [50, 50, 1600, 1400], 'Color', 'w', 'Visible', 'on');
set(fig_flow, 'DefaultAxesFontName', 'Arial', ...
              'DefaultAxesFontSize', 10, ...
              'DefaultAxesLineWidth', 0.8);

body_color = [0.20 0.65 0.33];
body_alpha = 0.40;

pili_colors = [0.20 0.40 0.80;
               0.85 0.33 0.10;
               0.47 0.67 0.19;
               0.93 0.69 0.13;
               0.49 0.18 0.56;
               0.30 0.75 0.93;
               0.64 0.08 0.18;
               0.00 0.45 0.74;
               0.80 0.40 0.60;
               0.50 0.50 0.50];

all_Umag = [Umag_xz(~isnan(Umag_xz)); Umag_yz(~isnan(Umag_yz))];
Umag_max = prctile(all_Umag, 99);
arrow_scale = 0.55 * grid_spacing_um * quiver_skip;
arrow_scale_2d = 1.5 * arrow_scale;

% ------------------------------------------------------------------
%  Subplot 1 (top-left): 3D view + XZ flow-field surface
% ------------------------------------------------------------------
ax1 = subplot(2, 2, 1, 'Parent', fig_flow);
hold(ax1, 'on');

alpha_xz_3d = ones(Nz_xz, Nx_xz);
alpha_xz_3d(mask_xz) = 0;
Umag_xz_s = Umag_xz;
Umag_xz_s(mask_xz) = 0;

surface(ax1, XG_xz, YG_xz, ZG_xz, Umag_xz_s, ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
    'FaceAlpha', 'texturemap', 'AlphaData', alpha_xz_3d);

render_body_mesh(ax1, body_mesh, body_color, body_alpha);

for k = 1:M_pili
    Xp   = geom.X_pili{k} * L_to_um;
    pk_c = pili_colors(mod(k-1, size(pili_colors,1)) + 1, :);
    plot3(ax1, Xp(:,1), Xp(:,2), Xp(:,3), '-', ...
          'Color', pk_c, 'LineWidth', 1.5);
end

colormap(ax1, jet(256));
caxis(ax1, [0 Umag_max]);
axis(ax1, 'equal');
xlim(ax1, X_cm_um(1) + [-half_range, half_range]);
ylim(ax1, X_cm_um(2) + [-half_range, half_range]);
zlim(ax1, X_cm_um(3) + [-half_range, half_range]);
view(ax1, 45, 25);
camlight(ax1, 'headlight');
camlight(ax1, 'left');
lighting(ax1, 'gouraud');
material(ax1, 'dull');
grid(ax1, 'on');
box(ax1, 'on');
set(ax1, 'TickDir', 'out');
xlabel(ax1, 'X (μm)', 'FontName', 'Arial', 'FontSize', 11);
ylabel(ax1, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 11);
zlabel(ax1, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
title(ax1, '3D — XZ Plane', ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
hold(ax1, 'off');

% ------------------------------------------------------------------
%  Subplot 2 (top-right): 3D view + YZ flow-field surface
% ------------------------------------------------------------------
ax2 = subplot(2, 2, 2, 'Parent', fig_flow);
hold(ax2, 'on');

alpha_yz_3d = ones(Nz_yz, Ny_yz);
alpha_yz_3d(mask_yz) = 0;
Umag_yz_s = Umag_yz;
Umag_yz_s(mask_yz) = 0;

surface(ax2, XG_yz, YG_yz, ZG_yz, Umag_yz_s, ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
    'FaceAlpha', 'texturemap', 'AlphaData', alpha_yz_3d);

render_body_mesh(ax2, body_mesh, body_color, body_alpha);

for k = 1:M_pili
    Xp   = geom.X_pili{k} * L_to_um;
    pk_c = pili_colors(mod(k-1, size(pili_colors,1)) + 1, :);
    plot3(ax2, Xp(:,1), Xp(:,2), Xp(:,3), '-', ...
          'Color', pk_c, 'LineWidth', 1.5);
end

colormap(ax2, jet(256));
caxis(ax2, [0 Umag_max]);
axis(ax2, 'equal');
xlim(ax2, X_cm_um(1) + [-half_range, half_range]);
ylim(ax2, X_cm_um(2) + [-half_range, half_range]);
zlim(ax2, X_cm_um(3) + [-half_range, half_range]);
view(ax2, 45, 25);
camlight(ax2, 'headlight');
camlight(ax2, 'left');
lighting(ax2, 'gouraud');
material(ax2, 'dull');
grid(ax2, 'on');
box(ax2, 'on');
set(ax2, 'TickDir', 'out');
xlabel(ax2, 'X (μm)', 'FontName', 'Arial', 'FontSize', 11);
ylabel(ax2, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 11);
zlabel(ax2, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
title(ax2, '3D — YZ Plane', ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
hold(ax2, 'off');

% ------------------------------------------------------------------
%  Subplot 3 (bottom-left): 2D XZ flow field
% ------------------------------------------------------------------
ax3 = subplot(2, 2, 3, 'Parent', fig_flow);
hold(ax3, 'on');

alpha_2d_xz = ones(Nz_xz, Nx_xz);
alpha_2d_xz(mask_xz) = 0;
Umag_xz_2d = Umag_xz;
Umag_xz_2d(mask_xz) = 0;

imagesc(ax3, x_vec_xz, z_vec_xz, Umag_xz_2d, 'AlphaData', alpha_2d_xz);
set(ax3, 'YDir', 'normal');
colormap(ax3, jet(256));
caxis(ax3, [0 Umag_max]);

if ~isempty(cap_xz_h)
    fill(ax3, cap_xz_h, cap_xz_v, 'w', ...
         'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1.5);
end

for k = 1:M_pili
    Xp = geom.X_pili{k} * L_to_um;
    dY = Xp(:,2) - X_cm_um(2);
    cx = [];  cz = [];
    for j = 1:size(Xp,1)-1
        if dY(j) * dY(j+1) < 0
            t  = dY(j) / (dY(j) - dY(j+1));
            cx(end+1) = Xp(j,1) + t*(Xp(j+1,1)-Xp(j,1)); %#ok<AGROW>
            cz(end+1) = Xp(j,3) + t*(Xp(j+1,3)-Xp(j,3)); %#ok<AGROW>
        end
    end
    if ~isempty(cx)
        plot(ax3, cx, cz, 'o', 'MarkerSize', 5, 'LineWidth', 1.2, ...
             'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.25 0.25 0.25]);
    end
end

q_ix = 1:quiver_skip:Nx_xz;
q_iz = 1:quiver_skip:Nz_xz;
Ux_sub = Ux_xz(q_iz, q_ix);
Uz_sub = Uz_xz(q_iz, q_ix);
XG_sub = XG_xz(q_iz, q_ix);
ZG_sub = ZG_xz(q_iz, q_ix);

spd = sqrt(Ux_sub.^2 + Uz_sub.^2);
ok  = spd > 1e-15 * max(spd(:), [], 'omitnan');
Ux_n = NaN(size(Ux_sub));  Uz_n = NaN(size(Uz_sub));
Ux_n(ok) = Ux_sub(ok) ./ spd(ok);
Uz_n(ok) = Uz_sub(ok) ./ spd(ok);

quiver(ax3, XG_sub, ZG_sub, Ux_n*arrow_scale_2d, Uz_n*arrow_scale_2d, 0, ...
    'Color', [0.12 0.12 0.12], 'LineWidth', 0.75, 'MaxHeadSize', 0.45, ...
    'AutoScale', 'off');

axis(ax3, 'equal');
xlim(ax3, X_cm_um(1) + [-half_range, half_range]);
ylim(ax3, X_cm_um(3) + [-half_range, half_range]);
xlabel(ax3, 'X (μm)', 'FontName', 'Arial', 'FontSize', 11);
ylabel(ax3, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
title(ax3, '2D — XZ Plane', ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
set(ax3, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top');
hold(ax3, 'off');

% ------------------------------------------------------------------
%  Subplot 4 (bottom-right): 2D YZ flow field
% ------------------------------------------------------------------
ax4 = subplot(2, 2, 4, 'Parent', fig_flow);
hold(ax4, 'on');

alpha_2d_yz = ones(Nz_yz, Ny_yz);
alpha_2d_yz(mask_yz) = 0;
Umag_yz_2d = Umag_yz;
Umag_yz_2d(mask_yz) = 0;

imagesc(ax4, y_vec_yz, z_vec_yz, Umag_yz_2d, 'AlphaData', alpha_2d_yz);
set(ax4, 'YDir', 'normal');
colormap(ax4, jet(256));
caxis(ax4, [0 Umag_max]);

if ~isempty(cap_yz_h)
    fill(ax4, cap_yz_h, cap_yz_v, 'w', ...
         'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1.5);
end

for k = 1:M_pili
    Xp = geom.X_pili{k} * L_to_um;
    dX = Xp(:,1) - X_cm_um(1);
    cy = [];  cz = [];
    for j = 1:size(Xp,1)-1
        if dX(j) * dX(j+1) < 0
            t  = dX(j) / (dX(j) - dX(j+1));
            cy(end+1) = Xp(j,2) + t*(Xp(j+1,2)-Xp(j,2)); %#ok<AGROW>
            cz(end+1) = Xp(j,3) + t*(Xp(j+1,3)-Xp(j,3)); %#ok<AGROW>
        end
    end
    if ~isempty(cy)
        plot(ax4, cy, cz, 'o', 'MarkerSize', 5, 'LineWidth', 1.2, ...
             'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.25 0.25 0.25]);
    end
end

q_iy  = 1:quiver_skip:Ny_yz;
q_iz2 = 1:quiver_skip:Nz_yz;
Uy_sub  = Uy_yz(q_iz2, q_iy);
Uz_sub2 = Uz_yz(q_iz2, q_iy);
YG_sub  = YG_yz(q_iz2, q_iy);
ZG_sub2 = ZG_yz(q_iz2, q_iy);

spd2 = sqrt(Uy_sub.^2 + Uz_sub2.^2);
ok2  = spd2 > 1e-15 * max(spd2(:), [], 'omitnan');
Uy_n  = NaN(size(Uy_sub));  Uz_n2 = NaN(size(Uz_sub2));
Uy_n(ok2)  = Uy_sub(ok2)  ./ spd2(ok2);
Uz_n2(ok2) = Uz_sub2(ok2) ./ spd2(ok2);

quiver(ax4, YG_sub, ZG_sub2, Uy_n*arrow_scale_2d, Uz_n2*arrow_scale_2d, 0, ...
    'Color', [0.12 0.12 0.12], 'LineWidth', 0.75, 'MaxHeadSize', 0.45, ...
    'AutoScale', 'off');

axis(ax4, 'equal');
xlim(ax4, X_cm_um(2) + [-half_range, half_range]);
ylim(ax4, X_cm_um(3) + [-half_range, half_range]);
xlabel(ax4, 'Y (μm)', 'FontName', 'Arial', 'FontSize', 11);
ylabel(ax4, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 11);
title(ax4, '2D — YZ Plane', ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
set(ax4, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top');
hold(ax4, 'off');

% ---- Shared colour bar (anchored to bottom-right subplot) ----
cb = colorbar(ax4, 'eastoutside');
cb.Label.String   = 'Velocity magnitude (μm/s)';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 11;

% ---- Global title ----
sgtitle(fig_flow, sprintf('Flow Field Panels  [%s]', run_id), ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');

drawnow;

%% ======================================================================
%  Step 14: Save
%  ======================================================================
png_file  = fullfile(output_dir, 'flowfield_panels.png');
pdf_file  = fullfile(output_dir, 'flowfield_panels.pdf');
html_file = fullfile(output_dir, 'flowfield_panels.html');
export_figure_png_pdf_html(fig_flow, png_file, pdf_file, html_file);
close(fig_flow);
logf('[FLOW] Flow field panel saved: %s (600 dpi)  |  %s (vector)\n', png_file, pdf_file);
if isfile(html_file)
    logf('[FLOW] Interactive HTML: %s\n', html_file);
end

end


%% =====================================================================
%  LOCAL FUNCTION: evaluate_plane_velocity
%  =====================================================================
function [Ux, Uy, Uz, Umag, inside_mask, dist_field] = evaluate_plane_velocity( ...
    grid_pts_dl, grid_size, X_cm, R_body, body_vis, ...
    body_data, seg_data, um_s_factor, chunk_size, progress_tag)
%EVALUATE_PLANE_VELOCITY  Stokes velocity on a planar evaluation grid
%
%  Performs body-interior masking and evaluates the RSM (body) plus
%  RSS+Rotlet (pili) velocity kernels at all exterior grid points.
%
%  Inputs:
%    grid_pts_dl   [N x 3]  Grid points in dimensionless units
%    grid_size     [1 x 2]  [Nrows, Ncols] of the grid
%    X_cm          [3 x 1]  Body CM (dimensionless)
%    R_body        [3 x 3]  Body rotation matrix
%    body_vis      struct   Body visualisation descriptor (see
%                             build_body_visualization): .type,
%                             .b_dl, .c_dl  (parametric) OR
%                             .alpha_shape_dl, .verts_local_dl, .faces (custom)
%    body_data     struct   .X, .fx, .fy, .fz, .e2
%    seg_data      struct   Pili segment data (.active, .y0, .y1, ...)
%    um_s_factor   scalar   Velocity unit conversion factor
%    chunk_size    scalar   Memory-management chunk size
%    progress_tag  char     (optional) Progress bar label, e.g. '[FLOW-XZ] '
%
%  Outputs:
%    Ux, Uy, Uz  [Nrows x Ncols]  Velocity components [um/s] (NaN inside)
%    Umag         [Nrows x Ncols]  Velocity magnitude [um/s]
%    inside_mask  [Nrows x Ncols]  Logical mask (true = interior)
%    dist_field   [N x 1]          Signed-distance-like field used for
%                                   contour extraction.  Level-set convention
%                                   depends on body_vis.type (see
%                                   extract_body_outline).

if nargin < 10 || isempty(progress_tag), progress_tag = ''; end
show_progress = ~isempty(progress_tag);

Nrows = grid_size(1);
Ncols = grid_size(2);

% ---- Body interior masking (mode-aware) ----
[inside_mask_vec, dist_field] = body_inside_mask(grid_pts_dl, X_cm, R_body, body_vis);
inside_mask = reshape(inside_mask_vec, [Nrows, Ncols]);
outside_idx = find(~inside_mask(:));
N_eval      = length(outside_idx);
X_eval      = grid_pts_dl(outside_idx, :);

% ---- Initialise output ----
ux_ev = zeros(N_eval, 1);
uy_ev = zeros(N_eval, 1);
uz_ev = zeros(N_eval, 1);

% ---- Progress bar setup ----
if show_progress
    n_chunks_body = ceil(N_eval / chunk_size);
    n_chunks_pili = 0;
    if seg_data.active
        n_chunks_pili = ceil(N_eval / chunk_size);
    end
    n_chunks_total = n_chunks_body + n_chunks_pili;
    flow_tick      = 0;
    t_flow_bar     = tic;
end

% ---- RSM body contribution ----
X_body = body_data.X;
bfx = body_data.fx;  bfy = body_data.fy;  bfz = body_data.fz;
e2_b = body_data.e2;

for ci = 1:chunk_size:N_eval
    ce  = min(ci + chunk_size - 1, N_eval);
    idx = ci:ce;

    rx = bsxfun(@minus, X_eval(idx,1), X_body(:,1)');
    ry = bsxfun(@minus, X_eval(idx,2), X_body(:,2)');
    rz = bsxfun(@minus, X_eval(idx,3), X_body(:,3)');
    r2 = rx.^2 + ry.^2 + rz.^2;

    denom = (r2 + e2_b).^1.5;
    H1    = (r2 + 2*e2_b) ./ denom;
    H2    = 1 ./ denom;
    rdotf = rx .* bfx' + ry .* bfy' + rz .* bfz';

    ux_ev(idx) = ux_ev(idx) + H1*bfx + sum(H2 .* rx .* rdotf, 2);
    uy_ev(idx) = uy_ev(idx) + H1*bfy + sum(H2 .* ry .* rdotf, 2);
    uz_ev(idx) = uz_ev(idx) + H1*bfz + sum(H2 .* rz .* rdotf, 2);

    if show_progress
        flow_tick = flow_tick + 1;
        render_progress_bar(progress_tag, flow_tick, n_chunks_total, t_flow_bar);
    end
end

% ---- RSS pili + Rotlet contribution ----
if seg_data.active
    e2_p = seg_data.e2;
    L_r  = seg_data.L';   L2_r = seg_data.L2';
    v1_r = seg_data.v(:,1)';
    v2_r = seg_data.v(:,2)';
    v3_r = seg_data.v(:,3)';

    f0x = seg_data.f0(:,1)';  f0y = seg_data.f0(:,2)';  f0z = seg_data.f0(:,3)';
    f1x = seg_data.f1(:,1)';  f1y = seg_data.f1(:,2)';  f1z = seg_data.f1(:,3)';

    for ci = 1:chunk_size:N_eval
        ce  = min(ci + chunk_size - 1, N_eval);
        idx = ci:ce;

        x01 = bsxfun(@minus, X_eval(idx,1), seg_data.y0(:,1)');
        x02 = bsxfun(@minus, X_eval(idx,2), seg_data.y0(:,2)');
        x03 = bsxfun(@minus, X_eval(idx,3), seg_data.y0(:,3)');
        r0_sq = x01.^2 + x02.^2 + x03.^2;

        x11 = x01 + v1_r;  x12 = x02 + v2_r;  x13 = x03 + v3_r;
        r1_sq = x11.^2 + x12.^2 + x13.^2;

        R0 = sqrt(r0_sq + e2_p);
        R1 = sqrt(r1_sq + e2_p);
        d0 = x01.*v1_r + x02.*v2_r + x03.*v3_r;

        H0_v = max(L_r.*R0 + d0,        1e-300);
        H1_v = max(L_r.*R1 + d0 + L2_r,  1e-300);

        inv_L  = 1 ./ L_r;
        inv_L2 = 1 ./ L2_r;
        d0_L2  = d0 .* inv_L2;

        T0m1 = (log(H1_v) - log(H0_v)) .* inv_L;
        T0m3 = 1./(R0.*H0_v) - 1./(R1.*H1_v);
        T1m1 = (R1 - R0) .* inv_L2 - d0_L2 .* T0m1;
        T1m3 = (1./R0 - 1./R1) .* inv_L2 - d0_L2 .* T0m3;
        T2m3 = -1./(R1.*L2_r) + T0m1.*inv_L2 - d0_L2 .* T1m3;
        T3m3 = -1./(R1.*L2_r) + 2*T1m1.*inv_L2 - d0_L2 .* T2m3;

        ci2 = L_r .* (T1m1 + e2_p * T1m3);
        c12 = L_r .* T1m3;
        c22 = L_r .* T2m3;
        c32 = L_r .* T3m3;

        ci1 = L_r .* (T0m1 + e2_p * T0m3) - ci2;
        c11 = L_r .* T0m3 - c12;
        c21 = L_r .* T1m3 - c22;
        c31 = L_r .* T2m3 - c32;

        % M1 * f0 (start-node kernel)
        ux_M1 = sum( ...
            (ci1+c11.*x01.*x01+c21.*2.*x01.*v1_r+c31.*v1_r.*v1_r).*f0x + ...
            (c11.*x01.*x02+c21.*(x01.*v2_r+x02.*v1_r)+c31.*v1_r.*v2_r).*f0y + ...
            (c11.*x01.*x03+c21.*(x01.*v3_r+x03.*v1_r)+c31.*v1_r.*v3_r).*f0z, 2);
        uy_M1 = sum( ...
            (c11.*x02.*x01+c21.*(x02.*v1_r+x01.*v2_r)+c31.*v2_r.*v1_r).*f0x + ...
            (ci1+c11.*x02.*x02+c21.*2.*x02.*v2_r+c31.*v2_r.*v2_r).*f0y + ...
            (c11.*x02.*x03+c21.*(x02.*v3_r+x03.*v2_r)+c31.*v2_r.*v3_r).*f0z, 2);
        uz_M1 = sum( ...
            (c11.*x03.*x01+c21.*(x03.*v1_r+x01.*v3_r)+c31.*v3_r.*v1_r).*f0x + ...
            (c11.*x03.*x02+c21.*(x03.*v2_r+x02.*v3_r)+c31.*v3_r.*v2_r).*f0y + ...
            (ci1+c11.*x03.*x03+c21.*2.*x03.*v3_r+c31.*v3_r.*v3_r).*f0z, 2);

        % M2 * f1 (end-node kernel)
        ux_M2 = sum( ...
            (ci2+c12.*x01.*x01+c22.*2.*x01.*v1_r+c32.*v1_r.*v1_r).*f1x + ...
            (c12.*x01.*x02+c22.*(x01.*v2_r+x02.*v1_r)+c32.*v1_r.*v2_r).*f1y + ...
            (c12.*x01.*x03+c22.*(x01.*v3_r+x03.*v1_r)+c32.*v1_r.*v3_r).*f1z, 2);
        uy_M2 = sum( ...
            (c12.*x02.*x01+c22.*(x02.*v1_r+x01.*v2_r)+c32.*v2_r.*v1_r).*f1x + ...
            (ci2+c12.*x02.*x02+c22.*2.*x02.*v2_r+c32.*v2_r.*v2_r).*f1y + ...
            (c12.*x02.*x03+c22.*(x02.*v3_r+x03.*v2_r)+c32.*v2_r.*v3_r).*f1z, 2);
        uz_M2 = sum( ...
            (c12.*x03.*x01+c22.*(x03.*v1_r+x01.*v3_r)+c32.*v3_r.*v1_r).*f1x + ...
            (c12.*x03.*x02+c22.*(x03.*v2_r+x02.*v3_r)+c32.*v3_r.*v2_r).*f1y + ...
            (ci2+c12.*x03.*x03+c22.*2.*x03.*v3_r+c32.*v3_r.*v3_r).*f1z, 2);

        ux_ev(idx) = ux_ev(idx) + ux_M1 + ux_M2;
        uy_ev(idx) = uy_ev(idx) + uy_M1 + uy_M2;
        uz_ev(idx) = uz_ev(idx) + uz_M1 + uz_M2;

        % Rotlet contribution.  The stored formulas are Cortez's full-curl
        % rotlet; physical torque-induced velocity uses half of this kernel.
        if seg_data.has_rotlet
            R0_3 = R0.^3;  R1_3 = R1.^3;

            d_sq     = max(L2_r.*(r0_sq + e2_p) - d0.^2, e2_p.*L2_r);
            inv_d_sq = 1 ./ d_sq;
            d1_v     = d0 + L2_r;
            inv3L2   = inv_L2 / 3;

            T0m5 = inv_d_sq/3 .* (d1_v./R1_3 - d0./R0_3) ...
                 + 2*L2_r.*inv_d_sq/3 .* T0m3;
            T1m5 = (1./R0_3 - 1./R1_3) .* inv3L2 - d0_L2 .* T0m5;

            ca  = 2*T0m3 + 3*e2_p*T0m5;
            cb  = 2*T1m3 + 3*e2_p*T1m5;
            Rxy =  L_r .* (ca.*x03 + cb.*v3_r);
            Rxz = -L_r .* (ca.*x02 + cb.*v2_r);
            Ryz =  L_r .* (ca.*x01 + cb.*v1_r);

            tx = seg_data.tau_x';
            ty = seg_data.tau_y';
            tz = seg_data.tau_z';

            ux_ev(idx) = ux_ev(idx) + 0.5 * sum( Rxy.*ty + Rxz.*tz, 2);
            uy_ev(idx) = uy_ev(idx) + 0.5 * sum(-Rxy.*tx + Ryz.*tz, 2);
            uz_ev(idx) = uz_ev(idx) + 0.5 * sum(-Rxz.*tx - Ryz.*ty, 2);
        end

        if show_progress
            flow_tick = flow_tick + 1;
            render_progress_bar(progress_tag, flow_tick, n_chunks_total, t_flow_bar);
        end
    end
end

% ---- Map to grid and convert units ----
Ux = NaN(Nrows, Ncols);
Uy = NaN(Nrows, Ncols);
Uz = NaN(Nrows, Ncols);
Ux(outside_idx) = ux_ev * um_s_factor;
Uy(outside_idx) = uy_ev * um_s_factor;
Uz(outside_idx) = uz_ev * um_s_factor;
Umag = sqrt(Ux.^2 + Uy.^2 + Uz.^2);

end


%% =====================================================================
%  LOCAL FUNCTION: build_body_visualization
%  =====================================================================
function body_vis = build_body_visualization(params, L_0)
%BUILD_BODY_VISUALIZATION  Mode-aware body descriptor for flow-field plot.
%
%  Returns a struct with the data that evaluate_plane_velocity,
%  build_body_mesh, and extract_body_outline need to render the body and
%  mask its interior.  For parametric mode, stores capsule SDF params.
%  For custom mode, stores dimensionless vertices/faces and an alphaShape
%  reconstructed in the body-local frame for inShape() masking.

    body_vis = struct('type', '');

    body_type = 'capsule';
    if isfield(params, 'body_type') && ~isempty(params.body_type)
        body_type = params.body_type;
    end
    body_vis.type = body_type;

    switch body_type
        case 'capsule'
            body_vis.b_dl = params.dim.minor_axis / L_0;
            body_vis.c_dl = params.dim.major_axis / L_0 - body_vis.b_dl;
            body_vis.a_um = params.dim.major_axis * 1e6;
            body_vis.b_um = params.dim.minor_axis * 1e6;

        case 'custom'
            if ~isfield(params, 'body_verts_local') || isempty(params.body_verts_local)
                error('compute_plot_flowfield:MissingCustomVerts', ...
                    ['Custom geometry mode expects params.body_verts_local ', ...
                     '(set in init_coupled_params).']);
            end
            body_vis.verts_local_dl = params.body_verts_local;     % [L_0]
            body_vis.faces          = params.body_faces;           % index triangulation
            body_vis.verts_m        = params.body_verts;           % [m]

            % Build alphaShape in dimensionless body-local frame (used for
            % inShape() interior masking).  If construction fails, we fall
            % back to a bounding-sphere approximation at evaluation time.
            body_vis.alpha_shape_dl = [];
            try
                shp_local = alphaShape(body_vis.verts_local_dl);
                % Prefer the critical alpha that yields a single watertight
                % region — matches the construction in custom_body_building.
                alpha_crit = criticalAlpha(shp_local, 'all-points');
                if isfinite(alpha_crit) && alpha_crit > 0
                    shp_local.Alpha = alpha_crit;
                end
                body_vis.alpha_shape_dl = shp_local;
            catch
                body_vis.alpha_shape_dl = [];
            end

            % Cache a bounding sphere for fallback masking
            centroid = mean(body_vis.verts_local_dl, 1);
            body_vis.bound_center_dl = centroid(:);
            body_vis.bound_radius_dl = max(sqrt(sum( ...
                (body_vis.verts_local_dl - centroid).^2, 2)));

        otherwise
            error('compute_plot_flowfield:UnknownBodyType', ...
                'Unknown params.body_type = ''%s''.', body_type);
    end
end


%% =====================================================================
%  LOCAL FUNCTION: body_inside_mask
%  =====================================================================
function [inside_mask_vec, dist_field] = body_inside_mask( ...
    grid_pts_dl, X_cm, R_body, body_vis)
%BODY_INSIDE_MASK  Return interior mask + SDF-like field for outline plotting.
%
%  The SDF sign convention is: values <= 0 indicate interior, values > 0
%  indicate exterior.  For parametric capsule, dist_field is the classic
%  signed distance (radius - b_dl).  For custom mode, dist_field is the
%  double() of the inside mask minus 0.5, so contourc at level 0 gives
%  the body silhouette on the plane.

    % Transform grid to body-local dimensionless frame
    X_local = (R_body.' * (grid_pts_dl - X_cm(:).').').';

    switch body_vis.type
        case 'capsule'
            x_cl = max(-body_vis.c_dl, min(body_vis.c_dl, X_local(:,1)));
            d_cap = sqrt((X_local(:,1) - x_cl).^2 + ...
                          X_local(:,2).^2 + X_local(:,3).^2);
            inside_mask_vec = d_cap <= body_vis.b_dl;
            dist_field      = d_cap - body_vis.b_dl;     % SDF: < 0 interior

        case 'custom'
            if ~isempty(body_vis.alpha_shape_dl)
                try
                    inside_mask_vec = inShape(body_vis.alpha_shape_dl, ...
                        X_local(:,1), X_local(:,2), X_local(:,3));
                catch
                    inside_mask_vec = bound_sphere_inside(X_local, body_vis);
                end
            else
                inside_mask_vec = bound_sphere_inside(X_local, body_vis);
            end
            % Binary SDF-like field centred on 0.5 -> 0 at boundary
            dist_field = double(inside_mask_vec) - 0.5;   % < 0 exterior
            dist_field = -dist_field;                     % < 0 interior (match capsule convention)

        otherwise
            error('compute_plot_flowfield:UnknownBodyType', ...
                'Unknown body_vis.type = ''%s''.', body_vis.type);
    end

    inside_mask_vec = logical(inside_mask_vec);
end


function mask_vec = bound_sphere_inside(X_local, body_vis)
    d = sqrt(sum((X_local - body_vis.bound_center_dl.').^2, 2));
    mask_vec = d <= body_vis.bound_radius_dl;
end


%% =====================================================================
%  LOCAL FUNCTION: build_body_mesh
%  =====================================================================
function body_mesh = build_body_mesh(body_vis, params, R_body, X_cm_um)
%BUILD_BODY_MESH  Construct body surface mesh in lab-frame micrometres.
%
%  Returns a struct with fields that render_body_mesh dispatches on:
%    .type = 'capsule'  : .Xsurf, .Ysurf, .Zsurf   (surf-style grid)
%    .type = 'custom'   : .X, .Y, .Z, .faces       (trisurf-style)

    body_mesh.type = body_vis.type;

    switch body_vis.type
        case 'capsule'
            a_um = body_vis.a_um;
            b_um = body_vis.b_um;

            n_circ = 36;
            n_cap  = 12;
            n_cyl  = 10;

            theta_cap = linspace(0, 2*pi, n_circ);

            phi_L = linspace(pi, pi/2, n_cap);
            [TH_L, PH_L] = meshgrid(theta_cap, phi_L);
            XL = b_um * cos(PH_L) - (a_um - b_um);
            YL = b_um * sin(PH_L) .* cos(TH_L);
            ZL = b_um * sin(PH_L) .* sin(TH_L);

            x_cyl = linspace(-(a_um - b_um), (a_um - b_um), n_cyl);
            [TH_C, XC_grid] = meshgrid(theta_cap, x_cyl);
            XC = XC_grid;
            YC = b_um * cos(TH_C);
            ZC = b_um * sin(TH_C);

            phi_R = linspace(pi/2, 0, n_cap);
            [TH_R, PH_R] = meshgrid(theta_cap, phi_R);
            XR = b_um * cos(PH_R) + (a_um - b_um);
            YR = b_um * sin(PH_R) .* cos(TH_R);
            ZR = b_um * sin(PH_R) .* sin(TH_R);

            capsule_X0 = [XL; XC; XR];
            capsule_Y0 = [YL; YC; YR];
            capsule_Z0 = [ZL; ZC; ZR];
            mesh_size  = size(capsule_X0);

            capsule_flat = [capsule_X0(:).'; capsule_Y0(:).'; capsule_Z0(:).'];
            pts_3d = R_body * capsule_flat + X_cm_um;
            body_mesh.Xsurf = reshape(pts_3d(1,:), mesh_size);
            body_mesh.Ysurf = reshape(pts_3d(2,:), mesh_size);
            body_mesh.Zsurf = reshape(pts_3d(3,:), mesh_size);

        case 'custom'
            verts_um = body_vis.verts_m * 1e6;          % body-local [um]
            pts_lab  = R_body * verts_um.' + X_cm_um;   % lab frame [um]
            body_mesh.X     = pts_lab(1, :).';
            body_mesh.Y     = pts_lab(2, :).';
            body_mesh.Z     = pts_lab(3, :).';
            body_mesh.faces = body_vis.faces;

        otherwise
            error('compute_plot_flowfield:UnknownBodyType', ...
                'Unknown body_vis.type = ''%s''.', body_vis.type);
    end

    % Suppress "unused" lint warning when body_mesh.type = 'custom' but
    % params field is unused in that branch.
    assert(isstruct(params));
end


%% =====================================================================
%  LOCAL FUNCTION: render_body_mesh
%  =====================================================================
function render_body_mesh(ax, body_mesh, body_color, body_alpha)
%RENDER_BODY_MESH  Draw the body surface on the given axes.

    switch body_mesh.type
        case 'capsule'
            surf(ax, body_mesh.Xsurf, body_mesh.Ysurf, body_mesh.Zsurf, ...
                'FaceColor', body_color, 'FaceAlpha', body_alpha, ...
                'EdgeColor', 'none', 'FaceLighting', 'gouraud');
        case 'custom'
            if isempty(body_mesh.faces) || size(body_mesh.faces, 2) ~= 3
                % Degenerate triangulation -> fall back to scatter
                scatter3(ax, body_mesh.X, body_mesh.Y, body_mesh.Z, ...
                    8, body_color, 'filled', 'MarkerFaceAlpha', body_alpha);
            else
                trisurf(body_mesh.faces, body_mesh.X, body_mesh.Y, body_mesh.Z, ...
                    'FaceColor', body_color, 'FaceAlpha', body_alpha, ...
                    'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
                    'Parent', ax);
            end
        otherwise
            error('compute_plot_flowfield:UnknownBodyType', ...
                'Unknown body_mesh.type = ''%s''.', body_mesh.type);
    end
end


%% =====================================================================
%  LOCAL FUNCTION: extract_body_outline
%  =====================================================================
function [cap_h, cap_v] = extract_body_outline(dist_field, grid_size, h_vec, v_vec, body_vis)
%EXTRACT_BODY_OUTLINE  Largest body-silhouette contour on a 2D plane.
%
%  For parametric capsule, dist_field is the classic SDF (< 0 interior),
%  so we use level 0.  For custom mode, dist_field is binary-derived;
%  using level 0 still produces a usable outline because dist_field was
%  shifted to (0.5 - mask).

    dist_grid = reshape(dist_field, grid_size);
    C = contourc(h_vec, v_vec, dist_grid, [0 0]);

    max_pts = 0;
    cap_h   = [];
    cap_v   = [];
    kk = 1;
    while kk < size(C, 2)
        n_pts = C(2, kk);
        if n_pts > max_pts
            max_pts = n_pts;
            cap_h = C(1, kk+1 : kk+n_pts);
            cap_v = C(2, kk+1 : kk+n_pts);
        end
        kk = kk + n_pts + 1;
    end

    % Suppress unused-parameter warning (body_vis reserved for future
    % custom-mode outline refinements, e.g. projected convex hull).
    if isempty(body_vis.type)
        error('compute_plot_flowfield:BadBodyVis', 'body_vis.type is empty.');
    end
end


function v = get_flow_opt(s, fname, default)
%GET_FLOW_OPT  Read field from options struct with fallback default
    if isfield(s, fname), v = s.(fname); else, v = default; end
end
