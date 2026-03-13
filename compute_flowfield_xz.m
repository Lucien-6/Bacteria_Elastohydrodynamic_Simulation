function compute_flowfield_xz(Z_final, params, output_dir, run_id, logf)
%COMPUTE_FLOWFIELD_XZ  Compute and plot XZ-plane flow field at final state
%
%  Description:
%    Computes the Stokes flow velocity field around the bacterium in the
%    XZ plane (Y = Y_cm) using force coefficients extracted from the coupled
%    body-pili system at the final simulation state.
%
%    Evaluation uses the same hydrodynamic kernels as the simulation:
%      Body : RSM regularised Stokeslet (Cortez 2005)
%      Pili : RSS regularised Stokeslet segments + Rotlet (Cortez 2018)
%
%    Output: a figure with velocity magnitude (jet colourmap) as background
%    and unit-length quiver arrows showing flow direction.  The bacterial
%    body interior is masked.
%
%  Syntax:
%    compute_flowfield_xz(Z_final, params, output_dir, run_id, logf)
%
%  Inputs:
%    Z_final    - Final state vector from ODE integration
%    params     - Parameter structure (from init_coupled_params)
%    output_dir - Output directory path
%    run_id     - Run identifier string
%    logf       - Logging function handle (dual_log)
%
%  References:
%    [1] Cortez et al. (2005), Phys. Fluids 17, 031504.
%    [2] Cortez (2018), J. Comput. Phys. 375, 783-796.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-06

%% ---- Configuration ----
grid_spacing_um = 0.05;          % Grid spacing [um]
grid_range_um   = 10.0;         % Total grid extent [um] (square domain)
quiver_skip     = 5;            % Subsample interval for quiver arrows
chunk_size      = 2000;         % Eval-point chunk size (memory management)

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
geom   = reconstruct_geometry(Z_final, params);
X_cm   = geom.X_cm;
X_all  = geom.X_all;
R_body = geom.R_body;
N_total = size(X_all, 1);

%% ======================================================================
%  Step 2: Get state rates from ODE right-hand side
%  ======================================================================
logf('[FLOW] Computing ODE right-hand side for force extraction...\n');
t_flow_start = tic;
dZ = ode_coupled_rhs(0, Z_final, params);

%% ======================================================================
%  Step 3: Extract kinematic quantities from dZ
%  ======================================================================
V    = dZ(1:3);
dqdt = dZ(4:7);
q    = geom.q;

% Recover angular velocity Omega from the quaternion rate dq/dt.
% Derivation: dq/dt = 0.5 * [0; Omega] ⊗ q
%   => [0; Omega] = 2 * (dq/dt) ⊗ q^{-1}
%                = 2 * (dq/dt) ⊗ q*          (unit quaternion: q^{-1} = q*)
% Take the vector part (components 2-4) to obtain Omega in the lab frame.
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
[M_weighted, A_twist_cells, M_rotlet] = build_global_resistance(X_all, params, geom);

A_twist_total = [];
for k = 1:M_pili
    A_twist_total = blkdiag(A_twist_total, A_twist_cells{k});
end

%% ======================================================================
%  Step 7: Compute torque density and solve for force coefficients
%  ======================================================================
if ~isempty(A_twist_total)
    tau = A_twist_total * v_twist;
    if ~isempty(M_rotlet)
        v_trans_corr = v_trans - M_rotlet * tau;
    else
        v_trans_corr = v_trans;
    end
else
    tau = zeros(0, 1);
    v_trans_corr = v_trans;
end

logf('[FLOW] Solving for force coefficients (%d DOFs)...\n', 3*N_total);
f_coeffs = M_weighted \ v_trans_corr;

%% ======================================================================
%  Step 8: Build evaluation grid (XZ plane, centred on body CM)
%  ======================================================================
X_cm_um    = X_cm * L_to_um;
half_range = grid_range_um / 2;

x_vec_um = (X_cm_um(1) - half_range) : grid_spacing_um : (X_cm_um(1) + half_range);
z_vec_um = (X_cm_um(3) - half_range) : grid_spacing_um : (X_cm_um(3) + half_range);
[XG_um, ZG_um] = meshgrid(x_vec_um, z_vec_um);
YG_um = X_cm_um(2) * ones(size(XG_um));

Nx = length(x_vec_um);
Nz = length(z_vec_um);

XG = XG_um * um_to_L;
YG = YG_um * um_to_L;
ZG = ZG_um * um_to_L;

logf('[FLOW] Evaluation grid: %d x %d = %d points\n', Nx, Nz, Nx*Nz);

%% ======================================================================
%  Step 9: Mask interior points (capsule body)
%  ======================================================================
b_dl = params.dim.minor_axis / L_0;
c_dl = params.dim.major_axis / L_0 - b_dl;

X_grid_lab   = [XG(:), YG(:), ZG(:)];
X_grid_local = (R_body' * (X_grid_lab - X_cm(:)').')';

x_clamped    = max(-c_dl, min(c_dl, X_grid_local(:,1)));
dist_capsule = sqrt((X_grid_local(:,1) - x_clamped).^2 + ...
                     X_grid_local(:,2).^2 + X_grid_local(:,3).^2);
inside_mask  = reshape(dist_capsule <= b_dl, [Nz, Nx]);

outside_idx = find(~inside_mask(:));
N_eval      = length(outside_idx);
X_eval      = [XG(outside_idx), YG(outside_idx), ZG(outside_idx)];

logf('[FLOW] Interior points masked: %d / %d\n', Nx*Nz - N_eval, Nx*Nz);

%% ======================================================================
%  Step 10: Evaluate velocity at exterior grid points
%  ======================================================================
logf('[FLOW] Evaluating velocity field at %d points...\n', N_eval);

ux_eval = zeros(N_eval, 1);
uy_eval = zeros(N_eval, 1);
uz_eval = zeros(N_eval, 1);

% ---- 10a: RSM body contribution ----
e2_body = params.eps_body^2;
X_body  = X_all(1:N_body, :);
fx_body = f_coeffs(1:3:3*N_body);
fy_body = f_coeffs(2:3:3*N_body);
fz_body = f_coeffs(3:3:3*N_body);

for ci = 1:chunk_size:N_eval
    ci_end = min(ci + chunk_size - 1, N_eval);
    idx    = ci:ci_end;

    rx = bsxfun(@minus, X_eval(idx,1), X_body(:,1)');
    ry = bsxfun(@minus, X_eval(idx,2), X_body(:,2)');
    rz = bsxfun(@minus, X_eval(idx,3), X_body(:,3)');
    r2 = rx.^2 + ry.^2 + rz.^2;

    denom = (r2 + e2_body).^1.5;
    H1    = (r2 + 2*e2_body) ./ denom;
    H2    = 1 ./ denom;

    rdotf = rx .* fx_body' + ry .* fy_body' + rz .* fz_body';

    ux_eval(idx) = ux_eval(idx) + H1 * fx_body + sum(H2 .* rx .* rdotf, 2);
    uy_eval(idx) = uy_eval(idx) + H1 * fy_body + sum(H2 .* ry .* rdotf, 2);
    uz_eval(idx) = uz_eval(idx) + H1 * fz_body + sum(H2 .* rz .* rdotf, 2);
end

% ---- 10b + 10c: RSS pili + Rotlet contribution ----
if M_pili > 0
    eps_pili  = params.epsilon;
    e2_pili   = eps_pili^2;
    N_seg_all = sum(N_segs);

    y0_seg = zeros(N_seg_all, 3);
    y1_seg = zeros(N_seg_all, 3);
    f0_seg = zeros(N_seg_all, 3);
    f1_seg = zeros(N_seg_all, 3);

    s_count     = 0;
    pili_offset = N_body;
    for k = 1:M_pili
        N_k = N_segs(k);
        for s = 1:N_k
            s_count = s_count + 1;
            g0 = pili_offset + s;
            g1 = pili_offset + s + 1;
            y0_seg(s_count, :) = X_all(g0, :);
            y1_seg(s_count, :) = X_all(g1, :);
            f0_seg(s_count, :) = f_coeffs(3*(g0-1)+1 : 3*g0)';
            f1_seg(s_count, :) = f_coeffs(3*(g1-1)+1 : 3*g1)';
        end
        pili_offset = pili_offset + N_k + 1;
    end

    v_seg  = y0_seg - y1_seg;
    L_seg  = sqrt(sum(v_seg.^2, 2));
    L2_seg = L_seg.^2;

    L_r  = L_seg';    L2_r = L2_seg';
    v1_r = v_seg(:,1)';  v2_r = v_seg(:,2)';  v3_r = v_seg(:,3)';

    has_rotlet = ~isempty(tau);
    if has_rotlet
        tau_x = tau(1:3:end);
        tau_y = tau(2:3:end);
        tau_z = tau(3:3:end);
    end

    for ci = 1:chunk_size:N_eval
        ci_end = min(ci + chunk_size - 1, N_eval);
        idx    = ci:ci_end;

        x01 = bsxfun(@minus, X_eval(idx,1), y0_seg(:,1)');
        x02 = bsxfun(@minus, X_eval(idx,2), y0_seg(:,2)');
        x03 = bsxfun(@minus, X_eval(idx,3), y0_seg(:,3)');

        r0_sq = x01.^2 + x02.^2 + x03.^2;

        x11 = x01 + v1_r;
        x12 = x02 + v2_r;
        x13 = x03 + v3_r;
        r1_sq = x11.^2 + x12.^2 + x13.^2;

        R0 = sqrt(r0_sq + e2_pili);
        R1 = sqrt(r1_sq + e2_pili);
        d0 = x01.*v1_r + x02.*v2_r + x03.*v3_r;

        H0_v = L_r.*R0 + d0;
        H1_v = L_r.*R1 + d0 + L2_r;
        H0_v = max(H0_v, 1e-300);
        H1_v = max(H1_v, 1e-300);

        inv_L  = 1 ./ L_r;
        inv_L2 = 1 ./ L2_r;
        d0_L2  = d0 .* inv_L2;

        T0m1 = (log(H1_v) - log(H0_v)) .* inv_L;
        T0m3 = 1./(R0.*H0_v) - 1./(R1.*H1_v);
        T1m1 = (R1 - R0) .* inv_L2 - d0_L2 .* T0m1;
        T1m3 = (1./R0 - 1./R1) .* inv_L2 - d0_L2 .* T0m3;
        T2m3 = -1./(R1.*L2_r) + T0m1.*inv_L2 - d0_L2 .* T1m3;
        T3m3 = -1./(R1.*L2_r) + 2*T1m1.*inv_L2 - d0_L2 .* T2m3;

        % M2 coefficients (end-node kernel)
        ci2 = L_r .* (T1m1 + e2_pili * T1m3);
        c12 = L_r .* T1m3;
        c22 = L_r .* T2m3;
        c32 = L_r .* T3m3;

        % M1 = T0-block - M2 coefficients (start-node kernel)
        ci1 = L_r .* (T0m1 + e2_pili * T0m3) - ci2;
        c11 = L_r .* T0m3                     - c12;
        c21 = L_r .* T1m3                     - c22;
        c31 = L_r .* T2m3                     - c32;

        % ---- M1 * f0  (start node contribution) ----
        ux_M1 = sum( ...
            (ci1 + c11.*x01.*x01 + c21.*2.*x01.*v1_r + c31.*v1_r.*v1_r) .* f0_seg(:,1)' + ...
            (      c11.*x01.*x02 + c21.*(x01.*v2_r+x02.*v1_r) + c31.*v1_r.*v2_r) .* f0_seg(:,2)' + ...
            (      c11.*x01.*x03 + c21.*(x01.*v3_r+x03.*v1_r) + c31.*v1_r.*v3_r) .* f0_seg(:,3)', 2);

        uy_M1 = sum( ...
            (      c11.*x02.*x01 + c21.*(x02.*v1_r+x01.*v2_r) + c31.*v2_r.*v1_r) .* f0_seg(:,1)' + ...
            (ci1 + c11.*x02.*x02 + c21.*2.*x02.*v2_r + c31.*v2_r.*v2_r) .* f0_seg(:,2)' + ...
            (      c11.*x02.*x03 + c21.*(x02.*v3_r+x03.*v2_r) + c31.*v2_r.*v3_r) .* f0_seg(:,3)', 2);

        uz_M1 = sum( ...
            (      c11.*x03.*x01 + c21.*(x03.*v1_r+x01.*v3_r) + c31.*v3_r.*v1_r) .* f0_seg(:,1)' + ...
            (      c11.*x03.*x02 + c21.*(x03.*v2_r+x02.*v3_r) + c31.*v3_r.*v2_r) .* f0_seg(:,2)' + ...
            (ci1 + c11.*x03.*x03 + c21.*2.*x03.*v3_r + c31.*v3_r.*v3_r) .* f0_seg(:,3)', 2);

        % ---- M2 * f1  (end node contribution) ----
        ux_M2 = sum( ...
            (ci2 + c12.*x01.*x01 + c22.*2.*x01.*v1_r + c32.*v1_r.*v1_r) .* f1_seg(:,1)' + ...
            (      c12.*x01.*x02 + c22.*(x01.*v2_r+x02.*v1_r) + c32.*v1_r.*v2_r) .* f1_seg(:,2)' + ...
            (      c12.*x01.*x03 + c22.*(x01.*v3_r+x03.*v1_r) + c32.*v1_r.*v3_r) .* f1_seg(:,3)', 2);

        uy_M2 = sum( ...
            (      c12.*x02.*x01 + c22.*(x02.*v1_r+x01.*v2_r) + c32.*v2_r.*v1_r) .* f1_seg(:,1)' + ...
            (ci2 + c12.*x02.*x02 + c22.*2.*x02.*v2_r + c32.*v2_r.*v2_r) .* f1_seg(:,2)' + ...
            (      c12.*x02.*x03 + c22.*(x02.*v3_r+x03.*v2_r) + c32.*v2_r.*v3_r) .* f1_seg(:,3)', 2);

        uz_M2 = sum( ...
            (      c12.*x03.*x01 + c22.*(x03.*v1_r+x01.*v3_r) + c32.*v3_r.*v1_r) .* f1_seg(:,1)' + ...
            (      c12.*x03.*x02 + c22.*(x03.*v2_r+x02.*v3_r) + c32.*v3_r.*v2_r) .* f1_seg(:,2)' + ...
            (ci2 + c12.*x03.*x03 + c22.*2.*x03.*v3_r + c32.*v3_r.*v3_r) .* f1_seg(:,3)', 2);

        ux_eval(idx) = ux_eval(idx) + ux_M1 + ux_M2;
        uy_eval(idx) = uy_eval(idx) + uy_M1 + uy_M2;
        uz_eval(idx) = uz_eval(idx) + uz_M1 + uz_M2;

        % ---- Rotlet contribution (Cortez 2018 Sec. 2.5) ----
        if has_rotlet
            R0_3 = R0.^3;
            R1_3 = R1.^3;

            d_sq     = L2_r .* (r0_sq + e2_pili) - d0.^2;
            d_sq     = max(d_sq, e2_pili .* L2_r);
            inv_d_sq = 1 ./ d_sq;
            d1_v     = d0 + L2_r;
            inv3L2   = inv_L2 / 3;

            T0m5 = inv_d_sq/3 .* (d1_v./R1_3 - d0./R0_3) ...
                 + 2*L2_r.*inv_d_sq/3 .* T0m3;
            T1m5 = (1./R0_3 - 1./R1_3) .* inv3L2 - d0_L2 .* T0m5;

            ca  = 2*T0m3 + 3*e2_pili*T0m5;
            cb  = 2*T1m3 + 3*e2_pili*T1m5;
            Rxy =  L_r .* (ca.*x03 + cb.*v3_r);
            Rxz = -L_r .* (ca.*x02 + cb.*v2_r);
            Ryz =  L_r .* (ca.*x01 + cb.*v1_r);

            ux_eval(idx) = ux_eval(idx) + sum( Rxy.*tau_y' + Rxz.*tau_z', 2);
            uy_eval(idx) = uy_eval(idx) + sum(-Rxy.*tau_x' + Ryz.*tau_z', 2);
            uz_eval(idx) = uz_eval(idx) + sum(-Rxz.*tau_x' - Ryz.*tau_y', 2);
        end
    end
end

%% ======================================================================
%  Step 11: Map velocity to grid and convert to physical units [um/s]
%  ======================================================================
um_s_factor = U_0 * 1e6;

Ux_um = NaN(Nz, Nx);
Uy_um = NaN(Nz, Nx);
Uz_um = NaN(Nz, Nx);

Ux_um(outside_idx) = ux_eval * um_s_factor;
Uy_um(outside_idx) = uy_eval * um_s_factor;
Uz_um(outside_idx) = uz_eval * um_s_factor;

Umag = sqrt(Ux_um.^2 + Uy_um.^2 + Uz_um.^2);

t_flow_elapsed = toc(t_flow_start);
logf('[FLOW] Velocity range: [%.4e, %.4e] um/s  (computed in %.1f s)\n', ...
    min(Umag(:), [], 'omitnan'), max(Umag(:), [], 'omitnan'), t_flow_elapsed);

%% ======================================================================
%  Step 12: Extract capsule outline from distance field
%  ======================================================================
dist_grid = reshape(dist_capsule, [Nz, Nx]);

C_cap   = contourc(x_vec_um, z_vec_um, dist_grid, [b_dl b_dl]);
max_pts = 0;
cap_x   = [];
cap_z   = [];
kk = 1;
while kk < size(C_cap, 2)
    n_pts = C_cap(2, kk);
    if n_pts > max_pts
        max_pts = n_pts;
        cap_x = C_cap(1, kk+1 : kk+n_pts);
        cap_z = C_cap(2, kk+1 : kk+n_pts);
    end
    kk = kk + n_pts + 1;
end

%% ======================================================================
%  Step 13: Plot
%  ======================================================================
logf('[FLOW] Generating flow field figure...\n');

fig_flow = figure('Position', [100, 100, 900, 800], 'Color', 'w', 'Visible', 'on');
set(fig_flow, 'DefaultAxesFontName', 'Arial', ...
              'DefaultAxesFontSize', 11, ...
              'DefaultAxesLineWidth', 0.8);

ax_flow = axes(fig_flow);
hold(ax_flow, 'on');

% ---- Background: velocity magnitude (jet) ----
alpha_data = ones(Nz, Nx);
alpha_data(inside_mask) = 0;

Umag_plot = Umag;
Umag_plot(inside_mask) = 0;

imagesc(ax_flow, x_vec_um, z_vec_um, Umag_plot, 'AlphaData', alpha_data);
set(ax_flow, 'YDir', 'normal');
colormap(ax_flow, jet(256));
caxis(ax_flow, [0 max(Umag(:), [], 'omitnan')]);

cb = colorbar(ax_flow);
cb.Label.String   = 'Velocity magnitude (μm/s)';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 12;

% ---- Capsule body fill (white) + outline ----
if ~isempty(cap_x)
    fill(ax_flow, cap_x, cap_z, 'w', ...
         'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1.5);
end

% ---- Pili overlay ----
pili_palette = [0.85 0.85 0.85;
                0.75 0.75 0.75;
                0.90 0.90 0.90;
                0.70 0.70 0.70;
                0.80 0.80 0.80];
for k = 1:M_pili
    Xp   = geom.X_pili{k} * L_to_um;
    pk_c = pili_palette(mod(k-1, size(pili_palette,1)) + 1, :);
    plot(ax_flow, Xp(:,1), Xp(:,3), '-', 'Color', pk_c, 'LineWidth', 2.0);
end

% ---- Quiver: unit-length direction arrows (subsampled) ----
q_idx_x = 1:quiver_skip:Nx;
q_idx_z = 1:quiver_skip:Nz;

Ux_sub = Ux_um(q_idx_z, q_idx_x);
Uz_sub = Uz_um(q_idx_z, q_idx_x);
XG_sub = XG_um(q_idx_z, q_idx_x);
ZG_sub = ZG_um(q_idx_z, q_idx_x);

speed_in = sqrt(Ux_sub.^2 + Uz_sub.^2);
valid    = speed_in > 1e-15 * max(speed_in(:), [], 'omitnan');
Ux_norm  = NaN(size(Ux_sub));
Uz_norm  = NaN(size(Uz_sub));
Ux_norm(valid) = Ux_sub(valid) ./ speed_in(valid);
Uz_norm(valid) = Uz_sub(valid) ./ speed_in(valid);

arrow_scale = 0.55 * grid_spacing_um * quiver_skip;

% White halo for arrow visibility on dark backgrounds
quiver(ax_flow, XG_sub, ZG_sub, ...
    Ux_norm * arrow_scale, Uz_norm * arrow_scale, 0, ...
    'Color', 'w', 'LineWidth', 1.2, 'MaxHeadSize', 0.4, ...
    'AutoScale', 'off');

quiver(ax_flow, XG_sub, ZG_sub, ...
    Ux_norm * arrow_scale, Uz_norm * arrow_scale, 0, ...
    'Color', [0.12 0.12 0.12], 'LineWidth', 0.5, 'MaxHeadSize', 0.3, ...
    'AutoScale', 'off');

% ---- Formatting ----
axis(ax_flow, 'equal');
xlim(ax_flow, [X_cm_um(1) - half_range, X_cm_um(1) + half_range]);
ylim(ax_flow, [X_cm_um(3) - half_range, X_cm_um(3) + half_range]);
xlabel(ax_flow, 'X (μm)', 'FontName', 'Arial', 'FontSize', 13);
ylabel(ax_flow, 'Z (μm)', 'FontName', 'Arial', 'FontSize', 13);
title(ax_flow, sprintf('Flow Field — XZ Plane  [%s]', run_id), ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
set(ax_flow, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top');

hold(ax_flow, 'off');
drawnow;

%% ======================================================================
%  Step 14: Save
%  ======================================================================
png_file = fullfile(output_dir, 'flowfield_xz.png');
exportgraphics(fig_flow, png_file, 'Resolution', 600);
logf('[FLOW] Flow field saved: %s (600 dpi)\n', png_file);

end
