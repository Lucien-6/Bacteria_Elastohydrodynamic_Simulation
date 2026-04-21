function [lin_sys, geom, hydro] = build_coupled_lin_sys_v0(Z, params)
%BUILD_COUPLED_LIN_SYS_V0  Frozen reference copy (pre-optimization) used by test_optimization_invariance
%
%  Description:
%    Builds the full coupled linear system for the body-pili
%    elastohydrodynamic problem:
%      lin_sys * z_dot = R
%    where lin_sys encodes the hydrodynamic resistance (mobility
%    inversion + force/moment balance) and clamped boundary conditions.
%
%    This function extracts Steps 1-7 of ode_coupled_rhs.m into a
%    reusable subroutine, shared by ode_coupled_rhs (deterministic
%    RHS) and compute_effective_body_brownian (Brownian diffusion).
%
%    The elastic restoring moments (RHS) and state-rate packing are
%    NOT included here — those remain in ode_coupled_rhs.m.
%
%  Syntax:
%    [lin_sys, geom] = build_coupled_lin_sys(Z, params)
%    [lin_sys, geom, hydro] = build_coupled_lin_sys(Z, params)
%
%  Inputs:
%    Z      - [n x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    lin_sys - [n_state x n_state] Coupled linear system matrix
%              (with clamped BC applied).  Represents -Xi_gen in
%              dimensionless units.
%    geom    - Geometry struct from reconstruct_geometry (needed by
%              caller for elastic moments, director frames, etc.)
%    hydro   - (Optional, 3rd output) Hydrodynamic intermediate struct
%              for background flow computation.  Fields:
%                .B_global, .M_weighted, .A_twist_total, .M_rotlet, .n_trans
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] Cortez et al. (2005), Phys. Fluids 17, 031504.
%    [3] Cortez (2018), J. Comput. Phys. 375, 783-796.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-31

N_body = params.N_body;
M      = params.M_pili;
N_segs = params.N_segments;

%% Step 1: Reconstruct geometry
geom    = reconstruct_geometry(Z, params);
X_all   = geom.X_all;
X_cm    = geom.X_cm;
N_total = size(X_all, 1);

%% Step 2: Build global mobility / resistance matrices
[M_weighted, A_twist_cells, M_rotlet] = build_global_resistance_v0(X_all, params, geom);

%% Step 3-4: Build B, Q matrices
n_state = 6 + 3 * sum(N_segs);
n_vel   = 3 * N_total + sum(N_segs);
n_force = 3 * N_total + 3 * sum(N_segs);

%% Build Q_global: [n_vel x n_state]
Q_global = zeros(n_vel, n_state);

for i = 1:N_body
    ii  = 3*(i-1) + 1;
    r_i = X_all(i,:)' - X_cm;
    Q_global(ii:ii+2, 1:3) = eye(3);
    Q_global(ii:ii+2, 4:6) = -skew_matrix(r_i);
end

pt_offset  = N_body;
col_offset = 6;
twist_row  = 3 * N_total;

for k = 1:M
    N_k   = N_segs(k);
    rot_k = params.rot{k};
    ell_k = params.ell(k);

    theta_k = geom.pili_angles{k}.theta;
    phi_k   = geom.pili_angles{k}.phi;

    r_root = X_all(pt_offset + 1, :)' - X_cm;

    Q_k_spatial = zeros(3*(N_k+1), 2*N_k + 3);

    Q_k_spatial(1:N_k+1, 1)                     = 1;
    Q_k_spatial(N_k+2 : 2*(N_k+1), 2)           = 1;
    Q_k_spatial(2*(N_k+1)+1 : 3*(N_k+1), 3)     = 1;

    Q12 = zeros(N_k+1, N_k);  Q13 = zeros(N_k+1, N_k);
    Q22 = zeros(N_k+1, N_k);  Q23 = zeros(N_k+1, N_k);
    Q32 = zeros(N_k+1, N_k);

    for j = 1:N_k
        cT = cos(theta_k(j));  sT = sin(theta_k(j));
        cP = cos(phi_k(j));    sP = sin(phi_k(j));
        Q12(j+1:end, j) =  cT * cP;
        Q13(j+1:end, j) = -sT * sP;
        Q22(j+1:end, j) =  cT * sP;
        Q23(j+1:end, j) =  sT * cP;
        Q32(j+1:end, j) = -sT;
    end

    Q_k_spatial(:, 4:end) = [Q12, Q13; Q22, Q23; Q32, zeros(N_k+1,N_k)] ...
                            * (ell_k / N_k);

    row_order = zeros(3*(N_k+1), 1);
    for j = 1:N_k+1
        row_order(3*(j-1)+1 : 3*j) = [j, N_k+1+j, 2*(N_k+1)+j];
    end
    Q_k_spatial = Q_k_spatial(row_order, :);

    Q_k_ext = blkdiag(Q_k_spatial, eye(N_k));
    Q_k_ext(3*(N_k+1)+1 : 3*(N_k+1)+N_k, 3+N_k+1 : 3+2*N_k) = diag(cos(theta_k));

    rot_k_inv    = rot_k';
    R_block_trans = kron(eye(N_k+1), rot_k_inv);
    R_block_full  = blkdiag(R_block_trans, eye(N_k));
    Q_k_lab       = R_block_full * Q_k_ext;

    J_base = zeros(3, 6);
    J_base(1:3, 1:3) = eye(3);
    J_base(1:3, 4:6) = -skew_matrix(r_root);
    Q_k_body = Q_k_lab(:, 1:3) * rot_k * J_base;

    trans_rows = 3*pt_offset + (1 : 3*(N_k+1));
    twist_rows = twist_row + (1 : N_k);
    all_rows   = [trans_rows, twist_rows];
    body_cols  = 1:6;
    pili_cols  = col_offset + (1 : 3*N_k);

    Q_global(all_rows, body_cols) = Q_k_body;
    Q_global(all_rows, pili_cols) = Q_k_lab(:, 4:end);

    pt_offset  = pt_offset + N_k + 1;
    col_offset = col_offset + 3 * N_k;
    twist_row  = twist_row + N_k;
end

%% Build torsional resistance block
A_twist_total = [];
for k = 1:M
    A_twist_total = blkdiag(A_twist_total, A_twist_cells{k});
end

%% Build B_global: [n_state x n_force]
B_global = zeros(n_state, n_force);

for i = 1:N_body
    ii = 3*(i-1) + 1;
    B_global(1, ii)   = 1;
    B_global(2, ii+1) = 1;
    B_global(3, ii+2) = 1;
end

pt_offset = N_body;
for k = 1:M
    N_k   = N_segs(k);
    ell_k = params.ell(k);
    temp  = zeros(3*(N_k+1), 1);
    temp(1:3:end) = 2;
    temp(1)       = 1;
    temp(end-2)   = 1;
    temp = temp * ell_k / (2 * N_k);

    base_col = 3 * pt_offset;
    B_global(1, base_col+1 : base_col+3*(N_k+1)) = temp';
    B_global(2, base_col+1 : base_col+3*(N_k+1)) = circshift(temp, 1)';
    B_global(3, base_col+1 : base_col+3*(N_k+1)) = circshift(temp, 2)';

    pt_offset = pt_offset + N_k + 1;
end

for i = 1:N_total
    r_i = X_all(i,:)' - X_cm;
    ii  = 3*(i-1) + 1;

    if i <= N_body
        w_i = 1;
    else
        w_i = get_trapezoidal_weight(i, N_body, N_segs, params.ell);
    end

    B_global(4, ii:ii+2) = w_i * [0, -r_i(3),  r_i(2)];
    B_global(5, ii:ii+2) = w_i * [r_i(3),  0, -r_i(1)];
    B_global(6, ii:ii+2) = w_i * [-r_i(2), r_i(1),  0];
end

row_offset       = 6;
col_offset_trans = 3 * N_body;
col_offset_twist = 3 * N_total;

for k = 1:M
    N_k    = N_segs(k);
    Nk2    = N_k^2;
    rot_k  = params.rot{k};
    ell_k_B = params.ell(k);

    X_k_lab = geom.X_pili{k};
    X_k     = (rot_k * X_k_lab')' / ell_k_B;

    d1_k = geom.d1_basis{k};
    d2_k = geom.d2_basis{k};
    d3_k = geom.d3_basis{k};

    [d1N, d2N, d3N] = compute_node_directors(d1_k, d2_k, d3_k, N_k);

    B_k   = zeros(3*N_k, 3*(N_k+1) + 3*N_k);
    inter = (d3_k(2:end,:) + 2*d3_k(1:end-1,:)) / (6*Nk2) + ...
            (X_k(2:end-1,:) + X_k(1:end-2,:)) / (2*N_k);

    for i = 1:N_k
        i_row = 3*(i-1) + 1;
        j = i;  j_col = 3*(j-1);
        B_k(i_row, j_col+1:j_col+3) = cross(d1N(i,:), ...
            (X_k(j,:) - X_k(i,:))/(2*N_k) + d3_k(j,:)/(6*Nk2));
        if i < N_k
            K_mat = [0, -d1N(i,3), d1N(i,2); d1N(i,3), 0, -d1N(i,1); ...
                     -d1N(i,2), d1N(i,1), 0];
            d_mat   = inter(i:end,:) - X_k(i,:)/N_k;
            crosses = K_mat * d_mat';
            B_k(i_row, 3*i+1:3*N_k) = crosses(:)';
        end
        j = N_k + 1;  j_col = 3*(j-1);
        B_k(i_row, j_col+1:j_col+3) = cross(d1N(i,:), ...
            (X_k(j-1,:) - X_k(i,:))/(2*N_k) + d3_k(j-1,:)/(3*Nk2));
        j_col_tau = 3*(N_k+1) + 3*(i-1) + 1;
        for j = i:N_k
            B_k(i_row, j_col_tau)   = d1N(i,1) / N_k;
            B_k(i_row, j_col_tau+1) = d1N(i,2) / N_k;
            B_k(i_row, j_col_tau+2) = d1N(i,3) / N_k;
            j_col_tau = j_col_tau + 3;
        end
    end

    for i = 1:N_k
        i_row = 3*(i-1) + 2;
        j = i;  j_col = 3*(j-1);
        B_k(i_row, j_col+1:j_col+3) = cross(d2N(i,:), ...
            (X_k(j,:) - X_k(i,:))/(2*N_k) + d3_k(j,:)/(6*Nk2));
        if i < N_k
            K_mat = [0, -d2N(i,3), d2N(i,2); d2N(i,3), 0, -d2N(i,1); ...
                     -d2N(i,2), d2N(i,1), 0];
            d_mat   = inter(i:end,:) - X_k(i,:)/N_k;
            crosses = K_mat * d_mat';
            B_k(i_row, 3*i+1:3*N_k) = crosses(:)';
        end
        j = N_k + 1;  j_col = 3*(j-1);
        B_k(i_row, j_col+1:j_col+3) = cross(d2N(i,:), ...
            (X_k(j-1,:) - X_k(i,:))/(2*N_k) + d3_k(j-1,:)/(3*Nk2));
        j_col_tau = 3*(N_k+1) + 3*(i-1) + 1;
        for j = i:N_k
            B_k(i_row, j_col_tau)   = d2N(i,1) / N_k;
            B_k(i_row, j_col_tau+1) = d2N(i,2) / N_k;
            B_k(i_row, j_col_tau+2) = d2N(i,3) / N_k;
            j_col_tau = j_col_tau + 3;
        end
    end

    for i = 1:N_k
        i_row = 3*(i-1) + 3;
        j = i;  j_col = 3*(j-1);
        B_k(i_row, j_col+1:j_col+3) = cross(d3N(i,:), ...
            (X_k(j,:) - X_k(i,:))/(2*N_k) + d3_k(j,:)/(6*Nk2));
        if i < N_k
            K_mat = [0, -d3N(i,3), d3N(i,2); d3N(i,3), 0, -d3N(i,1); ...
                     -d3N(i,2), d3N(i,1), 0];
            d_mat   = inter(i:end,:) - X_k(i,:)/N_k;
            crosses = K_mat * d_mat';
            B_k(i_row, 3*i+1:3*N_k) = crosses(:)';
        end
        j = N_k + 1;  j_col = 3*(j-1);
        B_k(i_row, j_col+1:j_col+3) = cross(d3N(i,:), ...
            (X_k(j-1,:) - X_k(i,:))/(2*N_k) + d3_k(j-1,:)/(3*Nk2));
        j_col_tau = 3*(N_k+1) + 3*(i-1) + 1;
        for j = i:N_k
            B_k(i_row, j_col_tau)   = d3N(i,1) / N_k;
            B_k(i_row, j_col_tau+1) = d3N(i,2) / N_k;
            B_k(i_row, j_col_tau+2) = d3N(i,3) / N_k;
            j_col_tau = j_col_tau + 3;
        end
    end

    R_f_trans = kron(eye(N_k+1), rot_k);
    R_f_twist = kron(eye(N_k), rot_k);
    R_f_block = blkdiag(R_f_trans, R_f_twist);
    B_k_lab   = ell_k_B^3 * (B_k * R_f_block);

    f_cols_start = col_offset_trans + 1;
    f_cols_end   = col_offset_trans + 3*(N_k+1);
    t_cols_start = col_offset_twist + 1;
    t_cols_end   = col_offset_twist + 3*N_k;

    B_global(row_offset+1 : row_offset+3*N_k, f_cols_start:f_cols_end) = ...
        B_k_lab(:, 1:3*(N_k+1));
    B_global(row_offset+1 : row_offset+3*N_k, t_cols_start:t_cols_end) = ...
        B_k_lab(:, 3*(N_k+1)+1:end);

    row_offset       = row_offset + 3 * N_k;
    col_offset_trans = col_offset_trans + 3 * (N_k + 1);
    col_offset_twist = col_offset_twist + 3 * N_k;
end

%% Step 6: Assemble lin_sys = -B * AQ
n_trans  = 3 * N_total;
Q_trans  = Q_global(1:n_trans, :);

if ~isempty(A_twist_total)
    Q_twist  = Q_global(n_trans+1:end, :);
    AQ_twist = A_twist_total * Q_twist;

    if ~isempty(M_rotlet)
        Q_trans_corr = Q_trans - M_rotlet * AQ_twist;
    else
        Q_trans_corr = Q_trans;
    end
    AQ_trans = M_weighted \ Q_trans_corr;
    AQ       = [AQ_trans; AQ_twist];
else
    AQ_trans = M_weighted \ Q_trans;
    AQ = AQ_trans;
end

lin_sys = -B_global * AQ;

%% Step 6b: Optionally export hydrodynamic intermediates (for background flow)
if nargout >= 3
    hydro.B_global      = B_global;
    hydro.M_weighted    = M_weighted;
    hydro.A_twist_total = A_twist_total;
    hydro.M_rotlet      = M_rotlet;
    hydro.n_trans        = n_trans;
end

%% Step 7: Apply clamped boundary conditions
col_offset = 6;
row_offset = 6;

for k = 1:M
    N_k    = N_segs(k);
    rot_k  = params.rot{k};
    theta1 = geom.pili_angles{k}.theta(1);
    phi1   = geom.pili_angles{k}.phi(1);

    sT = sin(theta1);  cT = cos(theta1);
    sP = sin(phi1);    cP = cos(phi1);

    if abs(sT) < 1e-12
        warning('build_coupled_lin_sys:SingularTheta', ...
            'sin(theta_1) = %.2e for pili %d -- Euler angle singularity.', sT, k);
        sT = sign(sT + 1e-20) * 1e-12;
    end

    r1 = row_offset + 1;
    r2 = row_offset + 2;
    r3 = row_offset + 3;

    omega_cols = 4:6;
    theta1_col = col_offset + 1;
    phi1_col   = col_offset + N_k + 1;
    psi1_col   = col_offset + 2*N_k + 1;

    J_euler_basis = [-sP,         cP,        0;
                     -cT*cP/sT,  -cT*sP/sT, 1;
                      cP/sT,      sP/sT,     0];
    J_euler = J_euler_basis * rot_k;

    lin_sys(r1, :) = 0;
    lin_sys(r2, :) = 0;
    lin_sys(r3, :) = 0;

    lin_sys(r1, theta1_col) = 1;   lin_sys(r1, omega_cols) = -J_euler(1,:);
    lin_sys(r2, phi1_col)   = 1;   lin_sys(r2, omega_cols) = -J_euler(2,:);
    lin_sys(r3, psi1_col)   = 1;   lin_sys(r3, omega_cols) = -J_euler(3,:);

    col_offset = col_offset + 3 * N_k;
    row_offset = row_offset + 3 * N_k;
end

end


%% ==================== Local Helper Function ====================

function w = get_trapezoidal_weight(global_idx, N_body, N_segs, ell)
%GET_TRAPEZOIDAL_WEIGHT  Trapezoidal quadrature weight for a pili point
    local_idx  = global_idx - N_body;
    cumulative = 0;
    for k = 1:length(N_segs)
        n_pts_k = N_segs(k) + 1;
        if local_idx <= cumulative + n_pts_k
            j_in_pili = local_idx - cumulative;
            ds_k      = ell(k) / N_segs(k);
            if j_in_pili == 1 || j_in_pili == n_pts_k
                w = ds_k / 2;
            else
                w = ds_k;
            end
            return;
        end
        cumulative = cumulative + n_pts_k;
    end
    w = 0;
end
