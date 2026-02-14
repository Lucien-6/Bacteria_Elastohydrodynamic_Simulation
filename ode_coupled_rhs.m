function dZ = ode_coupled_rhs(t, Z, params)
%ODE_COUPLED_RHS  Core ODE right-hand side for coupled body-pili system
%
%  Description:
%    Computes dZ/dt for the coupled bacterial body + elastic pili system.
%    Evaluation steps:
%      1. Reconstruct geometry from state vector Z
%      2. Build global mobility (RSM/RSS) + torsional resistance + rotlet
%      3. Build Q matrix (state rates -> velocities)
%      4. Build B matrix (forces -> balance equations)
%      5. Compute elastic restoring moments for each pilus
%      6. Assemble and solve the coupled linear system  L * z_dot = R
%         (with rotlet coupling: Q_trans_corr = Q_trans - M_rotlet*AQ_twist)
%      7. Apply clamped boundary conditions
%      8. Pack state rates into dZ
%
%  Syntax:
%    dZ = ode_coupled_rhs(t, Z, params)
%
%  Inputs:
%    t      - Current dimensionless time
%    Z      - [n x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    dZ - [n x 1] Time derivative of state vector
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-14 -- RSS rotlet segment coupling + vectorised assembly
%           2026-02-14 -- BUG FIX: added rot_k in Q base-velocity mapping

%% Extract frequently used parameters
N_body = params.N_body;
M      = params.M_pili;
N_segs = params.N_segments;
sigma  = params.sigma;
Eg     = params.Eg;
g_hat  = params.g_hat(:);

%% Step 1: Reconstruct geometry
geom    = reconstruct_geometry(Z, params);
X_all   = geom.X_all;
X_cm    = geom.X_cm;
N_total = size(X_all, 1);

%% Step 2: Build global mobility / resistance matrices
[M_weighted, A_twist_cells, M_rotlet] = build_global_resistance(X_all, params, geom);

%% Step 3-4: Build B, Q matrices
n_state = 6 + 3 * sum(N_segs);
n_vel   = 3 * N_total + sum(N_segs);
n_force = 3 * N_total + 3 * sum(N_segs);

%% Build Q_global: [n_vel x n_state]
Q_global = zeros(n_vel, n_state);

% --- Body point velocities: v_i = V + Omega x (x_i - X_cm) ---
for i = 1:N_body
    ii  = 3*(i-1) + 1;
    r_i = X_all(i,:)' - X_cm;
    Q_global(ii:ii+2, 1:3) = eye(3);
    Q_global(ii:ii+2, 4:6) = -skew_matrix(r_i);
end

% --- Pili node velocities and twist rates ---
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

    % === Walker's Q_k (spatial) in rot_k basis ===
    Q_k_spatial = zeros(3*(N_k+1), 2*N_k + 3);

    % Base-velocity columns (1-3)
    Q_k_spatial(1:N_k+1, 1)                     = 1;
    Q_k_spatial(N_k+2 : 2*(N_k+1), 2)           = 1;
    Q_k_spatial(2*(N_k+1)+1 : 3*(N_k+1), 3)     = 1;

    % Theta and Phi rate blocks
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

    % Reorder to interleaved [x1,y1,z1, x2,y2,z2, ...]
    row_order = zeros(3*(N_k+1), 1);
    for j = 1:N_k+1
        row_order(3*(j-1)+1 : 3*j) = [j, N_k+1+j, 2*(N_k+1)+j];
    end
    Q_k_spatial = Q_k_spatial(row_order, :);

    % Extend with psi (twist)
    Q_k_ext = blkdiag(Q_k_spatial, eye(N_k));
    Q_k_ext(3*(N_k+1)+1 : 3*(N_k+1)+N_k, 3+N_k+1 : 3+2*N_k) = diag(cos(theta_k));

    % === Transform to lab frame ===
    rot_k_inv    = rot_k';
    R_block_trans = kron(eye(N_k+1), rot_k_inv);
    R_block_full  = blkdiag(R_block_trans, eye(N_k));
    Q_k_lab       = R_block_full * Q_k_ext;

    % === Map base velocity to body DOFs (V, Omega) ===
    %  J_base gives the root velocity in LAB frame:
    %    v_root_lab = V + Omega x r_root = [I, -skew(r)] * [V; Omega]
    %
    %  Q_k_lab(:,1:3) was built from Q_k_spatial (pili basis) then
    %  transformed by rot_k' per node.  Its columns 1:3 expect the base
    %  velocity in the PILI BASIS, not the lab frame.
    %
    %  Correct mapping:  v_base_basis = rot_k * v_base_lab
    %    => Q_k_body = Q_k_lab(:,1:3) * rot_k * J_base
    %
    %  BUG FIX 2026-02-14: Previously rot_k was missing, which is benign
    %  when rot_k == I (no adaptive basis rotation yet) but introduces a
    %  spurious O(V_sink) velocity error once rot_k != I, causing the body
    %  to rotate away from its equilibrium orientation.
    J_base = zeros(3, 6);
    J_base(1:3, 1:3) = eye(3);
    J_base(1:3, 4:6) = -skew_matrix(r_root);
    Q_k_body = Q_k_lab(:, 1:3) * rot_k * J_base;

    % Place into Q_global
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

% === Global force balance (rows 1-3) ===
for i = 1:N_body
    ii = 3*(i-1) + 1;
    B_global(1, ii)   = 1;
    B_global(2, ii+1) = 1;
    B_global(3, ii+2) = 1;
end

% Pili points: trapezoidal quadrature
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

% === Global torque balance (rows 4-6) ===
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

% === Pili moment balance (rows 7 onward) ===
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

    % Interpolate directors at nodes
    segBases = zeros(3, 3, N_k);
    for j = 1:N_k
        segBases(:,:,j) = [d1_k(j,:)', d2_k(j,:)', d3_k(j,:)'];
    end

    meanRots = zeros(3, 3, N_k+1);
    meanRots(:,:,1) = 1.5*segBases(:,:,1) - 0.5*segBases(:,:,2);
    if N_k >= 2
        meanRots(:,:,2:end-1) = 0.5*(segBases(:,:,1:end-1) + segBases(:,:,2:end));
    end
    meanRots(:,:,end) = 1.5*segBases(:,:,end) - 0.5*segBases(:,:,max(1,end-1));

    for j = 1:N_k+1
        [U_svd, ~, V_svd] = svd(meanRots(:,:,j));
        meanRots(:,:,j) = U_svd * V_svd';
    end

    d1N = reshape(meanRots(:,1,:), [3, N_k+1])';
    d2N = reshape(meanRots(:,2,:), [3, N_k+1])';
    d3N = reshape(meanRots(:,3,:), [3, N_k+1])';

    % Build B_k (Walker convention, pili basis)
    B_k   = zeros(3*N_k, 3*(N_k+1) + 3*N_k);
    inter = (d3_k(2:end,:) + 2*d3_k(1:end-1,:)) / (6*Nk2) + ...
            (X_k(2:end-1,:) + X_k(1:end-2,:)) / (2*N_k);

    % d1 moment balance
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

    % d2 moment balance
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

    % d3 moment balance (twist)
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

    % Transform B_k from pili basis to lab frame
    R_f_trans = kron(eye(N_k+1), rot_k);
    R_f_twist = kron(eye(N_k), rot_k);
    R_f_block = blkdiag(R_f_trans, R_f_twist);
    B_k_lab   = ell_k_B^3 * (B_k * R_f_block);

    % Place into B_global
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

%% Step 5: Compute elastic restoring moments (RHS)
R_global = zeros(n_state, 1);
R_global(1:3) = -Eg * g_hat;

row_offset = 6;
for k = 1:M
    N_k = N_segs(k);
    [kappa1, kappa2, kappa3] = compute_pili_elastic_moments( ...
        geom.d1_basis{k}, geom.d2_basis{k}, geom.d3_basis{k}, ...
        N_k, sigma, params.kb0{k});
    for j = 1:N_k
        R_global(row_offset + 3*(j-1) + 1) = kappa1(j);
        R_global(row_offset + 3*(j-1) + 2) = kappa2(j);
        R_global(row_offset + 3*(j-1) + 3) = kappa3(j);
    end
    row_offset = row_offset + 3 * N_k;
end

%% Step 6: Assemble and solve  L * z_dot = R
n_trans  = 3 * N_total;
Q_trans  = Q_global(1:n_trans, :);

if ~isempty(A_twist_total)
    Q_twist  = Q_global(n_trans+1:end, :);
    AQ_twist = A_twist_total * Q_twist;

    % Rotlet coupling: torque density generates translational velocity.
    % Correct Q_trans before solving the Stokeslet system.
    %   v_trans = M_stokeslet * f + M_rotlet * tau
    %   => f = M_stokeslet \ (Q_trans - M_rotlet * AQ_twist) * z_dot
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

    % Guard against singularity (defensive -- should be prevented by
    % adaptive basis rotation)
    if abs(sT) < 1e-12
        warning('ode_coupled_rhs:SingularTheta', ...
            'sin(theta_1) = %.2e for pili %d -- Euler angle singularity. Check adaptive_basis_check.', sT, k);
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

    R_global(r1) = 0;
    R_global(r2) = 0;
    R_global(r3) = 0;

    col_offset = col_offset + 3 * N_k;
    row_offset = row_offset + 3 * N_k;
end

%% Solve
z_dot = lin_sys \ R_global;

%% Step 8: Pack state rates
V     = z_dot(1:3);
Omega = z_dot(4:6);

dqdt       = quat_kinematics(geom.q, Omega);
q_curr     = geom.q;
dqdt       = dqdt - dot(dqdt, q_curr) * q_curr;

pili_rates = cell(params.M_pili, 1);
col_offset = 6;
for k = 1:params.M_pili
    N_k = N_segs(k);
    pili_rates{k}.dtheta = z_dot(col_offset + 1         : col_offset + N_k);
    pili_rates{k}.dphi   = z_dot(col_offset + N_k + 1   : col_offset + 2*N_k);
    pili_rates{k}.dpsi   = z_dot(col_offset + 2*N_k + 1 : col_offset + 3*N_k);
    col_offset = col_offset + 3 * N_k;
end

dZ = pack_state_rates(V, dqdt, pili_rates, N_segs);

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
