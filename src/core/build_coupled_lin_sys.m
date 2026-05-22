function [lin_sys, geom, hydro] = build_coupled_lin_sys(Z, params)
%BUILD_COUPLED_LIN_SYS  Assemble the coupled body-pili linear system
%
%  Description:
%    Builds the full coupled linear system for the body-pili
%    elastohydrodynamic problem:
%      lin_sys * z_dot = R
%    where lin_sys encodes the hydrodynamic resistance (mobility
%    inversion + force/moment balance) and root boundary conditions
%    (clamped or free-torsion per pili).
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
%                .B_global, .M_weighted, .dM_weighted,
%                .A_twist_total, .M_rotlet, .n_trans
%              .dM_weighted is a reusable LU decomposition object for
%              M_weighted \ rhs solves (lets apply_Ainv_u skip a second
%              full LU factorisation in the background-flow branch of
%              ode_coupled_rhs).
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
[M_weighted, A_twist_cells, M_rotlet] = build_global_resistance(X_all, params, geom);

%% Step 3-4: Build B, Q matrices
%  NOTE: Step 5 (elastic restoring moments RHS) is intentionally omitted
%  from this function — it remains in ode_coupled_rhs.m where it is
%  assembled after the linear system solve.  Hence the step numbering
%  jumps from 4 to 6 here.
n_state = 6 + 3 * sum(N_segs);
n_vel   = 3 * N_total + sum(N_segs);
n_force = 3 * N_total + 3 * sum(N_segs);

%% Build Q_global: [n_vel x n_state]
Q_global = zeros(n_vel, n_state);

%  Vectorised body block of Q_global.  For i = 1..N_body:
%    rows (3i-2, 3i-1, 3i), cols 1..3 are identity 3x3
%    rows (3i-2, 3i-1, 3i), cols 4..6 are -skew(r_i)
%  where r_i = X_all(i,:) - X_cm.  -skew(r) = [0 z -y; -z 0 x; y -x 0].
r_body = X_all(1:N_body, :) - X_cm';              % [N_body x 3]
rbx    = r_body(:,1);  rby = r_body(:,2);  rbz = r_body(:,3);

idx_bx = (1 : 3 : 3*N_body)';
idx_by = (2 : 3 : 3*N_body)';
idx_bz = (3 : 3 : 3*N_body)';

Q_global(idx_bx, 1) = 1;  Q_global(idx_by, 2) = 1;  Q_global(idx_bz, 3) = 1;
Q_global(idx_bx, 5) =  rbz;   Q_global(idx_bx, 6) = -rby;
Q_global(idx_by, 4) = -rbz;   Q_global(idx_by, 6) =  rbx;
Q_global(idx_bz, 4) =  rby;   Q_global(idx_bz, 5) = -rbx;

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

    %  Vectorised lower-triangular fill of Q12/Q13/Q22/Q23/Q32.  Each block
    %  has entry (i,j) equal to the per-j trig product for i > j, else 0.
    cT = cos(theta_k(:))';     sT = sin(theta_k(:))';     % 1 x N_k
    cP = cos(phi_k(:))';       sP = sin(phi_k(:))';
    maskL = tril(true(N_k+1, N_k), -1);                    % (N_k+1) x N_k
    Q12 = maskL .* ( cT .* cP);
    Q13 = maskL .* (-sT .* sP);
    Q22 = maskL .* ( cT .* sP);
    Q23 = maskL .* ( sT .* cP);
    Q32 = maskL .* (-sT);

    Q_k_spatial(:, 4:end) = [Q12, Q13; Q22, Q23; Q32, zeros(N_k+1,N_k)] ...
                            * (ell_k / N_k);

    %  Vectorised row-permutation: interleave x/y/z rows.  For j=1..N_k+1,
    %  row_order(3(j-1)+1 : 3j) = [j, (N_k+1)+j, 2(N_k+1)+j].
    j_vec     = (1 : N_k+1)';
    row_order = reshape([j_vec, (N_k+1)+j_vec, 2*(N_k+1)+j_vec]', [], 1);
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

%% Build torsional resistance block (single blkdiag call; O(M) vs O(M^2))
if M > 0
    A_twist_total = blkdiag(A_twist_cells{:});
else
    A_twist_total = [];
end

%% Build B_global: [n_state x n_force]
B_global = zeros(n_state, n_force);

%  Vectorised identity stripe for body delta-force rows
B_global(1, 1:3:3*N_body) = 1;
B_global(2, 2:3:3*N_body) = 1;
B_global(3, 3:3:3*N_body) = 1;

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

%  Vectorised angular-moment rows.  Uses the precomputed quadrature
%  weights from init_coupled_params (w=1 for body points, trapezoidal
%  ds/2 or ds for pili node points -- identical to get_trapezoidal_weight).
r_all = X_all - X_cm';                                % [N_total x 3]
w_all = params.quad_weights(:)';                      % 1 x N_total
rax   = r_all(:,1)';  ray = r_all(:,2)';  raz = r_all(:,3)';

B_global(4, 2:3:3*N_total) = -w_all .* raz;
B_global(4, 3:3:3*N_total) =  w_all .* ray;
B_global(5, 1:3:3*N_total) =  w_all .* raz;
B_global(5, 3:3:3*N_total) = -w_all .* rax;
B_global(6, 1:3:3*N_total) = -w_all .* ray;
B_global(6, 2:3:3*N_total) =  w_all .* rax;

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

    %  -------- Vectorised B_k assembly (4 sub-blocks) --------
    %
    %  Block layout (row index a = 1/2/3 ↔ d1/d2/d3):
    %    Off-diag : i=1..N_k-1, j=i+1..N_k    cols 3*(j-1)+(1:3)
    %    Diagonal : i=1..N_k,   j=i            cols 3*(i-1)+(1:3)
    %    Tip      : i=1..N_k,   j=N_k+1        cols 3*N_k+(1:3)
    %    Tau      : i=1..N_k,   j=i..N_k       cols 3(N_k+1)+3(j-1)+(1:3)
    %
    %  Ordering: we write Off-diag first, then Diagonal, then Tip, then Tau.
    %  The masked Off-diag writes may put zeros into the diagonal 3x3 cells
    %  (when jm = i-1); those zeros are subsequently overwritten by the
    %  Diagonal block with the correct cross-product values.  The Tip and
    %  Tau blocks occupy disjoint column ranges, so their order is free.

    row_d1 = 1:3:3*N_k;
    row_d2 = 2:3:3*N_k;
    row_d3 = 3:3:3*N_k;

    %  NOTE: compute_node_directors returns d1N/d2N/d3N of size [(N_k+1) x 3].
    %  The original loops use only rows i = 1..N_k.  We slice here so that
    %  broadcasting dimensions match the arg_off matrix [N_k x (N_k-1)].
    d1N_i = d1N(1:N_k, :);    % [N_k x 3]  (rows used by row i = 1..N_k)
    d2N_i = d2N(1:N_k, :);
    d3N_i = d3N(1:N_k, :);

    %% Block 1: Off-diagonal (j > i, j <= N_k)
    if N_k > 1
        %  arg_off(i, jm, c) = inter(jm, c) - X_k(i, c)/N_k  for jm = j-1 = 1..N_k-1
        Xi_N      = X_k(1:N_k, :) / N_k;                % [N_k x 3]
        arg_off_1 = inter(:,1)' - Xi_N(:,1);            % [N_k x (N_k-1)]
        arg_off_2 = inter(:,2)' - Xi_N(:,2);
        arg_off_3 = inter(:,3)' - Xi_N(:,3);

        %  Component-wise cross product  (c-th component of  dN x arg_off)
        %    cross_c1 = dN_2 .* arg3 - dN_3 .* arg2
        %    cross_c2 = dN_3 .* arg1 - dN_1 .* arg3
        %    cross_c3 = dN_1 .* arg2 - dN_2 .* arg1
        %  Mathematically identical (and bit-exact in IEEE-754 forward
        %  accumulation) to the original  K_mat * d_mat'  formulation.
        d1_1 = d1N_i(:,1);  d1_2 = d1N_i(:,2);  d1_3 = d1N_i(:,3);
        d2_1 = d2N_i(:,1);  d2_2 = d2N_i(:,2);  d2_3 = d2N_i(:,3);
        d3_1 = d3N_i(:,1);  d3_2 = d3N_i(:,2);  d3_3 = d3N_i(:,3);

        co_d1_1 = d1_2 .* arg_off_3 - d1_3 .* arg_off_2;
        co_d1_2 = d1_3 .* arg_off_1 - d1_1 .* arg_off_3;
        co_d1_3 = d1_1 .* arg_off_2 - d1_2 .* arg_off_1;

        co_d2_1 = d2_2 .* arg_off_3 - d2_3 .* arg_off_2;
        co_d2_2 = d2_3 .* arg_off_1 - d2_1 .* arg_off_3;
        co_d2_3 = d2_1 .* arg_off_2 - d2_2 .* arg_off_1;

        co_d3_1 = d3_2 .* arg_off_3 - d3_3 .* arg_off_2;
        co_d3_2 = d3_3 .* arg_off_1 - d3_1 .* arg_off_3;
        co_d3_3 = d3_1 .* arg_off_2 - d3_2 .* arg_off_1;

        mask_UT = triu(true(N_k, N_k-1), 0);            % true when jm >= i

        co_d1_1 = co_d1_1 .* mask_UT;  co_d1_2 = co_d1_2 .* mask_UT;  co_d1_3 = co_d1_3 .* mask_UT;
        co_d2_1 = co_d2_1 .* mask_UT;  co_d2_2 = co_d2_2 .* mask_UT;  co_d2_3 = co_d2_3 .* mask_UT;
        co_d3_1 = co_d3_1 .* mask_UT;  co_d3_2 = co_d3_2 .* mask_UT;  co_d3_3 = co_d3_3 .* mask_UT;

        cols_co_1 = 3*(1:N_k-1) + 1;
        cols_co_2 = 3*(1:N_k-1) + 2;
        cols_co_3 = 3*(1:N_k-1) + 3;

        B_k(row_d1, cols_co_1) = co_d1_1;
        B_k(row_d1, cols_co_2) = co_d1_2;
        B_k(row_d1, cols_co_3) = co_d1_3;
        B_k(row_d2, cols_co_1) = co_d2_1;
        B_k(row_d2, cols_co_2) = co_d2_2;
        B_k(row_d2, cols_co_3) = co_d2_3;
        B_k(row_d3, cols_co_1) = co_d3_1;
        B_k(row_d3, cols_co_2) = co_d3_2;
        B_k(row_d3, cols_co_3) = co_d3_3;
    end

    %% Block 2: Diagonal (j = i)
    %  arg = (X_k(i,:)-X_k(i,:))/(2*N_k) + d3_k(i,:)/(6*Nk2) = d3_k(i,:)/(6*Nk2)
    %  The zero subtraction term is bit-exact (x + 0 = x), so dropping it is
    %  IEEE-754 equivalent to the original.
    arg_diag = d3_k(1:N_k, :) / (6*Nk2);                % [N_k x 3]
    cd_d1    = cross(d1N_i, arg_diag, 2);
    cd_d2    = cross(d2N_i, arg_diag, 2);
    cd_d3    = cross(d3N_i, arg_diag, 2);

    cc_diag  = 3*(0:N_k-1)' + (1:3);                    % [N_k x 3] column indices
    rr_d1    = repmat((3*(0:N_k-1)' + 1), 1, 3);
    lin_d1   = sub2ind(size(B_k), rr_d1,     cc_diag);
    lin_d2   = sub2ind(size(B_k), rr_d1 + 1, cc_diag);
    lin_d3   = sub2ind(size(B_k), rr_d1 + 2, cc_diag);
    B_k(lin_d1) = cd_d1;
    B_k(lin_d2) = cd_d2;
    B_k(lin_d3) = cd_d3;

    %% Block 3: Tip (j = N_k + 1)
    arg_tip = (X_k(N_k, :) - X_k(1:N_k, :)) / (2*N_k) + d3_k(N_k, :) / (3*Nk2);
    ct_d1   = cross(d1N_i, arg_tip, 2);
    ct_d2   = cross(d2N_i, arg_tip, 2);
    ct_d3   = cross(d3N_i, arg_tip, 2);

    tip_cols = 3*N_k + (1:3);
    B_k(row_d1, tip_cols) = ct_d1;
    B_k(row_d2, tip_cols) = ct_d2;
    B_k(row_d3, tip_cols) = ct_d3;

    %% Block 4: Tau (upper triangular, j >= i)
    tau_offset = 3*(N_k+1);
    maskU_tau  = triu(true(N_k, N_k));
    for c = 1:3
        cols_tau = tau_offset + c : 3 : tau_offset + 3*N_k;
        B_k(row_d1, cols_tau) = maskU_tau .* (d1N_i(:, c) / N_k);
        B_k(row_d2, cols_tau) = maskU_tau .* (d2N_i(:, c) / N_k);
        B_k(row_d3, cols_tau) = maskU_tau .* (d3N_i(:, c) / N_k);
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

% ---- Direction 1: factor M_weighted ONCE, reuse across solves.            ----
% M_weighted mixes RSM (body point-Stokeslets) and RSS (pili segment
% integrals); the resulting block matrix is non-symmetric by construction,
% so LU is the right factorisation.
%
% Build a reusable decomposition so the SAME LU serves:
%   (a) this function's  AQ_trans = M \ Q_trans  solve;
%   (b) ode_coupled_rhs.apply_Ainv_u's  M \ u_t  solve (background-flow
%       branch).  Case (b) would otherwise trigger a second full O(N^3)
%       LU on every RHS call; reuse saves it entirely.
%
% CheckCondition=false suppresses MATLAB's automatic rcond estimate on
% each solve (O(N^2) per invocation); overall condition health is
% validated externally by the regression tests.  Measurement also shows
% the explicit LU dispatch is marginally faster than the generic
% backslash path even in the single-solve (bg=off) case, because it
% skips backslash's per-call symmetry / triangular heuristics.
dM_weighted = decomposition(M_weighted, 'lu', 'CheckCondition', false);

if ~isempty(A_twist_total)
    Q_twist  = Q_global(n_trans+1:end, :);
    AQ_twist = A_twist_total * Q_twist;

    if ~isempty(M_rotlet)
        Q_trans_corr = Q_trans - M_rotlet * AQ_twist;
    else
        Q_trans_corr = Q_trans;
    end
    AQ_trans = dM_weighted \ Q_trans_corr;
    AQ       = [AQ_trans; AQ_twist];
else
    AQ_trans = dM_weighted \ Q_trans;
    AQ = AQ_trans;
end

lin_sys = -B_global * AQ;

%% Step 6b: Optionally export hydrodynamic intermediates (for background flow)
if nargout >= 3
    hydro.B_global      = B_global;
    hydro.M_weighted    = M_weighted;
    hydro.dM_weighted   = dM_weighted;
    hydro.A_twist_total = A_twist_total;
    hydro.M_rotlet      = M_rotlet;
    hydro.n_trans        = n_trans;
end

%% Step 7: Apply root boundary conditions (clamped or free-torsion)
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
        logf_h = get_params_logf(params);
        logf_h('[ODE] WARNING: sin(theta_1) = %.2e for pili %d -- Euler angle singularity.\n', sT, k);
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

    rc_k = 'clamped';
    if isfield(params, 'root_coupling') && numel(params.root_coupling) >= k
        rc_k = params.root_coupling{k};
    end

    if strcmp(rc_k, 'clamped')
        % Full clamped BC: replace 3 rows (theta_1, phi_1, psi_1 all locked)
        lin_sys(r1, :) = 0;
        lin_sys(r2, :) = 0;
        lin_sys(r3, :) = 0;

        lin_sys(r1, theta1_col) = 1;   lin_sys(r1, omega_cols) = -J_euler(1,:);
        lin_sys(r2, phi1_col)   = 1;   lin_sys(r2, omega_cols) = -J_euler(2,:);
        lin_sys(r3, psi1_col)   = 1;   lin_sys(r3, omega_cols) = -J_euler(3,:);

    elseif strcmp(rc_k, 'free_torsion')
        % Free-torsion BC: replace only 2 rows (theta_1, phi_1 locked; psi_1 free)
        lin_sys(r1, :) = 0;
        lin_sys(r2, :) = 0;
        % r3 row is kept: elastic moment balance equation drives psi_1 dynamics

        lin_sys(r1, theta1_col) = 1;   lin_sys(r1, omega_cols) = -J_euler(1,:);
        lin_sys(r2, phi1_col)   = 1;   lin_sys(r2, omega_cols) = -J_euler(2,:);

    else
        error('build_coupled_lin_sys:UnknownRootCoupling', ...
            ['Pili %d: root_coupling = ''%s''. ', ...
             'Supported: clamped, free_torsion.'], k, rc_k);
    end

    % Diagnostic export of J_euler(3,:), the omega -> dot_psi kinematic
    % row.  Production code (v7.3.0+) does NOT add any motor reaction
    % term to R_global(4:6) — the system moment balance must remain
    % free of internal motor torques, which cancel by Newton's third
    % law (see ode_coupled_rhs.m §"NOTE (v7.3.0)").  This field is
    % retained purely for offline diagnostics that compare the
    % omega -> dot_psi projection against d3_lab (see CHANGELOG v7.3.0).
    geom.J_euler3{k} = J_euler(3,:);

    col_offset = col_offset + 3 * N_k;
    row_offset = row_offset + 3 * N_k;
end

end


%  Note: the former local helper get_trapezoidal_weight is superseded by
%  params.quad_weights, precomputed once in init_coupled_params.  This
%  eliminates per-call iteration over all pili points at every ODE RHS.


function logf_h = get_params_logf(params)
%GET_PARAMS_LOGF  Retrieve logf handle from params, or fallback to fprintf
    if isstruct(params) && isfield(params, 'logf') && ~isempty(params.logf)
        logf_h = params.logf;
    else
        logf_h = @(varargin) fprintf(varargin{:});
    end
end
