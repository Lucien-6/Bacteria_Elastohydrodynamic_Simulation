function [B_global, twist_seg_lengths, B_blocks] = build_force_balance(geom, params)
%BUILD_FORCE_BALANCE  Force / moment-balance operator of the coupled system
%
%  Description:
%    Assembles B_global, the operator that maps the hydrodynamic force
%    unknowns (point forces on the body, force densities on the pilus
%    nodes, torque densities on the pilus segments) to the rows of the
%    coupled system: total force (rows 1:3), total torque about the centre
%    of mass (rows 4:6), and the integrated moment balance of every pilus
%    segment in the local director frame (three rows per segment).
%
%    B_global depends on the configuration only (node positions and
%    directors), not on the hydrodynamic mobility.  The assembly is the one
%    formerly inlined in build_coupled_lin_sys (Steps 3-4), moved here
%    unchanged so that it can also be evaluated without the mobility
%    assembly (compute_effective_coupled_brownian, held-load drift).
%
%  Syntax:
%    [B_global, twist_seg_lengths] = build_force_balance(geom, params)
%
%  Inputs:
%    geom   - Geometry struct from reconstruct_geometry
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    B_global          - [n_state x n_force] operator, n_state =
%                        6 + 3*sum(N_segments), n_force = 3*N_total +
%                        3*sum(N_segments)
%    twist_seg_lengths - [1 x sum(N_segments)] current segment lengths
%                        (weights of the net axial couple)
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] Cortez (2018), J. Comput. Phys. 375, 783-796.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-12

N_body = params.N_body;
M      = params.M_pili;
N_segs = params.N_segments;

X_all   = geom.X_all;
X_cm    = geom.X_cm;
N_total = size(X_all, 1);

n_state = 6 + 3 * sum(N_segs);
n_force = 3 * N_total + 3 * sum(N_segs);

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
%  ds/2 at the two end nodes and ds at interior nodes of each pilus).
r_all = X_all - X_cm';                                % [N_total x 3]
w_all = params.quad_weights(:)';                      % 1 x N_total
rax   = r_all(:,1)';  ray = r_all(:,2)';  raz = r_all(:,3)';

B_global(4, 2:3:3*N_total) = -w_all .* raz;
B_global(4, 3:3:3*N_total) =  w_all .* ray;
B_global(5, 1:3:3*N_total) =  w_all .* raz;
B_global(5, 3:3:3*N_total) = -w_all .* rax;
B_global(6, 1:3:3*N_total) = -w_all .* ray;
B_global(6, 2:3:3*N_total) =  w_all .* rax;

%  Net axial couple in the system angular-momentum balance.  Cortez (2018)
%  defines the rotlet source as a torque density on a line segment, so the
%  physical net couple is the segment line integral sum_s L_s*(D_twist*m)_s.
%  The twist block of the force vector stores these 3D torque densities.
twist_seg_lengths = zeros(1, sum(N_segs));
seg_count = 0;
for k = 1:M
    N_k = N_segs(k);
    X_k = geom.X_pili{k};
    seg_lengths_k = sqrt(sum(diff(X_k, 1, 1).^2, 2));
    twist_seg_lengths(seg_count + (1:N_k)) = seg_lengths_k(:)';
    seg_count = seg_count + N_k;
end

B_global(4, 3*N_total+1 : 3 : n_force) = twist_seg_lengths;
B_global(5, 3*N_total+2 : 3 : n_force) = twist_seg_lengths;
B_global(6, 3*N_total+3 : 3 : n_force) = twist_seg_lengths;

row_offset       = 6;
col_offset_trans = 3 * N_body;
col_offset_twist = 3 * N_total;

want_blocks = (nargout >= 3);
if want_blocks
    B_blocks = repmat(struct('row_start', 0, 'row_end', 0, ...
                             'f_start', 0, 'f_end', 0, 'B_force', [], ...
                             't_start', 0, 't_end', 0, 'B_tau', []), 1, M);
else
    B_blocks = [];
end

has_cached_dN = isfield(geom, 'd1N_basis') && numel(geom.d1N_basis) >= M;

for k = 1:M
    N_k    = N_segs(k);
    Nk2    = N_k^2;
    rot_k  = params.rot{k};
    ell_k_B = params.ell(k);

    X_k_lab = geom.X_pili{k};
    X_k     = (X_k_lab * rot_k') / ell_k_B;

    d1_k = geom.d1_basis{k};
    d2_k = geom.d2_basis{k};
    d3_k = geom.d3_basis{k};

    if has_cached_dN && ~isempty(geom.d1N_basis{k})
        d1N = geom.d1N_basis{k};
        d2N = geom.d2N_basis{k};
        d3N = geom.d3N_basis{k};
    else
        [d1N, d2N, d3N] = compute_node_directors(d1_k, d2_k, d3_k, N_k);
    end

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

    %  Right-multiply every 3-column block by rot_k via batched pagemtimes
    %  (replaces dense kron(eye, rot_k) matrix multiplication).
    n_blk_k = 2 * N_k + 1;
    B_k_rot = reshape(pagemtimes(reshape(B_k, 3*N_k, 3, n_blk_k), rot_k), ...
                      3*N_k, 3*n_blk_k);

    f_cols_start = col_offset_trans + 1;
    f_cols_end   = col_offset_trans + 3*(N_k+1);
    t_cols_start = col_offset_twist + 1;
    t_cols_end   = col_offset_twist + 3*N_k;

    n_force_cols_k = 3 * (N_k + 1);
    B_k_lab_force  = ell_k_B^3 * B_k_rot(:, 1:n_force_cols_k);
    B_k_lab_tau    = ell_k_B^2 * B_k_rot(:, n_force_cols_k+1:end);

    B_global(row_offset+1 : row_offset+3*N_k, f_cols_start:f_cols_end) = ...
        B_k_lab_force;
    B_global(row_offset+1 : row_offset+3*N_k, t_cols_start:t_cols_end) = ...
        B_k_lab_tau;

    if want_blocks
        B_blocks(k).row_start = row_offset + 1;
        B_blocks(k).row_end   = row_offset + 3 * N_k;
        B_blocks(k).f_start   = f_cols_start;
        B_blocks(k).f_end     = f_cols_end;
        B_blocks(k).B_force   = B_k_lab_force;
        B_blocks(k).t_start   = t_cols_start - 3 * N_total;
        B_blocks(k).t_end     = t_cols_end   - 3 * N_total;
        B_blocks(k).B_tau     = B_k_lab_tau;
    end

    row_offset       = row_offset + 3 * N_k;
    col_offset_trans = col_offset_trans + 3 * (N_k + 1);
    col_offset_twist = col_offset_twist + 3 * N_k;
end

end
