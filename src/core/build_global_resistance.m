function [M_weighted, M_twist, M_rotlet, D_twist] = build_global_resistance(X_all, params, geom)
%BUILD_GLOBAL_RESISTANCE  Global hydrodynamic mobility and twist mobility
%
%  Description:
%    Builds the translational mobility matrix for ALL points (body + pili)
%    and the full non-local axial twist mobility for pili segments.
%
%    Mobility method:
%      Body columns :  RSM point forces (Cortez 2005)
%      Pili columns :  RSS exact segment integrals (Cortez 2018)
%      Rotlet       :  RSS rotlet segment integrals (Cortez 2018)
%
%    Translational:  Returns M_weighted (NOT its inverse).  The caller
%    computes M_weighted \ Q directly, avoiding O(N^3) full inversion.
%
%    Torsional:  M_twist maps scalar axial torque densities on all pili
%                segments to scalar axial fluid angular velocities at all
%                segment midpoints.  It replaces the former diagonal local
%                RFT resistance with a full non-local segment-integral
%                mobility.
%
%    Rotlet segment mobility (Cortez 2018):
%      M_rotlet maps 3D torque densities (per segment) to translational
%      velocities (per eval point) using Cortez's full-curl rotlet.  Callers
%      use 0.5*M_rotlet for the physical standard rotlet so that the coupled
%      [u; omega] mobility is Lorentz-reciprocal.
%
%    Adaptive epsilon strategy:
%      Body points  :  eps_body = h_mean / 2.6667 (default; batch sweep optimum;
%                      see init_coupled_params.m)
%      Pili points  :  eps_pili = hydro_radius_factor * r_pili / L_0
%                      (params.epsilon; default factor sqrt(2), see init_coupled_params.m)
%      Pair rule    :  eps_ij = eps_j  (force-point-determined)
%
%    Matrix layout (block view; rows = evaluation points, columns = force
%    (source) points, block name = M_<eval><source>):
%      M_weighted = [ M_BB  M_BP ]   BB : eval body, source body  (RSM)
%                   [ M_PB  M_PP ]   BP : eval body, source pili  (RSS)
%                                    PB : eval pili, source body  (RSM)
%                                    PP : eval pili, source pili  (RSS)
%      Column method is determined by the force point (source j):
%        Body columns (j = body point) → RSM, eps = eps_body, w = 1
%        Pili columns (j = pili point) → RSS, eps = eps_pili, w = 1
%      (RSS integrals already incorporate the segment-length weight)
%
%  Syntax:
%    [M_weighted, M_twist, M_rotlet, D_twist] = build_global_resistance(X_all, params, geom)
%
%  Inputs:
%    X_all  - [N_total x 3] All point positions in lab frame
%    params - Parameter structure (from init_coupled_params)
%    geom   - Geometry structure (from reconstruct_geometry)
%
%  Outputs:
%    M_weighted    - [3*N_total x 3*N_total] Weighted mobility matrix
%    M_twist       - [N_seg_all x N_seg_all] Axial twist mobility
%    M_rotlet      - [3*N_total x 3*N_seg_all] Cortez rotlet segment mobility
%    D_twist       - [3*N_seg_all x N_seg_all] Axial torque embedding
%
%  References:
%    [1] Cortez et al. (2005), Phys. Fluids 17, 031504.
%    [2] Cortez (2018), J. Comput. Phys. 375, 783-796.
%    [3] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-11-04 — v8.7.2: comment and message-text corrections only; no behavioural change.
%  Updated: 2026-11-24 — exploit 3x3 tensor symmetry, reciprocal reuse, and in-place RSS scatter.

N_total = size(X_all, 1);
N_body  = params.N_body;
M       = params.M_pili;
N_segs  = params.N_segments;
dim     = 3 * N_total;

M_weighted = zeros(dim, dim);
e2_body    = params.eps_body^2;
eps_pili   = params.epsilon;
e2_pili    = eps_pili^2;

%% ---- Body columns (RSM vectorised, eps = eps_body) ----
%  All N_total evaluation points x N_body source points.
%  Force-point rule: source is body -> eps = eps_body everywhere.
%  Each 3x3 Stokeslet block G(i,j) is symmetric, so compute off-diagonals once.
rx_bb = X_all(:, 1) - X_all(1:N_body, 1)';
ry_bb = X_all(:, 2) - X_all(1:N_body, 2)';
rz_bb = X_all(:, 3) - X_all(1:N_body, 3)';
r2_bb = rx_bb.^2 + ry_bb.^2 + rz_bb.^2;

inv_R_bb = 1 ./ sqrt(r2_bb + e2_body);
H2_bb    = inv_R_bb .* inv_R_bb .* inv_R_bb;
H1_bb    = (r2_bb + 2 * e2_body) .* H2_bb;

row_x = 1:3:dim;           % N_total elements
row_y = 2:3:dim;
row_z = 3:3:dim;
col_x = 1:3:3*N_body;      % N_body elements
col_y = 2:3:3*N_body;
col_z = 3:3:3*N_body;

G_xy = H2_bb .* (rx_bb .* ry_bb);
G_xz = H2_bb .* (rx_bb .* rz_bb);
G_yz = H2_bb .* (ry_bb .* rz_bb);

M_weighted(row_x, col_x) = H1_bb + H2_bb .* (rx_bb .* rx_bb);
M_weighted(row_x, col_y) = G_xy;
M_weighted(row_x, col_z) = G_xz;
M_weighted(row_y, col_x) = G_xy;
M_weighted(row_y, col_y) = H1_bb + H2_bb .* (ry_bb .* ry_bb);
M_weighted(row_y, col_z) = G_yz;
M_weighted(row_z, col_x) = G_xz;
M_weighted(row_z, col_y) = G_yz;
M_weighted(row_z, col_z) = H1_bb + H2_bb .* (rz_bb .* rz_bb);

%% ---- Pre-collect all pili segment endpoints ----
N_seg_all = sum(N_segs);
if N_seg_all == 0
    M_twist  = [];
    M_rotlet = zeros(dim, 0);
    D_twist  = zeros(0, 0);
    return;
end

y0_all = zeros(N_seg_all, 3);     % Segment start (node k)
y1_all = zeros(N_seg_all, 3);     % Segment end   (node k+1)

s_count     = 0;
pili_offset = N_body;
for k = 1:M
    N_k = N_segs(k);
    idx_s = s_count + (1:N_k);
    y0_all(idx_s, :) = X_all(pili_offset + (1:N_k), :);
    y1_all(idx_s, :) = X_all(pili_offset + (2:N_k+1), :);
    s_count     = s_count + N_k;
    pili_offset = pili_offset + N_k + 1;
end

%% ---- Vectorised segment geometry  [1 x N_seg_all] ----
v_all  = y0_all - y1_all;                        % [N_seg x 3]
L2_all = sum(v_all.^2, 2);
L_all  = sqrt(L2_all);

% Transpose to row vectors for broadcasting [1 x N_seg]
L_r   = L_all';
L2_r  = L2_all';
v1_r  = v_all(:,1)';
v2_r  = v_all(:,2)';
v3_r  = v_all(:,3)';

%% ---- Displacements  [N_total x N_seg_all] ----
x01 = X_all(:,1) - y0_all(:,1)';
x02 = X_all(:,2) - y0_all(:,2)';
x03 = X_all(:,3) - y0_all(:,3)';

r0_sq = x01.^2 + x02.^2 + x03.^2;

x11 = x01 + v1_r;
x12 = x02 + v2_r;
x13 = x03 + v3_r;
r1_sq = x11.^2 + x12.^2 + x13.^2;

R0 = sqrt(r0_sq + e2_pili);
R1 = sqrt(r1_sq + e2_pili);
inv_R0 = 1 ./ R0;
inv_R1 = 1 ./ R1;

d0 = x01 .* v1_r + x02 .* v2_r + x03 .* v3_r;  % x0 . v

%% ---- Log arguments (strictly positive for eps > 0) ----
H0_v = max(L_r .* R0 + d0,        1e-300);
H1_v = max(L_r .* R1 + d0 + L2_r, 1e-300);

%% ---- T_{n,q} recursion: q = -1, -3  (Cortez 2018, Eq. 7) ----
inv_L  = 1 ./ L_r;
inv_L2 = 1 ./ L2_r;
d0_L2  = d0 .* inv_L2;                           % d0 / L^2

T0m1 = log(H1_v ./ H0_v) .* inv_L;
T0m3 = inv_R0 ./ H0_v - inv_R1 ./ H1_v;
T1m1 = (R1 - R0) .* inv_L2 - d0_L2 .* T0m1;
T1m3 = (inv_R0 - inv_R1) .* inv_L2 - d0_L2 .* T0m3;
neg_inv_R1_L2 = -inv_R1 .* inv_L2;
T2m3 = neg_inv_R1_L2 + T0m1 .* inv_L2 - d0_L2 .* T1m3;
T3m3 = neg_inv_R1_L2 + 2 * T1m1 .* inv_L2 - d0_L2 .* T2m3;

%% ---- T_{n,q} recursion: q = -5  (Cortez 2018, Sec. 2.4) ----
%  Needed for the regularised rotlet segment.
inv_R0_3 = inv_R0 .* inv_R0 .* inv_R0;
inv_R1_3 = inv_R1 .* inv_R1 .* inv_R1;

d_sq = max(L2_r .* (r0_sq + e2_pili) - d0.^2, e2_pili .* L2_r);
inv_d_sq_3 = (1 / 3) ./ d_sq;

d1_v   = d0 + L2_r;                              % x1 . v  =  x0.v + L^2
inv3L2 = inv_L2 / 3;

T0m5 = inv_d_sq_3 .* (d1_v .* inv_R1_3 - d0 .* inv_R0_3 + 2 * L2_r .* T0m3);
T1m5 = (inv_R0_3 - inv_R1_3) .* inv3L2 - d0_L2 .* T0m5;

%% ---- Stokeslet: M2 & M1 components (3x3 symmetric per pair) ----
%  Exploit symmetry M_21 == M_12, M_31 == M_13, M_32 == M_23 and share
%  the dyadic products (x0_a*x0_b, x0_a*v_b+x0_b*v_a, v_a*v_b) across M1 & M2.
c_iso = L_r .* (T1m1 + e2_pili * T1m3);
c1s   = L_r .* T1m3;
c2s   = L_r .* T2m3;
c3s   = L_r .* T3m3;

c_iso0 = L_r .* (T0m1 + e2_pili * T0m3) - c_iso;
c10    = L_r .* T0m3 - c1s;
c20    = c1s - c2s;          % L_r .* T1m3 - L_r .* T2m3
c30    = c2s - c3s;          % L_r .* T2m3 - L_r .* T3m3

xx11 = x01 .* x01;   xx12 = x01 .* x02;   xx13 = x01 .* x03;
xx22 = x02 .* x02;   xx23 = x02 .* x03;   xx33 = x03 .* x03;

xv11 = 2 * (x01 .* v1_r);
xv12 = x01 .* v2_r + x02 .* v1_r;
xv13 = x01 .* v3_r + x03 .* v1_r;
xv22 = 2 * (x02 .* v2_r);
xv23 = x02 .* v3_r + x03 .* v2_r;
xv33 = 2 * (x03 .* v3_r);

vv11 = v1_r .* v1_r;   vv12 = v1_r .* v2_r;   vv13 = v1_r .* v3_r;
vv22 = v2_r .* v2_r;   vv23 = v2_r .* v3_r;   vv33 = v3_r .* v3_r;

M2_11 = c_iso + c1s.*xx11 + c2s.*xv11 + c3s.*vv11;
M2_12 =         c1s.*xx12 + c2s.*xv12 + c3s.*vv12;
M2_13 =         c1s.*xx13 + c2s.*xv13 + c3s.*vv13;
M2_22 = c_iso + c1s.*xx22 + c2s.*xv22 + c3s.*vv22;
M2_23 =         c1s.*xx23 + c2s.*xv23 + c3s.*vv23;
M2_33 = c_iso + c1s.*xx33 + c2s.*xv33 + c3s.*vv33;

M1_11 = c_iso0 + c10.*xx11 + c20.*xv11 + c30.*vv11;
M1_12 =          c10.*xx12 + c20.*xv12 + c30.*vv12;
M1_13 =          c10.*xx13 + c20.*xv13 + c30.*vv13;
M1_22 = c_iso0 + c10.*xx22 + c20.*xv22 + c30.*vv22;
M1_23 =          c10.*xx23 + c20.*xv23 + c30.*vv23;
M1_33 = c_iso0 + c10.*xx33 + c20.*xv33 + c30.*vv33;

%% ---- Scatter Stokeslet M1/M2 into M_weighted (in-place per pilus) ----
%  Pili columns of M_weighted start at zero; write M1 into nodes 1..N_k
%  then accumulate M2 into nodes 2..N_k+1 (no concat scratch arrays).
seg_offset       = 0;
pili_node_offset = N_body;
for k = 1:M
    N_k     = N_segs(k);
    seg_idx = seg_offset + (1:N_k);

    col_base = 3 * pili_node_offset;
    cx1 = col_base + (1 : 3 : 3*N_k);
    cy1 = cx1 + 1;
    cz1 = cx1 + 2;
    cx2 = cx1 + 3;
    cy2 = cy1 + 3;
    cz2 = cz1 + 3;

    % Compute combined symmetric off-diagonal blocks once per pilus
    M_weighted(row_x, cx1) = M1_11(:, seg_idx);
    M_weighted(row_x, cx2) = M_weighted(row_x, cx2) + M2_11(:, seg_idx);

    M_weighted(row_y, cy1) = M1_22(:, seg_idx);
    M_weighted(row_y, cy2) = M_weighted(row_y, cy2) + M2_22(:, seg_idx);

    M_weighted(row_z, cz1) = M1_33(:, seg_idx);
    M_weighted(row_z, cz2) = M_weighted(row_z, cz2) + M2_33(:, seg_idx);

    M_weighted(row_x, cy1) = M1_12(:, seg_idx);
    M_weighted(row_x, cy2) = M_weighted(row_x, cy2) + M2_12(:, seg_idx);
    M_weighted(row_y, cx1) = M_weighted(row_x, cy1);
    M_weighted(row_y, cx2(end)) = M_weighted(row_x, cy2(end));

    M_weighted(row_x, cz1) = M1_13(:, seg_idx);
    M_weighted(row_x, cz2) = M_weighted(row_x, cz2) + M2_13(:, seg_idx);
    M_weighted(row_z, cx1) = M_weighted(row_x, cz1);
    M_weighted(row_z, cx2(end)) = M_weighted(row_x, cz2(end));

    M_weighted(row_y, cz1) = M1_23(:, seg_idx);
    M_weighted(row_y, cz2) = M_weighted(row_y, cz2) + M2_23(:, seg_idx);
    M_weighted(row_z, cy1) = M_weighted(row_y, cz1);
    M_weighted(row_z, cy2(end)) = M_weighted(row_y, cz2(end));

    pili_node_offset = pili_node_offset + N_k + 1;
    seg_offset       = seg_offset       + N_k;
end

%% ---- Rotlet segment mobility  (Cortez 2018, Sec. 2.5) ----
M_rotlet = zeros(dim, 3 * N_seg_all);

ca = 2 * T0m3 + 3 * e2_pili * T0m5;              % [N_total x N_seg_all]
cb = 2 * T1m3 + 3 * e2_pili * T1m5;

Rxy =  L_r .* (ca .* x03 + cb .* v3_r);          % R12
Rxz = -L_r .* (ca .* x02 + cb .* v2_r);          % R13
Ryz =  L_r .* (ca .* x01 + cb .* v1_r);          % R23

%  Vectorised scatter: each segment occupies its own 3 disjoint columns, so
%  the assignments are direct (no accumulation).  The strided writes below
%  address exactly the same (row, col) cells as the original loop.
M_rotlet(row_x, 2:3:end) =  Rxy;
M_rotlet(row_x, 3:3:end) =  Rxz;
M_rotlet(row_y, 1:3:end) = -Rxy;
M_rotlet(row_y, 3:3:end) =  Ryz;
M_rotlet(row_z, 1:3:end) = -Rxz;
M_rotlet(row_z, 2:3:end) = -Ryz;

%% ---- Axial torque embedding and full non-local twist mobility ----
D_twist = zeros(3 * N_seg_all, N_seg_all);
d3_all  = zeros(N_seg_all, 3);

seg_offset = 0;
for k = 1:M
    N_k = N_segs(k);
    idx = seg_offset + (1:N_k);
    d3_all(idx, :) = geom.d3{k};
    seg_offset = seg_offset + N_k;
end

if N_seg_all > 0
    rows_3xNs = reshape(1:3*N_seg_all, 3, N_seg_all);
    cols_3xNs = repmat(1:N_seg_all, 3, 1);
    lin_idx   = sub2ind([3*N_seg_all, N_seg_all], rows_3xNs, cols_3xNs);
    D_twist(lin_idx) = d3_all';

    M_twist = compute_nonlocal_twist_mobility(y0_all, y1_all, d3_all, eps_pili);
else
    M_twist = [];
end

end


function M_twist = compute_nonlocal_twist_mobility(y0_all, y1_all, d3_all, eps_pili)
%COMPUTE_NONLOCAL_TWIST_MOBILITY  Segment-integral axial rotation mobility
%
%  M_twist(i,j) maps a scalar axial torque density on source segment j
%  (torque vector = d3_j * m_j) to the scalar axial fluid angular velocity
%  at the midpoint of target segment i (omega_i = d3_i . Omega).
%
%  Cortez (2018) Eq. (3) defines the rotlet as the full curl of the
%  regularised Stokeslet.  The physical rod angular velocity is half the
%  fluid vorticity, and the reciprocal standard rotlet is half of Cortez's
%  full-curl rotlet; therefore the final prefactor is L/4.

mids  = 0.5 * (y0_all + y1_all);
v_all = y0_all - y1_all;
L_all = sqrt(sum(v_all.^2, 2));
L2_all = L_all .* L_all;

e2 = eps_pili^2;
e4 = e2^2;

% Source-segment data are row vectors for broadcasting over target rows.
L_r   = L_all';
L2_r  = L2_all';
v1_r  = v_all(:,1)';
v2_r  = v_all(:,2)';
v3_r  = v_all(:,3)';

% x(alpha) = target_mid - (y0 - alpha*v) = x0 + alpha*v.
x01 = mids(:,1) - y0_all(:,1)';
x02 = mids(:,2) - y0_all(:,2)';
x03 = mids(:,3) - y0_all(:,3)';

r0_sq = x01.^2 + x02.^2 + x03.^2;
x11 = x01 + v1_r;
x12 = x02 + v2_r;
x13 = x03 + v3_r;
r1_sq = x11.^2 + x12.^2 + x13.^2;

R0 = sqrt(r0_sq + e2);
R1 = sqrt(r1_sq + e2);

d0 = x01 .* v1_r + x02 .* v2_r + x03 .* v3_r;
d1 = d0 + L2_r;

H0 = max(L_r .* R0 + d0, 1e-300);
H1 = max(L_r .* R1 + d1, 1e-300);

inv_L2 = 1 ./ L2_r;
d0_L2  = d0 .* inv_L2;

T0m3 = 1./(R0 .* H0) - 1./(R1 .* H1);

d_sq = L2_r .* (r0_sq + e2) - d0.^2;
d_sq = max(d_sq, e2 .* L2_r);
inv_d_sq = 1 ./ d_sq;

R0_3 = R0.^3;
R1_3 = R1.^3;
T0m5 = inv_d_sq / 3 .* (d1 ./ R1_3 - d0 ./ R0_3) ...
     + 2 * L2_r .* inv_d_sq / 3 .* T0m3;
T1m5 = (1./R0_3 - 1./R1_3) .* (inv_L2 / 3) - d0_L2 .* T0m5;
T2m5 = -1./(3 * R1_3 .* L2_r) + T0m3 .* (inv_L2 / 3) - d0_L2 .* T1m5;

R0_5 = R0.^5;
R1_5 = R1.^5;
T0m7 = inv_d_sq / 5 .* (d1 ./ R1_5 - d0 ./ R0_5) ...
     + 4 * L2_r .* inv_d_sq / 5 .* T0m5;
T1m7 = (1./R0_5 - 1./R1_5) .* (inv_L2 / 5) - d0_L2 .* T0m7;
T2m7 = -1./(5 * R1_5 .* L2_r) + T0m5 .* (inv_L2 / 5) - d0_L2 .* T1m7;

P0 = -2 * T0m3 - 3 * e2 * T0m5 + 15 * e4 * T0m7;
S0 =  6 * T0m5 + 15 * e2 * T0m7;
S1 =  6 * T1m5 + 15 * e2 * T1m7;
S2 =  6 * T2m5 + 15 * e2 * T2m7;

ei_x = d3_all(:,1);  ei_y = d3_all(:,2);  ei_z = d3_all(:,3);
ej_x = d3_all(:,1)'; ej_y = d3_all(:,2)'; ej_z = d3_all(:,3)';

ei_dot_ej = ei_x * ej_x + ei_y * ej_y + ei_z * ej_z;
ei_dot_x0 = ei_x .* x01 + ei_y .* x02 + ei_z .* x03;
ej_dot_x0 = x01 .* ej_x + x02 .* ej_y + x03 .* ej_z;
ei_dot_v  = ei_x .* v1_r + ei_y .* v2_r + ei_z .* v3_r;
ej_dot_v  = ej_x .* v1_r + ej_y .* v2_r + ej_z .* v3_r;

M_twist = (L_r / 4) .* ( ...
    ei_dot_ej .* P0 ...
  + ei_dot_x0 .* ej_dot_x0 .* S0 ...
  + (ei_dot_x0 .* ej_dot_v + ei_dot_v .* ej_dot_x0) .* S1 ...
  + ei_dot_v .* ej_dot_v .* S2);

% The operator maps torque density to angular velocity, so reciprocity uses
% the actual segment-length inner product.  This reduces to ordinary
% symmetrisation when all segment lengths are equal.
M_twist = weighted_self_adjoint(M_twist, L_all);
end
