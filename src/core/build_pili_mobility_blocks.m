function [M_PB, M_BP, M_PP, C_twist, M_twist, L_all] = build_pili_mobility_blocks( ...
        X_all_bf, d3_all_bf, params)
%BUILD_PILI_MOBILITY_BLOCKS  Body-frame non-BB hydrodynamic mobility blocks
%
%  Description:
%    Builds only the configuration-dependent hydrodynamic blocks needed by
%    the body-frame block Schur complement solver in build_coupled_lin_sys:
%      - M_PB     [3*N_pili x 3*N_body]  RSM body -> pili (eps_body)
%      - M_BP     [3*N_body x 3*N_pili]  RSS pili -> body (eps_pili)
%      - M_PP     [3*N_pili x 3*N_pili]  RSS pili -> pili (eps_pili)
%      - C_twist  [3*N_total x N_seg]    0.5 * M_rotlet * D_twist (direct)
%      - M_twist  [N_seg x N_seg]        non-local axial twist mobility
%      - L_all    [N_seg x 1]            current segment lengths
%
%    Skips the rigid body-body block M_BB (precomputed once in body frame
%    via build_body_precomp) and contracts the rotlet directly against d3
%    without allocating the dense M_rotlet or D_twist matrices.
%
%  Syntax:
%    [M_PB, M_BP, M_PP, C_twist, M_twist, L_all] = ...
%        build_pili_mobility_blocks(X_all_bf, d3_all_bf, params)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

N_total  = size(X_all_bf, 1);
N_body   = params.N_body;
M        = params.M_pili;
N_segs   = params.N_segments;
N_pili   = N_total - N_body;
nB       = 3 * N_body;
nP       = 3 * N_pili;
dim      = 3 * N_total;

e2_body  = params.eps_body^2;
eps_pili = params.epsilon;
e2_pili  = eps_pili^2;

%% ---- M_PB: body sources -> pili eval points (RSM, eps = eps_body) ----
X_pili_bf = X_all_bf(N_body + 1 : end, :);
X_body_bf = X_all_bf(1 : N_body, :);

rx_pb = X_pili_bf(:, 1) - X_body_bf(:, 1)';
ry_pb = X_pili_bf(:, 2) - X_body_bf(:, 2)';
rz_pb = X_pili_bf(:, 3) - X_body_bf(:, 3)';
r2_pb = rx_pb.^2 + ry_pb.^2 + rz_pb.^2;

inv_R_pb = 1 ./ sqrt(r2_pb + e2_body);
H2_pb    = inv_R_pb .* inv_R_pb .* inv_R_pb;
H1_pb    = (r2_pb + 2 * e2_body) .* H2_pb;

M_PB  = zeros(nP, nB);
rpx   = 1:3:nP;   rpy = 2:3:nP;   rpz = 3:3:nP;
cbx   = 1:3:nB;   cby = 2:3:nB;   cbz = 3:3:nB;
G_xy  = H2_pb .* (rx_pb .* ry_pb);
G_xz  = H2_pb .* (rx_pb .* rz_pb);
G_yz  = H2_pb .* (ry_pb .* rz_pb);

M_PB(rpx, cbx) = H1_pb + H2_pb .* (rx_pb .* rx_pb);
M_PB(rpx, cby) = G_xy;
M_PB(rpx, cbz) = G_xz;
M_PB(rpy, cbx) = G_xy;
M_PB(rpy, cby) = H1_pb + H2_pb .* (ry_pb .* ry_pb);
M_PB(rpy, cbz) = G_yz;
M_PB(rpz, cbx) = G_xz;
M_PB(rpz, cby) = G_yz;
M_PB(rpz, cbz) = H1_pb + H2_pb .* (rz_pb .* rz_pb);

%% ---- Collect pili segment endpoints in body frame ----
N_seg_all = sum(N_segs);
y0_all = zeros(N_seg_all, 3);
y1_all = zeros(N_seg_all, 3);

s_count     = 0;
pili_offset = N_body;
for k = 1:M
    N_k   = N_segs(k);
    idx_s = s_count + (1:N_k);
    y0_all(idx_s, :) = X_all_bf(pili_offset + (1:N_k), :);
    y1_all(idx_s, :) = X_all_bf(pili_offset + (2:N_k+1), :);
    s_count     = s_count + N_k;
    pili_offset = pili_offset + N_k + 1;
end

%% ---- Segment geometry [1 x N_seg_all] ----
v_all  = y0_all - y1_all;
L2_all = sum(v_all.^2, 2);
L_all  = sqrt(L2_all);

L_r  = L_all';
L2_r = L2_all';
v1_r = v_all(:, 1)';
v2_r = v_all(:, 2)';
v3_r = v_all(:, 3)';

%% ---- Displacements [N_total x N_seg_all] ----
x01 = X_all_bf(:, 1) - y0_all(:, 1)';
x02 = X_all_bf(:, 2) - y0_all(:, 2)';
x03 = X_all_bf(:, 3) - y0_all(:, 3)';

r0_sq = x01.^2 + x02.^2 + x03.^2;

x11 = x01 + v1_r;
x12 = x02 + v2_r;
x13 = x03 + v3_r;
r1_sq = x11.^2 + x12.^2 + x13.^2;

R0 = sqrt(r0_sq + e2_pili);
R1 = sqrt(r1_sq + e2_pili);
inv_R0 = 1 ./ R0;
inv_R1 = 1 ./ R1;

d0 = x01 .* v1_r + x02 .* v2_r + x03 .* v3_r;

H0_v = max(L_r .* R0 + d0,        1e-300);
H1_v = max(L_r .* R1 + d0 + L2_r, 1e-300);

%% ---- T_{n,q} recursion: q = -1, -3, -5 (Cortez 2018) ----
inv_L  = 1 ./ L_r;
inv_L2 = 1 ./ L2_r;
d0_L2  = d0 .* inv_L2;

T0m1 = log(H1_v ./ H0_v) .* inv_L;
T0m3 = inv_R0 ./ H0_v - inv_R1 ./ H1_v;
T1m1 = (R1 - R0) .* inv_L2 - d0_L2 .* T0m1;
T1m3 = (inv_R0 - inv_R1) .* inv_L2 - d0_L2 .* T0m3;
neg_inv_R1_L2 = -inv_R1 .* inv_L2;
T2m3 = neg_inv_R1_L2 + T0m1 .* inv_L2 - d0_L2 .* T1m3;
T3m3 = neg_inv_R1_L2 + 2 * T1m1 .* inv_L2 - d0_L2 .* T2m3;

inv_R0_3 = inv_R0 .* inv_R0 .* inv_R0;
inv_R1_3 = inv_R1 .* inv_R1 .* inv_R1;
d_sq = max(L2_r .* (r0_sq + e2_pili) - d0.^2, e2_pili .* L2_r);
inv_d_sq_3 = (1 / 3) ./ d_sq;
d1_v   = d0 + L2_r;
inv3L2 = inv_L2 / 3;

T0m5 = inv_d_sq_3 .* (d1_v .* inv_R1_3 - d0 .* inv_R0_3 + 2 * L2_r .* T0m3);
T1m5 = (inv_R0_3 - inv_R1_3) .* inv3L2 - d0_L2 .* T0m5;

%% ---- Stokeslet M2 & M1 (3x3 symmetric per pair) ----
c_iso  = L_r .* (T1m1 + e2_pili * T1m3);
c1s    = L_r .* T1m3;
c2s    = L_r .* T2m3;
c3s    = L_r .* T3m3;

c_iso0 = L_r .* (T0m1 + e2_pili * T0m3) - c_iso;
c10    = L_r .* T0m3 - c1s;
c20    = c1s - c2s;
c30    = c2s - c3s;

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

%% ---- Scatter pili columns into [M_BP; M_PP] (in-place per pilus) ----
M_Pcols = zeros(dim, nP);
row_x = 1:3:dim;  row_y = 2:3:dim;  row_z = 3:3:dim;

seg_offset  = 0;
node_offset = 0;
for k = 1:M
    N_k     = N_segs(k);
    seg_idx = seg_offset + (1:N_k);
    col_base = 3 * node_offset;
    cx1 = col_base + (1 : 3 : 3*N_k);
    cy1 = cx1 + 1;
    cz1 = cx1 + 2;
    cx2 = cx1 + 3;
    cy2 = cy1 + 3;
    cz2 = cz1 + 3;

    M_Pcols(row_x, cx1) = M1_11(:, seg_idx);
    M_Pcols(row_x, cx2) = M_Pcols(row_x, cx2) + M2_11(:, seg_idx);

    M_Pcols(row_y, cy1) = M1_22(:, seg_idx);
    M_Pcols(row_y, cy2) = M_Pcols(row_y, cy2) + M2_22(:, seg_idx);

    M_Pcols(row_z, cz1) = M1_33(:, seg_idx);
    M_Pcols(row_z, cz2) = M_Pcols(row_z, cz2) + M2_33(:, seg_idx);

    M_Pcols(row_x, cy1) = M1_12(:, seg_idx);
    M_Pcols(row_x, cy2) = M_Pcols(row_x, cy2) + M2_12(:, seg_idx);
    M_Pcols(row_y, cx1) = M_Pcols(row_x, cy1);
    M_Pcols(row_y, cx2(end)) = M_Pcols(row_x, cy2(end));

    M_Pcols(row_x, cz1) = M1_13(:, seg_idx);
    M_Pcols(row_x, cz2) = M_Pcols(row_x, cz2) + M2_13(:, seg_idx);
    M_Pcols(row_z, cx1) = M_Pcols(row_x, cz1);
    M_Pcols(row_z, cx2(end)) = M_Pcols(row_x, cz2(end));

    M_Pcols(row_y, cz1) = M1_23(:, seg_idx);
    M_Pcols(row_y, cz2) = M_Pcols(row_y, cz2) + M2_23(:, seg_idx);
    M_Pcols(row_z, cy1) = M_Pcols(row_y, cz1);
    M_Pcols(row_z, cy2(end)) = M_Pcols(row_y, cz2(end));

    node_offset = node_offset + N_k + 1;
    seg_offset  = seg_offset  + N_k;
end

M_BP = M_Pcols(1:nB, :);
M_PP = M_Pcols(nB+1:end, :);

%% ---- Direct C_twist = 0.5 * M_rotlet * D_twist (no M_rotlet/D_twist alloc) ----
ca  = 2 * T0m3 + 3 * e2_pili * T0m5;
cb  = 2 * T1m3 + 3 * e2_pili * T1m5;
Rxy =  L_r .* (ca .* x03 + cb .* v3_r);
Rxz = -L_r .* (ca .* x02 + cb .* v2_r);
Ryz =  L_r .* (ca .* x01 + cb .* v1_r);

d3x = d3_all_bf(:, 1)';
d3y = d3_all_bf(:, 2)';
d3z = d3_all_bf(:, 3)';

C_twist = zeros(dim, N_seg_all);
C_twist(row_x, :) = 0.5 * ( Rxy .* d3y + Rxz .* d3z);
C_twist(row_y, :) = 0.5 * (-Rxy .* d3x + Ryz .* d3z);
C_twist(row_z, :) = 0.5 * (-Rxz .* d3x - Ryz .* d3y);

%% ---- Full non-local axial twist mobility M_twist ----
M_twist = compute_nonlocal_twist_mobility_bf(y0_all, y1_all, d3_all_bf, eps_pili);

end


function M_twist = compute_nonlocal_twist_mobility_bf(y0_all, y1_all, d3_all, eps_pili)
%COMPUTE_NONLOCAL_TWIST_MOBILITY_BF  Segment-integral axial rotation mobility
mids   = 0.5 * (y0_all + y1_all);
v_all  = y0_all - y1_all;
L2_all = sum(v_all.^2, 2);
L_all  = sqrt(L2_all);

e2 = eps_pili^2;
e4 = e2^2;

L_r  = L_all';
L2_r = L2_all';
v1_r = v_all(:, 1)';
v2_r = v_all(:, 2)';
v3_r = v_all(:, 3)';

x01 = mids(:, 1) - y0_all(:, 1)';
x02 = mids(:, 2) - y0_all(:, 2)';
x03 = mids(:, 3) - y0_all(:, 3)';

r0_sq = x01.^2 + x02.^2 + x03.^2;
x11 = x01 + v1_r;  x12 = x02 + v2_r;  x13 = x03 + v3_r;
r1_sq = x11.^2 + x12.^2 + x13.^2;

R0 = sqrt(r0_sq + e2);
R1 = sqrt(r1_sq + e2);
inv_R0 = 1 ./ R0;
inv_R1 = 1 ./ R1;

d0 = x01 .* v1_r + x02 .* v2_r + x03 .* v3_r;
d1 = d0 + L2_r;

H0 = max(L_r .* R0 + d0, 1e-300);
H1 = max(L_r .* R1 + d1, 1e-300);

inv_L2 = 1 ./ L2_r;
d0_L2  = d0 .* inv_L2;

T0m3 = inv_R0 ./ H0 - inv_R1 ./ H1;

d_sq = max(L2_r .* (r0_sq + e2) - d0.^2, e2 .* L2_r);
inv_d_sq = 1 ./ d_sq;

inv_R0_3 = inv_R0 .* inv_R0 .* inv_R0;
inv_R1_3 = inv_R1 .* inv_R1 .* inv_R1;
T0m5 = (inv_d_sq / 3) .* (d1 .* inv_R1_3 - d0 .* inv_R0_3 + 2 * L2_r .* T0m3);
T1m5 = (inv_R0_3 - inv_R1_3) .* (inv_L2 / 3) - d0_L2 .* T0m5;
T2m5 = -inv_R1_3 .* (inv_L2 / 3) + T0m3 .* (inv_L2 / 3) - d0_L2 .* T1m5;

inv_R0_5 = inv_R0_3 .* inv_R0 .* inv_R0;
inv_R1_5 = inv_R1_3 .* inv_R1 .* inv_R1;
T0m7 = (inv_d_sq / 5) .* (d1 .* inv_R1_5 - d0 .* inv_R0_5 + 4 * L2_r .* T0m5);
T1m7 = (inv_R0_5 - inv_R1_5) .* (inv_L2 / 5) - d0_L2 .* T0m7;
T2m7 = -inv_R1_5 .* (inv_L2 / 5) + T0m5 .* (inv_L2 / 5) - d0_L2 .* T1m7;

P0 = -2 * T0m3 - 3 * e2 * T0m5 + 15 * e4 * T0m7;
S0 =  6 * T0m5 + 15 * e2 * T0m7;
S1 =  6 * T1m5 + 15 * e2 * T1m7;
S2 =  6 * T2m5 + 15 * e2 * T2m7;

ei_x = d3_all(:, 1);  ei_y = d3_all(:, 2);  ei_z = d3_all(:, 3);
ej_x = d3_all(:, 1)'; ej_y = d3_all(:, 2)'; ej_z = d3_all(:, 3)';

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

M_twist = weighted_self_adjoint(M_twist, L_all);
end
