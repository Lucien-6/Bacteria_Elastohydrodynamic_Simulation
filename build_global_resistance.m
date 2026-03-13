function [M_weighted, A_twist_cells, M_rotlet] = build_global_resistance(X_all, params, geom)
%BUILD_GLOBAL_RESISTANCE  Global hydrodynamic mobility and torsional resistance
%
%  Description:
%    Builds the translational mobility matrix for ALL points (body + pili)
%    and the torsional RFT resistance for pili segments.
%
%    Mobility method:
%      Body columns :  RSM point forces (Cortez 2005)
%      Pili columns :  RSS exact segment integrals (Cortez 2018)
%      Rotlet       :  RSS rotlet segment integrals (Cortez 2018)
%
%    Translational:  Returns M_weighted (NOT its inverse).  The caller
%    computes M_weighted \ Q directly, avoiding O(N^3) full inversion.
%
%    Torsional:  A_twist{k} = ell_k * eps_k^2 * Atilde_k / 2
%                (diagonal per pili, RFT).
%
%    Rotlet segment mobility (Cortez 2018):
%      M_rotlet maps 3D torque densities (per segment) to translational
%      velocities (per eval point).  This captures the hydrodynamic
%      coupling between torsion and translation that RFT neglects.
%
%    Adaptive epsilon strategy:
%      Body points  :  eps_body = h_mean / 3   (mean nearest-neighbour / 3,
%                      empirically optimal; see init_coupled_params.m)
%      Pili points  :  eps_pili = r_pili / L_0
%      Pair rule    :  eps_ij = eps_j  (force-point-determined)
%
%    Matrix layout (block view):
%      M_weighted = [ M_BB  M_BP ]   BB = body-body   (RSM)
%                   [ M_PB  M_PP ]   BP = body-pili   (RSS)
%                                    PB = pili-body   (RSM)
%                                    PP = pili-pili   (RSS)
%
%    Body columns:  RSM point forces,  eps = eps_body,  w = 1
%    Pili columns:  RSS segment integrals, eps = eps_pili, w = 1
%      (RSS integrals already incorporate the segment-length weight)
%
%  Syntax:
%    [M_weighted, A_twist_cells, M_rotlet] = build_global_resistance(X_all, params, geom)
%
%  Inputs:
%    X_all  - [N_total x 3] All point positions in lab frame
%    params - Parameter structure (from init_coupled_params)
%    geom   - Geometry structure (from reconstruct_geometry)
%
%  Outputs:
%    M_weighted    - [3*N_total x 3*N_total] Weighted mobility matrix
%    A_twist_cells - {M x 1} Cell array, each [3*N_k x N_k] torsional map
%    M_rotlet      - [3*N_total x 3*N_seg_all] Rotlet segment mobility
%
%  References:
%    [1] Cortez et al. (2005), Phys. Fluids 17, 031504.
%    [2] Cortez (2018), J. Comput. Phys. 375, 783-796.
%    [3] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-14 — Pure RSS + Rotlet (RSM removed for pili)

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
X_eval = reshape(X_all,              [N_total, 1,      3]);
X_body = reshape(X_all(1:N_body, :), [1,       N_body, 3]);

r_bb  = X_eval - X_body;
rx_bb = r_bb(:,:,1);
ry_bb = r_bb(:,:,2);
rz_bb = r_bb(:,:,3);
r2_bb = rx_bb.^2 + ry_bb.^2 + rz_bb.^2;

denom_bb = (r2_bb + e2_body).^(1.5);
H1_bb    = (r2_bb + 2 * e2_body) ./ denom_bb;
H2_bb    = 1 ./ denom_bb;

row_x = 1:3:dim;           % N_total elements
row_y = 2:3:dim;
row_z = 3:3:dim;
col_x = 1:3:3*N_body;      % N_body elements
col_y = 2:3:3*N_body;
col_z = 3:3:3*N_body;

M_weighted(row_x, col_x) = H1_bb + H2_bb .* (rx_bb .* rx_bb);
M_weighted(row_x, col_y) =         H2_bb .* (rx_bb .* ry_bb);
M_weighted(row_x, col_z) =         H2_bb .* (rx_bb .* rz_bb);
M_weighted(row_y, col_x) =         H2_bb .* (ry_bb .* rx_bb);
M_weighted(row_y, col_y) = H1_bb + H2_bb .* (ry_bb .* ry_bb);
M_weighted(row_y, col_z) =         H2_bb .* (ry_bb .* rz_bb);
M_weighted(row_z, col_x) =         H2_bb .* (rz_bb .* rx_bb);
M_weighted(row_z, col_y) =         H2_bb .* (rz_bb .* ry_bb);
M_weighted(row_z, col_z) = H1_bb + H2_bb .* (rz_bb .* rz_bb);

%% ---- Pre-collect all pili segment endpoints ----
N_seg_all = sum(N_segs);
y0_all = zeros(N_seg_all, 3);     % Segment start (node k)
y1_all = zeros(N_seg_all, 3);     % Segment end   (node k+1)
g0_idx = zeros(N_seg_all, 1);     % Global node index (start)
g1_idx = zeros(N_seg_all, 1);     % Global node index (end)

s_count     = 0;
pili_offset = N_body;
for k = 1:M
    N_k = N_segs(k);
    for s = 1:N_k
        s_count = s_count + 1;
        g0 = pili_offset + s;
        g1 = pili_offset + s + 1;
        g0_idx(s_count) = g0;
        g1_idx(s_count) = g1;
        y0_all(s_count, :) = X_all(g0, :);
        y1_all(s_count, :) = X_all(g1, :);
    end
    pili_offset = pili_offset + N_k + 1;
end

%% ---- Vectorised segment geometry  [1 x N_seg_all] ----
v_all  = y0_all - y1_all;                        % [N_seg x 3]
L_all  = sqrt(sum(v_all.^2, 2));                 % [N_seg x 1]
L2_all = L_all .* L_all;

% Transpose to row vectors for broadcasting [1 x N_seg]
L_r   = L_all';
L2_r  = L2_all';
v1_r  = v_all(:,1)';
v2_r  = v_all(:,2)';
v3_r  = v_all(:,3)';

%% ---- Displacements  [N_total x N_seg_all] ----
x01 = bsxfun(@minus, X_all(:,1), y0_all(:,1)');
x02 = bsxfun(@minus, X_all(:,2), y0_all(:,2)');
x03 = bsxfun(@minus, X_all(:,3), y0_all(:,3)');

r0_sq = x01.^2 + x02.^2 + x03.^2;

x11 = x01 + v1_r;
x12 = x02 + v2_r;
x13 = x03 + v3_r;
r1_sq = x11.^2 + x12.^2 + x13.^2;

R0 = sqrt(r0_sq + e2_pili);
R1 = sqrt(r1_sq + e2_pili);

d0 = x01 .* v1_r + x02 .* v2_r + x03 .* v3_r;  % x0 . v

%% ---- Log arguments (strictly positive for eps > 0) ----
H0_v = L_r .* R0 + d0;
H1_v = L_r .* R1 + d0 + L2_r;
H0_v = max(H0_v, 1e-300);
H1_v = max(H1_v, 1e-300);

%% ---- T_{n,q} recursion: q = -1, -3  (Cortez 2018, Eq. 7) ----
inv_L  = 1 ./ L_r;
inv_L2 = 1 ./ L2_r;
d0_L2  = d0 .* inv_L2;                           % d0 / L^2

T0m1 = (log(H1_v) - log(H0_v)) .* inv_L;
T0m3 = 1./(R0 .* H0_v) - 1./(R1 .* H1_v);
T1m1 = (R1 - R0) .* inv_L2 - d0_L2 .* T0m1;
T1m3 = (1./R0 - 1./R1) .* inv_L2 - d0_L2 .* T0m3;
T2m3 = -1./(R1 .* L2_r) + T0m1 .* inv_L2 - d0_L2 .* T1m3;
T3m3 = -1./(R1 .* L2_r) + 2 * T1m1 .* inv_L2 - d0_L2 .* T2m3;

%% ---- T_{n,q} recursion: q = -5  (Cortez 2018, Sec. 2.4) ----
%  Needed for the regularised rotlet segment.
R0_3 = R0.^3;
R1_3 = R1.^3;

d_sq = L2_r .* (r0_sq + e2_pili) - d0.^2;       % d^2
d_sq = max(d_sq, e2_pili .* L2_r);               % safeguard: d^2 >= L^2 eps^2
inv_d_sq = 1 ./ d_sq;

d1_v   = d0 + L2_r;                              % x1 . v  =  x0.v + L^2
inv3L2 = inv_L2 / 3;

T0m5 = inv_d_sq / 3 .* (d1_v ./ R1_3 - d0 ./ R0_3) ...
     + 2 * L2_r .* inv_d_sq / 3 .* T0m3;
T1m5 = (1./R0_3 - 1./R1_3) .* inv3L2 - d0_L2 .* T0m5;

%% ---- Stokeslet: M2 components  (Cortez 2018, Eq. 12) ----
c_iso = L_r .* (T1m1 + e2_pili * T1m3);
c1s   = L_r .* T1m3;
c2s   = L_r .* T2m3;
c3s   = L_r .* T3m3;

M2_11 = c_iso + c1s.*(x01.*x01) + c2s.*(2*x01.*v1_r)         + c3s.*(v1_r.*v1_r);
M2_12 =         c1s.*(x01.*x02) + c2s.*(x01.*v2_r+x02.*v1_r) + c3s.*(v1_r.*v2_r);
M2_13 =         c1s.*(x01.*x03) + c2s.*(x01.*v3_r+x03.*v1_r) + c3s.*(v1_r.*v3_r);
M2_21 =         c1s.*(x02.*x01) + c2s.*(x02.*v1_r+x01.*v2_r) + c3s.*(v2_r.*v1_r);
M2_22 = c_iso + c1s.*(x02.*x02) + c2s.*(2*x02.*v2_r)         + c3s.*(v2_r.*v2_r);
M2_23 =         c1s.*(x02.*x03) + c2s.*(x02.*v3_r+x03.*v2_r) + c3s.*(v2_r.*v3_r);
M2_31 =         c1s.*(x03.*x01) + c2s.*(x03.*v1_r+x01.*v3_r) + c3s.*(v3_r.*v1_r);
M2_32 =         c1s.*(x03.*x02) + c2s.*(x03.*v2_r+x02.*v3_r) + c3s.*(v3_r.*v2_r);
M2_33 = c_iso + c1s.*(x03.*x03) + c2s.*(2*x03.*v3_r)         + c3s.*(v3_r.*v3_r);

%% ---- Stokeslet: M1 = L*(T0-block) - M2  (Cortez 2018, Eq. 13) ----
c_iso0 = L_r .* (T0m1 + e2_pili * T0m3);
c10    = L_r .* T0m3;
c20    = L_r .* T1m3;
c30    = L_r .* T2m3;

M1_11 = c_iso0 + c10.*(x01.*x01) + c20.*(2*x01.*v1_r)         + c30.*(v1_r.*v1_r) - M2_11;
M1_12 =          c10.*(x01.*x02) + c20.*(x01.*v2_r+x02.*v1_r) + c30.*(v1_r.*v2_r) - M2_12;
M1_13 =          c10.*(x01.*x03) + c20.*(x01.*v3_r+x03.*v1_r) + c30.*(v1_r.*v3_r) - M2_13;
M1_21 =          c10.*(x02.*x01) + c20.*(x02.*v1_r+x01.*v2_r) + c30.*(v2_r.*v1_r) - M2_21;
M1_22 = c_iso0 + c10.*(x02.*x02) + c20.*(2*x02.*v2_r)         + c30.*(v2_r.*v2_r) - M2_22;
M1_23 =          c10.*(x02.*x03) + c20.*(x02.*v3_r+x03.*v2_r) + c30.*(v2_r.*v3_r) - M2_23;
M1_31 =          c10.*(x03.*x01) + c20.*(x03.*v1_r+x01.*v3_r) + c30.*(v3_r.*v1_r) - M2_31;
M1_32 =          c10.*(x03.*x02) + c20.*(x03.*v2_r+x02.*v3_r) + c30.*(v3_r.*v2_r) - M2_32;
M1_33 = c_iso0 + c10.*(x03.*x03) + c20.*(2*x03.*v3_r)         + c30.*(v3_r.*v3_r) - M2_33;

%% ---- Scatter Stokeslet M1/M2 into M_weighted ----
for s = 1:N_seg_all
    g0 = g0_idx(s);
    g1 = g1_idx(s);

    c0_1 = 3*(g0-1)+1;  c0_2 = 3*g0-1;  c0_3 = 3*g0;
    c1_1 = 3*(g1-1)+1;  c1_2 = 3*g1-1;  c1_3 = 3*g1;

    % M1 -> column block of node g0
    M_weighted(row_x, c0_1) = M_weighted(row_x, c0_1) + M1_11(:,s);
    M_weighted(row_x, c0_2) = M_weighted(row_x, c0_2) + M1_12(:,s);
    M_weighted(row_x, c0_3) = M_weighted(row_x, c0_3) + M1_13(:,s);
    M_weighted(row_y, c0_1) = M_weighted(row_y, c0_1) + M1_21(:,s);
    M_weighted(row_y, c0_2) = M_weighted(row_y, c0_2) + M1_22(:,s);
    M_weighted(row_y, c0_3) = M_weighted(row_y, c0_3) + M1_23(:,s);
    M_weighted(row_z, c0_1) = M_weighted(row_z, c0_1) + M1_31(:,s);
    M_weighted(row_z, c0_2) = M_weighted(row_z, c0_2) + M1_32(:,s);
    M_weighted(row_z, c0_3) = M_weighted(row_z, c0_3) + M1_33(:,s);

    % M2 -> column block of node g1
    M_weighted(row_x, c1_1) = M_weighted(row_x, c1_1) + M2_11(:,s);
    M_weighted(row_x, c1_2) = M_weighted(row_x, c1_2) + M2_12(:,s);
    M_weighted(row_x, c1_3) = M_weighted(row_x, c1_3) + M2_13(:,s);
    M_weighted(row_y, c1_1) = M_weighted(row_y, c1_1) + M2_21(:,s);
    M_weighted(row_y, c1_2) = M_weighted(row_y, c1_2) + M2_22(:,s);
    M_weighted(row_y, c1_3) = M_weighted(row_y, c1_3) + M2_23(:,s);
    M_weighted(row_z, c1_1) = M_weighted(row_z, c1_1) + M2_31(:,s);
    M_weighted(row_z, c1_2) = M_weighted(row_z, c1_2) + M2_32(:,s);
    M_weighted(row_z, c1_3) = M_weighted(row_z, c1_3) + M2_33(:,s);
end

%% ---- Rotlet segment mobility  (Cortez 2018, Sec. 2.5) ----
M_rotlet = zeros(dim, 3 * N_seg_all);

ca = 2 * T0m3 + 3 * e2_pili * T0m5;              % [N_total x N_seg_all]
cb = 2 * T1m3 + 3 * e2_pili * T1m5;

Rxy =  L_r .* (ca .* x03 + cb .* v3_r);          % R12
Rxz = -L_r .* (ca .* x02 + cb .* v2_r);          % R13
Ryz =  L_r .* (ca .* x01 + cb .* v1_r);          % R23

for s = 1:N_seg_all
    tc1 = 3*(s-1) + 1;
    tc2 = 3*(s-1) + 2;
    tc3 = 3*(s-1) + 3;

    M_rotlet(row_x, tc2) =  Rxy(:,s);
    M_rotlet(row_x, tc3) =  Rxz(:,s);
    M_rotlet(row_y, tc1) = -Rxy(:,s);
    M_rotlet(row_y, tc3) =  Ryz(:,s);
    M_rotlet(row_z, tc1) = -Rxz(:,s);
    M_rotlet(row_z, tc2) = -Ryz(:,s);
end

%% Build torsional resistance for each pilus (RFT)
A_twist_cells = cell(M, 1);

for k = 1:M
    N_k   = N_segs(k);
    eps_k = params.epsilon / params.ell(k);

    d3_k = geom.d3{k};                  % [N_k x 3] in lab frame

    Atilde_k = zeros(3*N_k, N_k);
    for j = 1:N_k
        Atilde_k(3*(j-1)+1 : 3*j, j) = d3_k(j, :)';
    end

    A_twist_cells{k} = params.ell(k) * eps_k^2 * Atilde_k / 2;
end

end
