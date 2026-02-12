function [M_weighted, A_twist_cells] = build_global_resistance(X_all, params, geom)
%BUILD_GLOBAL_RESISTANCE  Global hydrodynamic mobility and torsional resistance
%
%  Description:
%    Builds the translational RSM weighted mobility matrix for ALL points
%    (body + pili) and the torsional RFT resistance for pili segments.
%
%    Translational:  Returns M_weighted (NOT its inverse).  The caller
%    computes M_weighted \ Q directly, avoiding O(N^3) full inversion.
%
%    Torsional:  A_twist{k} = eps_k^2 * Atilde_k / 2  (diagonal per pili).
%
%    Adaptive epsilon strategy:
%      Body points  :  eps_body = min(NN_distance)/2
%      Pili points  :  eps_pili = r_pili / L_0
%      Pair rule    :  eps_ij = eps_j  (force-point-determined)
%
%    Walker (2020) consistency:  Uniform column weights ds_k for ALL pili
%    nodes.  Trapezoidal correction applied in B_global (ode_coupled_rhs).
%
%  Syntax:
%    [M_weighted, A_twist_cells] = build_global_resistance(X_all, params, geom)
%
%  Inputs:
%    X_all  - [N_total x 3] All point positions in lab frame
%    params - Parameter structure (from init_coupled_params)
%    geom   - Geometry structure (from reconstruct_geometry)
%
%  Outputs:
%    M_weighted    - [3*N_total x 3*N_total] Weighted mobility matrix
%    A_twist_cells - {M x 1} Cell array, each [3*N_k x N_k] torsional map
%
%  References:
%    [1] Cortez et al. (2005), Phys. Fluids 17, 031504.
%    [2] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Standardised formatting

N_total = size(X_all, 1);
N_body  = params.N_body;
M       = params.M_pili;
N_segs  = params.N_segments;

%% Build column weight vector (Walker convention)
w = ones(N_total, 1);
pt_idx = N_body + 1;
for k = 1:M
    N_k   = N_segs(k);
    ds_k  = params.ell(k) / N_k;
    n_pts = N_k + 1;
    w(pt_idx : pt_idx + n_pts - 1) = ds_k;
    pt_idx = pt_idx + n_pts;
end

%% Build per-point regularisation parameter vector
eps_vec = zeros(N_total, 1);
eps_vec(1:N_body) = params.eps_body;

pt_idx   = N_body + 1;
eps_pili = params.epsilon;
for k = 1:M
    N_k   = N_segs(k);
    n_pts = N_k + 1;
    eps_vec(pt_idx : pt_idx + n_pts - 1) = eps_pili;
    pt_idx = pt_idx + n_pts;
end

%% Pairwise regularisation: eps_ij = eps_j
eps_pair = repmat(eps_vec', [N_total, 1]);
e2_pair  = eps_pair.^2;

%% Build RSM mobility matrix S (vectorised)
dim = 3 * N_total;
S   = zeros(dim, dim);

X_i = reshape(X_all, [N_total, 1, 3]);
X_j = reshape(X_all, [1, N_total, 3]);
r_vec = X_i - X_j;
rx = r_vec(:,:,1);
ry = r_vec(:,:,2);
rz = r_vec(:,:,3);
r2 = rx.^2 + ry.^2 + rz.^2;

denom = (r2 + e2_pair).^(1.5);
H1    = (r2 + 2 * e2_pair) ./ denom;
H2    = 1 ./ denom;

idx_x = 1:3:dim;
idx_y = 2:3:dim;
idx_z = 3:3:dim;

S(idx_x, idx_x) = H1 + H2 .* (rx .* rx);
S(idx_x, idx_y) =       H2 .* (rx .* ry);
S(idx_x, idx_z) =       H2 .* (rx .* rz);
S(idx_y, idx_x) =       H2 .* (ry .* rx);
S(idx_y, idx_y) = H1 + H2 .* (ry .* ry);
S(idx_y, idx_z) =       H2 .* (ry .* rz);
S(idx_z, idx_x) =       H2 .* (rz .* rx);
S(idx_z, idx_y) =       H2 .* (rz .* ry);
S(idx_z, idx_z) = H1 + H2 .* (rz .* rz);

%% Apply column weighting: M_weighted = S * diag(w ⊗ [1;1;1])
W_vec      = repelem(w, 3);
M_weighted = S .* W_vec';

%% Build torsional resistance for each pilus
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
