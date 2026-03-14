function brown = compute_body_brownian_params(X_body_dim, mu, eps_body_dim, Temper)
%COMPUTE_BODY_BROWNIAN_PARAMS  Compute Brownian motion parameters for a rigid body
%
%  Computes the 6x6 friction matrix of a rigid capsule body from its
%  surface points using the Regularized Stokeslet Method (RSM, Cortez 2005),
%  then extracts decoupled translational and rotational diffusion tensors
%  via Schur complement, and returns their Cholesky factors for
%  Euler-Maruyama Brownian dynamics.
%
%  Inputs:
%    X_body_dim    [N_body x 3]  Body surface point coordinates [m]
%    mu            [scalar]      Dynamic viscosity [Pa*s]
%    eps_body_dim  [scalar]      RSM regularisation parameter [m]
%    Temper        [scalar]      Absolute temperature [K]
%
%  Output:
%    brown  struct with fields:
%      .L_trans    [3x3]  Lower Cholesky factor for translation [m/sqrt(s)]
%      .L_rot      [3x3]  Lower Cholesky factor for rotation [rad/sqrt(s)]
%      .D_trans    [3x3]  Translational diffusion tensor [m^2/s]
%      .D_rot      [3x3]  Rotational diffusion tensor [rad^2/s]
%      .FM         [6x6]  Full body friction matrix
%      .DFM_trans  [3x3]  Decoupled translational friction [N*s/m]
%      .DFM_rot    [3x3]  Decoupled rotational friction [N*m*s]
%
%  References:
%    Cortez et al. (2005), Physics of Fluids, 17, 031504.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-13

N_b = size(X_body_dim, 1);
e2  = eps_body_dim^2;

%% Step 1: Build RSM mobility matrix M_BB [3N_b x 3N_b]
X_eval = reshape(X_body_dim, [N_b, 1, 3]);
X_src  = reshape(X_body_dim, [1, N_b, 3]);
r_all  = X_eval - X_src;
rx = r_all(:,:,1);
ry = r_all(:,:,2);
rz = r_all(:,:,3);
r2 = rx.^2 + ry.^2 + rz.^2;

denom = (r2 + e2).^1.5;
H1    = (r2 + 2*e2) ./ denom;
H2    = 1 ./ denom;

prefac = 1 / (8 * pi * mu);

row_x = 1:3:3*N_b;  row_y = 2:3:3*N_b;  row_z = 3:3:3*N_b;
col_x = 1:3:3*N_b;  col_y = 2:3:3*N_b;  col_z = 3:3:3*N_b;

M_BB = zeros(3*N_b, 3*N_b);
M_BB(row_x, col_x) = prefac * (H1 + H2 .* rx.*rx);
M_BB(row_x, col_y) = prefac * (     H2 .* rx.*ry);
M_BB(row_x, col_z) = prefac * (     H2 .* rx.*rz);
M_BB(row_y, col_x) = prefac * (     H2 .* ry.*rx);
M_BB(row_y, col_y) = prefac * (H1 + H2 .* ry.*ry);
M_BB(row_y, col_z) = prefac * (     H2 .* ry.*rz);
M_BB(row_z, col_x) = prefac * (     H2 .* rz.*rx);
M_BB(row_z, col_y) = prefac * (     H2 .* rz.*ry);
M_BB(row_z, col_z) = prefac * (H1 + H2 .* rz.*rz);

%% Step 2: Build no-slip constraint matrix N_ns [3N_b x 6]
N_ns = zeros(3*N_b, 6);
for i = 1:N_b
    xi   = X_body_dim(i, :)';
    rows = (3*i-2):(3*i);
    N_ns(rows, 1:3) = eye(3);
    N_ns(rows, 4:6) = -skew_matrix(xi);
end

%% Step 3: Build force/torque extraction matrix E [6 x 3N_b]
E_mat = zeros(6, 3*N_b);
for j = 1:N_b
    xj   = X_body_dim(j, :)';
    cols = (3*j-2):(3*j);
    E_mat(1:3, cols) = eye(3);
    E_mat(4:6, cols) = skew_matrix(xj);
end

%% Step 4: Solve FM_body = E * (M_BB \ N_ns)
rc = rcond(M_BB);
if rc < 1e-12
    warning('compute_body_brownian_params:IllConditioned', ...
        'M_BB is ill-conditioned (rcond = %.2e). Friction matrix may be inaccurate.', rc);
end
F_coeff = M_BB \ N_ns;
FM      = E_mat * F_coeff;

%% Step 5: Schur complement decoupling
A = FM(1:3, 1:3);   B = FM(1:3, 4:6);
C = FM(4:6, 1:3);   D = FM(4:6, 4:6);

DFM_trans = A - B / D * C;
DFM_rot   = D - C / A * B;

DFM_trans = 0.5 * (DFM_trans + DFM_trans');
DFM_rot   = 0.5 * (DFM_rot   + DFM_rot');

%% Step 6: Diffusion tensors and Cholesky decomposition
kB_const = 1.380649e-23;

D_trans = kB_const * Temper * inv(DFM_trans);                      %#ok<MINV>
D_rot   = kB_const * Temper * inv(DFM_rot);                        %#ok<MINV>

D_trans = 0.5 * (D_trans + D_trans');
D_rot   = 0.5 * (D_rot   + D_rot');

[~, flag_t] = chol(2 * D_trans);
if flag_t ~= 0
    warning('compute_body_brownian_params:NotPosDef', ...
        'D_trans not positive definite. Regularising...');
    min_eig_t = min(eig(D_trans));
    reg_t     = max(eps(norm(D_trans, 'fro')), abs(min_eig_t) * (1 + 1e-6));
    D_trans   = D_trans + eye(3) * reg_t;
end
[~, flag_r] = chol(2 * D_rot);
if flag_r ~= 0
    warning('compute_body_brownian_params:NotPosDef', ...
        'D_rot not positive definite. Regularising...');
    min_eig_r = min(eig(D_rot));
    reg_r     = max(eps(norm(D_rot, 'fro')), abs(min_eig_r) * (1 + 1e-6));
    D_rot     = D_rot + eye(3) * reg_r;
end

L_trans = chol(2 * D_trans, 'lower');
L_rot   = chol(2 * D_rot,   'lower');

%% Pack output
brown.L_trans   = L_trans;
brown.L_rot     = L_rot;
brown.D_trans   = D_trans;
brown.D_rot     = D_rot;
brown.FM        = FM;
brown.DFM_trans = DFM_trans;
brown.DFM_rot   = DFM_rot;

end
