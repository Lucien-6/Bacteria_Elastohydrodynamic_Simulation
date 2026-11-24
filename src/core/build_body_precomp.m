function pre = build_body_precomp(body_local, eps_body)
%BUILD_BODY_PRECOMP  Precompute body-frame rigid-body RSM mobility factors
%
%  Description:
%    The cell body is rigid, so in the body-local frame its self-interaction
%    regularised Stokeslet block M_BB_loc (3*N_body x 3*N_body) is symmetric
%    positive definite and constant throughout the entire simulation.
%    Precomputing its Cholesky factorisation and its action on the 6 rigid-
%    body kinematic columns Q_B6_loc once during initialisation allows
%    build_coupled_lin_sys to eliminate all 3*N_body body unknowns via a
%    block Schur complement in the body frame, reducing the per-step LU
%    factorisation from O((3*N_body + 3*N_pili)^3) to O((3*N_pili)^3).
%
%    For bare-body simulations (M_pili == 0), the body-frame 6x6 resistance
%    matrix is also precomputed here, turning each RHS evaluation into a
%    single 6x6 rotation.
%
%  Syntax:
%    pre = build_body_precomp(body_local, eps_body)
%
%  Inputs:
%    body_local - [3 x N_body] Dimensionless body surface points in body frame
%    eps_body   - Positive scalar RSM regularisation length (dimensionless)
%
%  Outputs:
%    pre - Struct with fields (or empty [] if Cholesky fails):
%            .N_body, .nB, .eps_body
%            .dM_BB           Cholesky decomposition object of M_BB_loc
%            .U_QB_loc        [3*N_body x 6] = M_BB_loc \ Q_B6_loc
%            .lin_sys_bare_loc [6 x 6] = -B_B6_loc * U_QB_loc (body frame)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

X_b    = body_local';                     % [N_body x 3]
N_body = size(X_b, 1);
nB     = 3 * N_body;
e2     = eps_body^2;

% Assemble symmetric M_BB_loc in body frame
rx = X_b(:, 1) - X_b(:, 1)';
ry = X_b(:, 2) - X_b(:, 2)';
rz = X_b(:, 3) - X_b(:, 3)';
r2 = rx.^2 + ry.^2 + rz.^2;

inv_R = 1 ./ sqrt(r2 + e2);
H2    = inv_R .* inv_R .* inv_R;
H1    = (r2 + 2 * e2) .* H2;

M_BB = zeros(nB, nB);
ix = 1:3:nB;  iy = 2:3:nB;  iz = 3:3:nB;
G_xy = H2 .* (rx .* ry);
G_xz = H2 .* (rx .* rz);
G_yz = H2 .* (ry .* rz);

M_BB(ix, ix) = H1 + H2 .* (rx .* rx);
M_BB(ix, iy) = G_xy;
M_BB(ix, iz) = G_xz;
M_BB(iy, ix) = G_xy;
M_BB(iy, iy) = H1 + H2 .* (ry .* ry);
M_BB(iy, iz) = G_yz;
M_BB(iz, ix) = G_xz;
M_BB(iz, iy) = G_yz;
M_BB(iz, iz) = H1 + H2 .* (rz .* rz);

% Enforce exact symmetry before Cholesky
M_BB = 0.5 * (M_BB + M_BB');

% Verify strict positive definiteness via chol; fall back gracefully if ill-conditioned
[L_BB, p_flag] = chol(M_BB, 'lower');
if p_flag ~= 0 || any(~isfinite(L_BB(:)))
    pre = [];
    return;
end

dM_BB = decomposition(M_BB, 'chol', 'lower', 'CheckCondition', false);

% Rigid-body kinematic columns Q_B6_loc [nB x 6] and force-balance rows B_B6_loc [6 x nB]
% Note: B_B6_loc = Q_B6_loc' (transposed rigid-body velocity-force duality, w = 1 on body)
Q_B6 = zeros(nB, 6);
Q_B6(ix, 1) = 1;           Q_B6(iy, 2) = 1;           Q_B6(iz, 3) = 1;
Q_B6(ix, 5) =  X_b(:, 3);  Q_B6(ix, 6) = -X_b(:, 2);
Q_B6(iy, 4) = -X_b(:, 3);  Q_B6(iy, 6) =  X_b(:, 1);
Q_B6(iz, 4) =  X_b(:, 2);  Q_B6(iz, 5) = -X_b(:, 1);

U_QB_loc = dM_BB \ Q_B6;
lin_sys_bare_loc = -(Q_B6' * U_QB_loc);

pre = struct();
pre.N_body           = N_body;
pre.nB               = nB;
pre.eps_body         = eps_body;
pre.dM_BB            = dM_BB;
pre.U_QB_loc         = U_QB_loc;
pre.lin_sys_bare_loc = lin_sys_bare_loc;

end
