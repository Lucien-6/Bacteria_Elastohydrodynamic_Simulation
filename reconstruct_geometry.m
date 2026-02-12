function geom = reconstruct_geometry(Z, params)
%RECONSTRUCT_GEOMETRY  Rebuild all spatial positions from state vector
%
%  Description:
%    Given the global state vector Z and parameters, reconstructs the
%    lab-frame positions of all body surface points and pili nodes, as
%    well as the pili director vectors.  Called at the start of each
%    ODE RHS evaluation.
%
%  Syntax:
%    geom = reconstruct_geometry(Z, params)
%
%  Inputs:
%    Z      - [n x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    geom - Structure with fields:
%      .X_cm       [3 x 1]     Body centre of mass
%      .q          [4 x 1]     Body quaternion
%      .R_body     [3 x 3]     Body rotation matrix
%      .X_body     [N_body x 3] Body surface positions (lab)
%      .X_pili     {M x 1}     Pili node positions (lab)
%      .d1,.d2,.d3 {M x 1}     Directors in lab frame
%      .d1_basis,.d2_basis,.d3_basis {M x 1} Directors in rot_k basis
%      .X_all      [N_total x 3] All points concatenated
%      .pili_angles {M x 1}    Angle structs
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Standardised formatting

%% Unpack state
[X_cm, q, pili_angles] = unpack_state(Z, params.N_segments);
R_body = quat_to_rotmat(q);

%% Body surface points in lab frame
X_body_lab = (R_body * params.body_local)' + X_cm(:)';   % [N_body x 3]

%% Pili geometry
M = params.M_pili;
X_pili   = cell(M, 1);
d1_lab   = cell(M, 1);   d2_lab   = cell(M, 1);   d3_lab   = cell(M, 1);
d1_basis = cell(M, 1);   d2_basis = cell(M, 1);   d3_basis = cell(M, 1);

for k = 1:M
    N_k   = params.N_segments(k);
    rot_k = params.rot{k};
    ell_k = params.ell(k);

    theta_k = pili_angles{k}.theta;
    phi_k   = pili_angles{k}.phi;
    psi_k   = pili_angles{k}.psi;

    % Root position in lab frame
    x_root = X_cm + R_body * params.r_local{k};

    % Node positions
    X_k        = pili_spatial_coords(x_root, theta_k, phi_k, rot_k, ell_k, N_k);
    X_pili{k}  = X_k;

    % Directors in rot_k basis
    [d1b, d2b, d3b] = pili_directors(theta_k, phi_k, psi_k);
    d1_basis{k} = d1b;
    d2_basis{k} = d2b;
    d3_basis{k} = d3b;

    % Directors in lab frame
    d1_lab{k} = (rot_k' * d1b')';
    d2_lab{k} = (rot_k' * d2b')';
    d3_lab{k} = (rot_k' * d3b')';
end

%% Concatenate all points
N_body_pts = size(X_body_lab, 1);
N_pili_pts = sum(cellfun(@(x) size(x, 1), X_pili));
X_all = zeros(N_body_pts + N_pili_pts, 3);
X_all(1:N_body_pts, :) = X_body_lab;

pt_idx = N_body_pts + 1;
for k = 1:M
    n_pts = size(X_pili{k}, 1);
    X_all(pt_idx : pt_idx + n_pts - 1, :) = X_pili{k};
    pt_idx = pt_idx + n_pts;
end

%% Pack output
geom.X_cm        = X_cm;
geom.q           = q;
geom.R_body      = R_body;
geom.X_body      = X_body_lab;
geom.X_pili      = X_pili;
geom.d1          = d1_lab;
geom.d2          = d2_lab;
geom.d3          = d3_lab;
geom.d1_basis    = d1_basis;
geom.d2_basis    = d2_basis;
geom.d3_basis    = d3_basis;
geom.X_all       = X_all;
geom.pili_angles = pili_angles;

end
