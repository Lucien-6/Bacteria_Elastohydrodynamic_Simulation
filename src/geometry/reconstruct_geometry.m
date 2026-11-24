function geom = reconstruct_geometry(Z, params)
%RECONSTRUCT_GEOMETRY  Rebuild all spatial positions from state vector
%
%  Description:
%    Given the global state vector Z and parameters, reconstructs the
%    lab-frame positions of all body surface points and pili nodes, as
%    well as the pili director vectors.  Called at the start of each
%    ODE RHS evaluation.
%    All positions are nondimensional coordinates x/L_0; no conversion
%    back to SI units is performed here.  Multiply any positional field by
%    params.scales.L_0 to recover metres, for example
%      X_body_SI = geom.X_body * params.scales.L_0.
%
%  Syntax:
%    geom = reconstruct_geometry(Z, params)
%
%  Inputs:
%    Z      - [7 + 3*sum(params.N_segments) x 1] Global state vector;
%             Z(1:3) is the dimensionless lab-frame body COM position,
%             followed by the quaternion and per-segment Euler angles
%    params - Parameter structure from init_coupled_params; body_local,
%             r_local, and ell are nondimensionalised by L_0
%
%  Outputs:
%    geom - Structure with fields:
%      .X_cm       [3 x 1]     Body COM position in lab frame, x/L_0
%      .q          [4 x 1]     Body quaternion
%      .R_body     [3 x 3]     Rotation from body-local to lab frame
%      .X_body     [N_body x 3] Dimensionless body positions in lab frame
%      .X_pili     {M x 1}     Cell k: [(N_k+1) x 3] dimensionless pili
%                              node positions in lab frame
%      .d1,.d2,.d3 {M x 1}     Cell k: [N_k x 3] directors in lab frame
%      .d1_basis,.d2_basis,.d3_basis
%                    {M x 1}   Cell k: [N_k x 3] directors in rot_k basis
%      .X_all      [N_body+sum(N_k+1) x 3] Concatenated lab-frame positions
%      .pili_angles {M x 1}    Cell k: struct of [N_k x 1] angle arrays
%
%  Warnings:
%    Emits reconstruct_geometry:ZeroQuaternion and substitutes the identity
%    quaternion when the supplied body quaternion has zero norm.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-07-12
%  Updated: 2026-11-04 — v8.7.2: comment and message-text corrections only; no behavioural change.

%% Unpack state
[X_cm, q, pili_angles] = unpack_state(Z, params.N_segments);

% Ensure unit quaternion.  ode15s internal trial steps may produce states
% with |q| != 1; the rotation-matrix formula in quat_to_rotmat is only
% orthonormal for |q| = 1, so a non-unit quaternion would deform the body
% (non-rigid surface points, non-conserved volume).  Normalising here
% protects ALL downstream geometry without disturbing the solver's state
% vector or error control (we only touch the local copy used here).
q_norm = norm(q);
if q_norm > 0
    q = q / q_norm;
else
    % Degenerate zero quaternion — fall back to identity orientation.
    % This branch should never trigger in normal operation; the warning
    % makes any occurrence visible for diagnosis.
    warning('reconstruct_geometry:ZeroQuaternion', ...
        'Body quaternion has zero norm; falling back to identity orientation.');
    q = [1; 0; 0; 0];
end

R_body = quat_to_rotmat(q);

%% Body surface points in lab frame
X_body_lab = (R_body * params.body_local)' + X_cm(:)';   % [N_body x 3]

%% Pili geometry
M = params.M_pili;
X_pili    = cell(M, 1);
d1_lab    = cell(M, 1);   d2_lab    = cell(M, 1);   d3_lab    = cell(M, 1);
d1_basis  = cell(M, 1);   d2_basis  = cell(M, 1);   d3_basis  = cell(M, 1);
d1N_basis = cell(M, 1);   d2N_basis = cell(M, 1);   d3N_basis = cell(M, 1);

for k = 1:M
    N_k   = params.N_segments(k);
    rot_k = params.rot{k};
    ell_k = params.ell(k);

    theta_k = pili_angles{k}.theta;
    phi_k   = pili_angles{k}.phi;
    psi_k   = pili_angles{k}.psi;

    % Root position in lab frame
    x_root = X_cm + R_body * params.r_local{k};

    % Directors in rot_k basis (computed once and reused for positions + lab)
    [d1b, d2b, d3b] = pili_directors(theta_k, phi_k, psi_k);
    d1_basis{k} = d1b;
    d2_basis{k} = d2b;
    d3_basis{k} = d3b;

    % Directors in lab frame
    d1_lab{k} = d1b * rot_k;
    d2_lab{k} = d2b * rot_k;
    d3_lab{k} = d3b * rot_k;

    % Node positions directly from d3_lab (avoids duplicate trig in pili_spatial_coords)
    ds_k = ell_k / N_k;
    X_pili{k} = [x_root(:)'; x_root(:)' + ds_k * cumsum(d3_lab{k}, 1)];

    % Nodal directors in rot_k basis (cached once per RHS for B_global, elastic & active moments)
    [d1Nb, d2Nb, d3Nb] = compute_node_directors(d1b, d2b, d3b, N_k);
    d1N_basis{k} = d1Nb;
    d2N_basis{k} = d2Nb;
    d3N_basis{k} = d3Nb;
end

%% Concatenate all points
N_body_pts = size(X_body_lab, 1);
N_pili_pts = sum(params.N_segments(:)) + M;
X_all = zeros(N_body_pts + N_pili_pts, 3);
X_all(1:N_body_pts, :) = X_body_lab;

pt_idx = N_body_pts + 1;
for k = 1:M
    n_pts = params.N_segments(k) + 1;
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
geom.d1N_basis   = d1N_basis;
geom.d2N_basis   = d2N_basis;
geom.d3N_basis   = d3N_basis;
geom.X_all       = X_all;
geom.pili_angles = pili_angles;

end
