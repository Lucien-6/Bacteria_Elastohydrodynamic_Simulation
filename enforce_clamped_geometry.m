function Z = enforce_clamped_geometry(Z, params)
%ENFORCE_CLAMPED_GEOMETRY  Project first-segment Euler angles to match body
%
%  Description:
%    Corrects constraint drift at the pili-body junction by enforcing:
%
%        D_lab^{(1)} = R_body * D_local{k}
%
%    The ODE enforces clamped BCs at the rate (kinematic) level only.
%    This function acts as a geometric projection step (analogous to
%    quaternion renormalisation) to keep the constraint satisfied at
%    the position level.
%
%    Call after each ODE integration interval and after quaternion
%    renormalisation.
%
%  Syntax:
%    Z = enforce_clamped_geometry(Z, params)
%
%  Inputs:
%    Z      - [n x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    Z - [n x 1] State vector with corrected first-segment angles
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-07
%  Updated: 2026-02-11 — Standardised formatting

M = params.M_pili;
if M == 0
    return;
end

%% Unpack body orientation
q      = Z(4:7);
R_body = quat_to_rotmat(q);

%% Correct each pilus's first-segment Euler angles
N_segs = params.N_segments;
idx    = 8;

for k = 1:M
    N_k   = N_segs(k);
    rot_k = params.rot{k};

    % Expected director frame in lab frame
    D_lab   = R_body * params.D_local{k};
    D_basis = rot_k * D_lab;

    d1_basis = D_basis(:, 1)';
    d2_basis = D_basis(:, 2)';
    d3_basis = D_basis(:, 3)';

    [theta_c, phi_c, psi_c] = pili_euler_angles(d1_basis, d2_basis, d3_basis);

    % Overwrite first-segment angles
    Z(idx)          = theta_c;
    Z(idx + N_k)    = phi_c;
    Z(idx + 2*N_k)  = psi_c;

    idx = idx + 3 * N_k;
end

end
