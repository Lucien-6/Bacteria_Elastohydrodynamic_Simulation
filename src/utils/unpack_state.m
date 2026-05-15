function [X_cm, q, pili_angles] = unpack_state(Z, N_segments)
%UNPACK_STATE  Extract body state and pili Euler angles from state vector
%
%  Description:
%    Extracts the body centre-of-mass position, orientation quaternion,
%    and per-pili Euler angle arrays from the flat global state vector Z.
%    Inverse operation of pack_state_rates.
%
%    Expected Z layout:
%      Z(1:3)   = X_cm  (body position, dimensionless)
%      Z(4:7)   = q     (body quaternion)
%      Z(8:end) = per-pili Euler angles [theta_k; phi_k; psi_k], each Nk
%
%  Syntax:
%    [X_cm, q, pili_angles] = unpack_state(Z, N_segments)
%
%  Inputs:
%    Z          - [n x 1] Global state vector
%    N_segments - [M x 1] Number of segments per pilus
%
%  Outputs:
%    X_cm        - [3 x 1] Body centre-of-mass position
%    q           - [4 x 1] Body orientation quaternion
%    pili_angles - {M x 1} Cell array of structs (.theta, .phi, .psi)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06

X_cm = Z(1:3);
q    = Z(4:7);

M           = length(N_segments);
pili_angles = cell(M, 1);
idx         = 8;

for k = 1:M
    Nk = N_segments(k);
    pili_angles{k}.theta = Z(idx       : idx+Nk-1);
    pili_angles{k}.phi   = Z(idx+Nk    : idx+2*Nk-1);
    pili_angles{k}.psi   = Z(idx+2*Nk  : idx+3*Nk-1);
    idx = idx + 3*Nk;
end

end
