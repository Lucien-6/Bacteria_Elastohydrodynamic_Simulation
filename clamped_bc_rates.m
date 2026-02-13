function [dtheta1, dphi1, dpsi1] = clamped_bc_rates(Omega_lab, rot_k, theta_1, phi_1)
%CLAMPED_BC_RATES  First-segment Euler angle rates from body rotation
%
%  Description:
%    For a pilus clamped to the body, the first segment's Euler angle
%    rates are determined by the body angular velocity Omega.
%
%  Syntax:
%    [dtheta1, dphi1, dpsi1] = clamped_bc_rates(Omega_lab, rot_k, theta_1, phi_1)
%
%  Inputs:
%    Omega_lab - [3 x 1] Body angular velocity in lab frame (dimensionless)
%    rot_k     - [3 x 3] Adaptive basis rotation (lab -> basis)
%    theta_1   - First segment polar angle (in rot_k basis)
%    phi_1     - First segment azimuthal angle (in rot_k basis)
%
%  Outputs:
%    dtheta1, dphi1, dpsi1 - Euler angle rates for segment 1
%
%  Note:
%    Requires sin(theta_1) != 0 (ensured by adaptive basis rotation).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-13 -- Added sin(theta) singularity guard

%% Transform angular velocity to pili basis
Omega_basis = rot_k * Omega_lab(:);
ox = Omega_basis(1);
oy = Omega_basis(2);
oz = Omega_basis(3);

%% Trigonometric values
sT = sin(theta_1);  cT = cos(theta_1);
sP = sin(phi_1);    cP = cos(phi_1);

%% Guard against singularity (defensive -- should be prevented by adaptive basis)
if abs(sT) < 1e-12
    warning('clamped_bc_rates:SingularTheta', ...
        'sin(theta_1) = %.2e is near zero -- Euler angle singularity. Check adaptive_basis_check.', sT);
    sT = sign(sT + 1e-20) * 1e-12;    % Prevent division by zero
end

%% Euler angle rates (inverted ZYZ kinematics)
ox_cp_oy_sp = ox * cP + oy * sP;

dtheta1 = -ox * sP + oy * cP;
dphi1   =  oz - (cT / sT) * ox_cp_oy_sp;
dpsi1   =  ox_cp_oy_sp / sT;

end
