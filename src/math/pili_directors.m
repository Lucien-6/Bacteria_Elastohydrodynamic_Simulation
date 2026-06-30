function [d1, d2, d3] = pili_directors(theta, phi, psi)
%PILI_DIRECTORS  Material frame directors from Euler angles
%
%  Description:
%    Computes material frame directors (d1, d2, d3) for each segment from
%    Euler angles (theta, phi, psi).  d3 is the tangent direction.
%    Adapted from Walker (2020) directors.m.
%
%  Syntax:
%    [d1, d2, d3] = pili_directors(theta, phi, psi)
%
%  Inputs:
%    theta - [N x 1] Polar angles   (from z-axis)
%    phi   - [N x 1] Azimuthal angles (in xy-plane)
%    psi   - [N x 1] Twist angles
%
%  Outputs:
%    d1 - [N x 3] First director  (normal)
%    d2 - [N x 3] Second director (binormal)
%    d3 - [N x 3] Third director  (tangent)
%
%  Note:
%    Directors are in the pili's rotated basis (rot_k).
%    Lab frame: d_lab = (rot_k' * d_basis')'
%
%  Euler angle convention:
%    Uses the Walker (2020) convention (NOT textbook ZXZ).  The rotation
%    sequence is R = Rz(phi) * Ry(-theta) * Rz(psi), i.e. the transpose of
%    the standard ZYZ active rotation.  See pili_euler_angles.m for the full
%    convention specification and J_euler_basis derivation.
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06

THETA = theta(:);   PHI = phi(:);   PSI = psi(:);

cT = cos(THETA);  sT = sin(THETA);
cP = cos(PHI);    sP = sin(PHI);
cS = cos(PSI);    sS = sin(PSI);

d1 = [-sP.*cS - cT.*cP.*sS,   cP.*cS - cT.*sP.*sS,   sT.*sS];
d2 = [ sP.*sS - cT.*cP.*cS,  -cP.*sS - cT.*sP.*cS,   sT.*cS];
d3 = [ sT.*cP,                 sT.*sP,                 cT    ];

end
