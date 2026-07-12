function [d1, d2, d3] = pili_directors(theta, phi, psi)
%PILI_DIRECTORS  Construct material directors from Walker Euler angles
%
%  Description:
%    Computes the right-handed material frame (d1, d2, d3) for every
%    segment, with d3 the tangent.  theta, phi, and psi are in radians and
%    may be row or column vectors; outputs contain one director per row.
%
%    Directors are expressed in the pilus rotated basis.  If params.rot{k}
%    maps lab coordinates to basis coordinates, the corresponding lab-frame
%    directors are obtained as
%      d_lab = (params.rot{k}'*d_basis')'.
%
%    The Walker (2020) convention is
%      R = Rz(phi)*Ry(-theta)*Rz(psi),
%    not textbook ZXZ.  psi rotates d1 and d2 about d3 and therefore does
%    not change the tangent.  At theta = 0 or pi the Euler coordinates are
%    gimbal-locked, although the forward director frame remains defined.
%
%  Syntax:
%    [d1, d2, d3] = pili_directors(theta, phi, psi)
%
%  Inputs:
%    theta - [N x 1] Polar angles from the basis z-axis, in radians
%    phi   - [N x 1] Azimuthal angles in the basis xy-plane, in radians
%    psi   - [N x 1] Material twist angles, in radians
%
%  Outputs:
%    d1 - [N x 3] First material director in the rotated basis
%    d2 - [N x 3] Second material director in the rotated basis
%    d3 - [N x 3] Tangent director in the rotated basis
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-07-12

THETA = theta(:);   PHI = phi(:);   PSI = psi(:);

cT = cos(THETA);  sT = sin(THETA);
cP = cos(PHI);    sP = sin(PHI);
cS = cos(PSI);    sS = sin(PSI);

d1 = [-sP.*cS - cT.*cP.*sS,   cP.*cS - cT.*sP.*sS,   sT.*sS];
d2 = [ sP.*sS - cT.*cP.*cS,  -cP.*sS - cT.*sP.*cS,   sT.*cS];
d3 = [ sT.*cP,                 sT.*sP,                 cT    ];

end
