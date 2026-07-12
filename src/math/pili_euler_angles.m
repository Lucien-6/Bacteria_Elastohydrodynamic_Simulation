function [theta, phi, psi] = pili_euler_angles(d1, d2, d3)
%PILI_EULER_ANGLES  Recover Walker Euler angles from material directors
%
%  Description:
%    Inverts pili_directors for a right-handed orthonormal director frame.
%    All three input arrays must be expressed in the same coordinate basis;
%    in normal use this is the pilus rotated basis.  Returned angles are in
%    radians.
%
%    The Walker (2020) convention is
%      d3 = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
%    with frame rotation R = Rz(phi)*Ry(-theta)*Rz(psi).  Thus theta is
%    the polar angle, phi is the tangent azimuth, and psi is material twist
%    about d3.  This is not the textbook ZXZ convention.
%
%    d3(:,3) is clamped to [-1, 1] before acos to suppress round-off-domain
%    errors.  The routine otherwise assumes valid directors and does not
%    orthonormalise them.  At theta = 0 or pi, phi and psi are not
%    individually unique because of gimbal lock; atan2 returns one numerical
%    representative.  Callers requiring nonsingular coordinates should use
%    adaptive_basis_check before approaching either pole.
%
%  Syntax:
%    [theta, phi, psi] = pili_euler_angles(d1, d2, d3)
%
%  Inputs:
%    d1 - [N x 3] First material director
%    d2 - [N x 3] Second material director
%    d3 - [N x 3] Tangent director
%
%  Outputs:
%    theta - [N x 1] Polar angles in [0, pi], in radians
%    phi   - [N x 1] Azimuthal angles in [-pi, pi], in radians
%    psi   - [N x 1] Twist angles in [-pi, pi], in radians
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%        (Euler angle parametrisation; directors.m inverse)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-07-12

cos_theta = max(-1, min(1, d3(:, 3)));
theta = acos(cos_theta);
phi   = atan2(d3(:, 2), d3(:, 1));
psi   = atan2(d1(:, 3), d2(:, 3));

end
