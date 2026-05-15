function [theta, phi, psi] = pili_euler_angles(d1, d2, d3)
%PILI_EULER_ANGLES  Recover Euler angles from director vectors
%
%  Description:
%    Inverse of pili_directors: recovers (theta, phi, psi) from the
%    material frame directors (d1, d2, d3).  Adapted from Walker (2020).
%
%  Euler angle convention (ZXZ intrinsic rotation):
%    The directors are obtained by applying three successive rotations
%    starting from the reference basis {e1, e2, e3}:
%      1. R_z(phi):   rotate about z by phi      (azimuthal)
%      2. R_x(theta): rotate about new x by theta (polar / tilt)
%      3. R_z(psi):   rotate about new z by psi   (material twist)
%    This is the ZXZ convention (same as Walker 2020).  The tangent
%    d3 = R_z(phi) * R_x(theta) * e3, and the material frame (d1, d2)
%    are further rotated by psi around d3.
%
%  Syntax:
%    [theta, phi, psi] = pili_euler_angles(d1, d2, d3)
%
%  Inputs:
%    d1 - [N x 3] First director
%    d2 - [N x 3] Second director
%    d3 - [N x 3] Third director (tangent)
%
%  Outputs:
%    theta - [N x 1] Polar angles   (0 to pi)
%    phi   - [N x 1] Azimuthal angles (-pi to pi)
%    psi   - [N x 1] Twist angles    (-pi to pi)
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%         (Euler angle parametrisation; directors.m inverse)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06

cos_theta = max(-1, min(1, d3(:, 3)));
theta = acos(cos_theta);
phi   = atan2(d3(:, 2), d3(:, 1));
psi   = atan2(d1(:, 3), d2(:, 3));

end
