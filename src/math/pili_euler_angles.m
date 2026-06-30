function [theta, phi, psi] = pili_euler_angles(d1, d2, d3)
%PILI_EULER_ANGLES  Recover Euler angles from director vectors
%
%  Description:
%    Inverse of pili_directors: recovers (theta, phi, psi) from the
%    material frame directors (d1, d2, d3).  Adapted from Walker (2020).
%
%  Euler angle convention (Walker 2020, second-z rotation convention):
%    The tangent director is the spherical-coordinate parameterisation
%      d3 = (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
%    so theta is the polar angle from the z-axis and phi is the azimuth.
%    The material frame (d1, d2) is obtained from d3 by a rotation of
%    angle psi about d3:
%      d1 = (-sin(phi)*cos(psi) - cos(theta)*cos(phi)*sin(psi),
%             cos(phi)*cos(psi) - cos(theta)*sin(phi)*sin(psi),
%             sin(theta)*sin(psi))
%      d2 = ( sin(phi)*sin(psi) - cos(theta)*cos(phi)*cos(psi),
%            -cos(phi)*sin(psi) - cos(theta)*sin(phi)*cos(psi),
%             sin(theta)*cos(psi))
%    This is the convention used by Walker (2020) directors.m.  The
%    rotation sequence is R = Rz(phi) * Ry(-theta) * Rz(psi) (equivalently
%    the transpose of the standard ZYZ active rotation), so it is NOT the
%    textbook ZXZ convention despite the similar (theta, phi, psi) naming.
%    The clamped-BC Jacobian J_euler_basis in build_coupled_lin_sys.m is
%    derived consistently from this convention.
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
