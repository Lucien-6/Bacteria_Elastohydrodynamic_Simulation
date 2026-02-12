function [theta, phi, psi] = pili_euler_angles(d1, d2, d3)
%PILI_EULER_ANGLES  Recover Euler angles from director vectors
%
%  Description:
%    Inverse of pili_directors: recovers (theta, phi, psi) from the
%    material frame directors (d1, d2, d3).  Adapted from Walker (2020).
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
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Standardised formatting

theta = acos(d3(:, 3));
phi   = atan2(d3(:, 2), d3(:, 1));
psi   = atan2(d1(:, 3), d2(:, 3));

end
