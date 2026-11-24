function [force_over_mu_U, torque_over_mu_U] = slender_body_load(X, U, a_cyl, origin)
%SLENDER_BODY_LOAD  Johnson slender-body force and torque on a rigid centreline
%
%  Description:
%    Reference load for the helical-filament check.  The discrete operator is
%    Johnson's slender-body equation as used by Rodenborn et al., Proc. Natl.
%    Acad. Sci. U.S.A. 110, E338 (2013): a logarithmically singular local term
%    for a prolate cross-section, plus the finite-part Stokeslet.  The
%    separate a^2 dipole is left in the prolate logarithm: a raw neighbour
%    dipole without its analytic self-integral does not converge with the
%    segment length.
%
%    The prolate section vanishes at the ends.  Its maximum radius is
%    4*a_cyl/pi, so the arc-length average equals the cylinder radius a_cyl
%    used by the regularised segment line.  On that section the local
%    logarithm is constant,
%
%        Lambda = 2 * log(L_contour / a_max).
%
%    Unknowns are the force density at segment midpoints.  Outputs use the
%    same normalisation as rss_line_load.
%
%  Syntax:
%    [force_over_mu_U, torque_over_mu_U] = slender_body_load(X, U, a_cyl)
%    [force_over_mu_U, torque_over_mu_U] = slender_body_load(X, U, a_cyl, origin)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-15

if nargin < 4 || isempty(origin)
    origin = zeros(1, 3);
end
X = X - origin;
dX = diff(X, 1, 1);
ds = sqrt(sum(dX.^2, 2));
if any(ds <= 0)
    error('slender_body_load:ZeroSegment', 'A centreline segment has zero length.');
end
t = dX ./ ds;
mid = 0.5 * (X(1:end-1, :) + X(2:end, :));
Um = 0.5 * (U(1:end-1, :) + U(2:end, :));
n_seg = size(mid, 1);
s = [0; cumsum(ds)];
s_mid = 0.5 * (s(1:end-1) + s(2:end));
L_c = s(end);
a_max = 4 * a_cyl / pi;
lambda = 2 * log(L_c / a_max);

A = zeros(3 * n_seg, 3 * n_seg);
jx = 1:3:3 * n_seg;
jy = 2:3:3 * n_seg;
jz = 3:3:3 * n_seg;
for i = 1:n_seg
    ti = t(i, :).';
    tt = ti * ti.';
    ii = (3 * (i - 1) + 1):(3 * i);
    local = (eye(3) + tt) * lambda + (eye(3) - 3 * tt);

    R = mid(i, :) - mid;
    r2 = sum(R.^2, 2);
    r2(i) = inf;
    inv_r = 1 ./ sqrt(r2);
    dsj = ds;
    dsj(i) = 0;
    %  Leading slender-body operator: finite-part Stokeslet.  The a^2 dipole
    %  is already represented by the prolate logarithm; adding the raw
    %  neighbour dipole without its analytic self-integral moves the thrust
    %  by several percent when the segment length changes.
    iso = inv_r .* dsj;
    rr_c = (inv_r.^3) .* dsj;
    arc = abs(s_mid(i) - s_mid);
    arc(i) = inf;
    local = local - (eye(3) + tt) * sum(dsj ./ arc);

    Rx = R(:, 1); Ry = R(:, 2); Rz = R(:, 3);
    A(ii(1), jx) = (iso + rr_c .* Rx .* Rx).';
    A(ii(1), jy) = (rr_c .* Rx .* Ry).';
    A(ii(1), jz) = (rr_c .* Rx .* Rz).';
    A(ii(2), jx) = (rr_c .* Ry .* Rx).';
    A(ii(2), jy) = (iso + rr_c .* Ry .* Ry).';
    A(ii(2), jz) = (rr_c .* Ry .* Rz).';
    A(ii(3), jx) = (rr_c .* Rz .* Rx).';
    A(ii(3), jy) = (rr_c .* Rz .* Ry).';
    A(ii(3), jz) = (iso + rr_c .* Rz .* Rz).';
    A(ii, ii) = local;
end

f = A \ (8 * pi * reshape(Um.', [], 1));
f = reshape(f, 3, []).';
force_over_mu_U = (ds.' * f).';
moment = cross(mid, f, 2);
torque_over_mu_U = (ds.' * moment).';
end
