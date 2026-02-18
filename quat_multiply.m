function r = quat_multiply(p, q)
%QUAT_MULTIPLY  Hamilton product of two quaternions
%
%  Description:
%    Computes the Hamilton (left) product r = p (x) q using the
%    scalar-first convention q = [q0; qx; qy; qz].
%
%    Formula:
%        r0 = p0*q0 - pv . qv
%        rv = p0*qv + q0*pv + pv x qv
%
%    Rotation composition convention: p (x) q represents the rotation
%    obtained by first applying q, then p (right-to-left composition).
%
%  Syntax:
%    r = quat_multiply(p, q)
%
%  Inputs:
%    p, q - [4 x 1] Quaternions [scalar; vector] (scalar-first)
%
%  Outputs:
%    r - [4 x 1] Product quaternion p (x) q
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Standardised formatting

p0 = p(1);  pv = p(2:4);
q0 = q(1);  qv = q(2:4);

r = [p0*q0 - dot(pv, qv);
     p0*qv + q0*pv + cross(pv, qv)];

end
