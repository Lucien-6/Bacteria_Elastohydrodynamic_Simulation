function r = quat_multiply(p, q)
%QUAT_MULTIPLY  Compute the scalar-first Hamilton product
%
%  Description:
%    Computes r = p (x) q using scalar-first quaternions
%      p = [p0; pv], q = [q0; qv],
%      r0 = p0*q0 - dot(pv,qv),
%      rv = p0*qv + q0*pv + cross(pv,qv).
%    For rotation quaternions, p (x) q applies q first and p second.
%
%    This routine performs no normalisation, norm tolerance, or zero-norm
%    fallback.  The product is unit length only when both inputs are unit
%    quaternions, up to floating-point round-off.  Callers that store or
%    repeatedly compose rotations are responsible for renormalisation.
%    Non-unit inputs are accepted as algebraic quaternions.
%
%  Syntax:
%    r = quat_multiply(p, q)
%
%  Inputs:
%    p - [4 x 1] Left scalar-first quaternion [p0; px; py; pz]
%    q - [4 x 1] Right scalar-first quaternion [q0; qx; qy; qz]
%
%  Outputs:
%    r - [4 x 1] Unnormalised Hamilton product p (x) q
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-07-12

p0 = p(1);  pv = p(2:4);
q0 = q(1);  qv = q(2:4);

r = [p0*q0 - dot(pv, qv);
     p0*qv + q0*pv + cross(pv, qv)];

end
