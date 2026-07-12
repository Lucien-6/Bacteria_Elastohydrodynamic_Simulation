function R = quat_to_rotmat(q)
%QUAT_TO_ROTMAT  Convert a scalar-first quaternion to a rotation matrix
%
%  Description:
%    Converts q = [q0; q1; q2; q3] to the active rotation matrix satisfying
%      v_lab = R*v_body.
%    The input is reshaped to a column and normalised internally, so any
%    finite nonzero scalar multiple of a rotation quaternion represents the
%    same rotation.  The exact zero quaternion has no defined orientation
%    and is treated as the identity quaternion [1; 0; 0; 0].
%
%    No near-zero tolerance is applied: the implementation tests norm(q) > 0.
%    Callers should therefore supply a finite, meaningfully nonzero
%    quaternion; the fallback is only a defensive guard for the exact-zero
%    case.  Internal normalisation makes the returned matrix orthonormal up
%    to floating-point round-off.
%
%  Syntax:
%    R = quat_to_rotmat(q)
%
%  Inputs:
%    q - Four-element scalar-first quaternion [q0; q1; q2; q3]
%
%  Outputs:
%    R - [3 x 3] Active body-to-lab rotation matrix
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-07-12

%  Defensive: ensure unit quaternion.  Although reconstruct_geometry
%  normalises upstream, this guards any other callers (e.g. test scripts,
%  post-processing, visualisation helpers) against |q| != 1, which would
%  otherwise produce a non-orthonormal rotation matrix.  Cost: one 4-vector
%  norm + division per call — negligible.
q = q(:);
q_norm = norm(q);
if q_norm > 0
    q = q / q_norm;
else
    q = [1; 0; 0; 0];  % identity fallback
end

q0 = q(1);  q1 = q(2);  q2 = q(3);  q3 = q(4);

R = [1-2*(q2^2+q3^2),  2*(q1*q2-q0*q3),  2*(q1*q3+q0*q2);
     2*(q1*q2+q0*q3),  1-2*(q1^2+q3^2),  2*(q2*q3-q0*q1);
     2*(q1*q3-q0*q2),  2*(q2*q3+q0*q1),  1-2*(q1^2+q2^2)];

end
