function R = quat_to_rotmat(q)
%QUAT_TO_ROTMAT  Convert unit quaternion to 3x3 rotation matrix
%
%  Description:
%    Converts a unit quaternion q = [q0; q1; q2; q3] (scalar-first) to
%    its equivalent 3x3 rotation matrix R.  Active rotation convention:
%    v_lab = R * v_body  (maps body-frame to lab-frame vectors).
%
%    Explicit matrix:
%      R = [ 1-2(q2^2+q3^2)    2(q1*q2-q0*q3)    2(q1*q3+q0*q2) ]
%          [ 2(q1*q2+q0*q3)    1-2(q1^2+q3^2)    2(q2*q3-q0*q1) ]
%          [ 2(q1*q3-q0*q2)    2(q2*q3+q0*q1)    1-2(q1^2+q2^2) ]
%
%  Syntax:
%    R = quat_to_rotmat(q)
%
%  Inputs:
%    q - [4 x 1] Unit quaternion [q0; q1; q2; q3] (scalar-first)
%
%  Outputs:
%    R - [3 x 3] Rotation matrix such that v_lab = R * v_body
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06

q0 = q(1);  q1 = q(2);  q2 = q(3);  q3 = q(4);

R = [1-2*(q2^2+q3^2),  2*(q1*q2-q0*q3),  2*(q1*q3+q0*q2);
     2*(q1*q2+q0*q3),  1-2*(q1^2+q3^2),  2*(q2*q3-q0*q1);
     2*(q1*q3-q0*q2),  2*(q2*q3+q0*q1),  1-2*(q1^2+q2^2)];

end
