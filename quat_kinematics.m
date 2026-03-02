function dqdt = quat_kinematics(q, Omega)
%QUAT_KINEMATICS  Quaternion time derivative from angular velocity
%
%  Description:
%    Computes dq/dt for a rigid body rotating with angular velocity Omega
%    using the quaternion kinematic equation:
%
%        dq/dt = 0.5 * [0; Omega] (x) q
%
%    where (x) denotes the Hamilton left-product (scalar-first convention).
%    The angular velocity Omega is expressed in the lab frame.
%    The result is used directly by ode15s to advance the body orientation.
%
%  Syntax:
%    dqdt = quat_kinematics(q, Omega)
%
%  Inputs:
%    q     - [4 x 1] Current unit quaternion [q0; q1; q2; q3]
%    Omega - [3 x 1] Angular velocity in lab frame (dimensionless)
%
%  Outputs:
%    dqdt - [4 x 1] Time derivative  dq/dt = 0.5 * [0; Omega] (x) q
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Standardised formatting

omega_quat = [0; Omega(:)];
dqdt = 0.5 * quat_multiply(omega_quat, q);

end
