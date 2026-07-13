function tau = solve_motor_self_consistent(a, B, tau_motor, om_max, logf_h)
%SOLVE_MOTOR_SELF_CONSISTENT  Solve reduced motor stall-curve coupling
%
%  Solves omega = a + B*tau together with the elementwise stall law.  A
%  scalar motor uses the exact piecewise-linear root.  Multiple motors use
%  damped fixed-point iteration and must satisfy the true fixed-point
%  residual before returning.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13

if nargin < 5
    logf_h = [];
end

a = a(:);
tau_motor = tau_motor(:);
om_max = om_max(:);
n = numel(a);

if ~isequal(size(B), [n, n]) || numel(tau_motor) ~= n || numel(om_max) ~= n
    error('solve_motor_self_consistent:DimensionMismatch', ...
        'Motor vectors and sensitivity matrix dimensions are inconsistent.');
end
if any(~isfinite(a)) || any(~isfinite(B(:))) || ...
        any(~isfinite(tau_motor)) || any(~isfinite(om_max))
    error('solve_motor_self_consistent:NonFiniteInput', ...
        'Motor self-consistency inputs must be finite.');
end

if n == 1
    tau = motor_law_root_scalar(a, B(1,1), tau_motor(1), om_max(1));
    return;
end

tau = zeros(n, 1);
relax = 0.7;
tolerance = 1e-12;
for iteration = 1:500
    omega = a + B * tau;
    tau_target = motor_law_vec(omega, tau_motor, om_max);
    tau = (1 - relax) * tau + relax * tau_target;

    fixed_point_residual = tau - motor_law_vec( ...
        a + B * tau, tau_motor, om_max);
    if norm(fixed_point_residual, inf) <= ...
            tolerance * (1 + norm(tau, inf))
        return;
    end
end

fixed_point_residual = tau - motor_law_vec( ...
    a + B * tau, tau_motor, om_max);
relative_residual = norm(fixed_point_residual, inf) / ...
    (1 + norm(tau, inf));
if ~isempty(logf_h)
    logf_h(['[ODE] ERROR: motor self-consistent solve failed after 500 ', ...
            'iterations (n=%d, residual=%.3e).\n'], n, relative_residual);
end
error('solve_motor_self_consistent:NoConvergence', ...
    ['Motor self-consistent solve failed after 500 iterations ', ...
     '(n=%d, relative residual=%.3e).'], n, relative_residual);
end


function tau = motor_law_root_scalar(a, b, tau_motor, om_max)
%MOTOR_LAW_ROOT_SCALAR  Exact self-consistent torque for one motor
if ~(om_max > 0)
    tau = tau_motor;
    return;
end

omega_positive = (a + b*tau_motor) / ...
    (1 + b*tau_motor/om_max);
if omega_positive >= 0 && abs(omega_positive) < om_max
    tau = tau_motor * (1 - abs(omega_positive)/om_max);
    return;
end

omega_negative = (a + b*tau_motor) / ...
    (1 - b*tau_motor/om_max);
if omega_negative < 0 && abs(omega_negative) < om_max
    tau = tau_motor * (1 - abs(omega_negative)/om_max);
    return;
end

tau = 0;
end


function tau = motor_law_vec(omega, tau_motor, om_max)
%MOTOR_LAW_VEC  Evaluate the elementwise motor stall curve
n = numel(omega);
tau = zeros(n, 1);
for i = 1:n
    if om_max(i) > 0 && abs(omega(i)) < om_max(i)
        tau(i) = tau_motor(i) * (1 - abs(omega(i)) / om_max(i));
    else
        tau(i) = 0;
    end
end
end
