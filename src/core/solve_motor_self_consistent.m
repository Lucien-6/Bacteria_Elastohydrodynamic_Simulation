function tau = solve_motor_self_consistent(a, B, tau_motor, om_max, logf_h)
%SOLVE_MOTOR_SELF_CONSISTENT  Solve reduced motor torque-speed coupling
%
%  Description:
%    Solves the reduced motor problem
%        omega = a + B * tau                      (linear hydro/elastic response)
%        tau_i = tau_motor_i * (1 - omega_rel,i / om_max_i)   (torque-speed law)
%        omega_rel,i = sign(tau_motor_i) * omega_i
%    for the effective motor torques tau.  omega is the rotor speed of each
%    motorised pilus RELATIVE to the cell body (see ode_coupled_rhs Step 8d);
%    omega_rel is that speed measured positive in the motor's own driving
%    direction, so a single law serves CW and CCW motors alike.
%
%    The law is a signed, asymmetric linear torque-speed relation:
%      omega_rel = 0        -> tau = tau_motor            (stall torque)
%      0 < omega_rel < w_m  -> tau decreases linearly     (working branch)
%      omega_rel = w_m      -> tau = 0                    (zero-torque speed)
%      omega_rel > w_m      -> tau < 0                    (motor acts as a brake)
%      omega_rel < 0        -> |tau| > |tau_motor|        (forced backwards:
%                              torque rises above stall)
%    This follows the measured behaviour of the bacterial flagellar motor,
%    whose torque-speed curve continues linearly through zero torque and
%    rises for negative speeds (Berg & Turner 1993; Berg 2003), replacing the
%    former |omega| law, which was symmetric in the rotation direction and
%    therefore unphysical for reverse rotation.
%
%    Because the law is linear in omega the fixed point is the solution of
%    a small linear system:
%        (I + diag(|tau_motor|/om_max) * B) * tau = tau_motor - diag(|tau_motor|/om_max) * a
%    B (the response of the relative rotor speeds to unit motor torques) is
%    dissipative, so the system matrix is well conditioned and the solution
%    is unique for any number of motors -- no iteration is needed.
%    om_max <= 0 (or non-finite) disables speed limiting for that motor
%    (tau = tau_motor).
%
%  Syntax:
%    tau = solve_motor_self_consistent(a, B, tau_motor, om_max)
%    tau = solve_motor_self_consistent(a, B, tau_motor, om_max, logf_h)
%
%  Inputs:
%    a         - [n x 1] Relative rotor speeds with zero motor torque
%    B         - [n x n] d(omega)/d(tau) sensitivity (unit-torque responses)
%    tau_motor - [n x 1] Stall torques (signed; sign = driving direction)
%    om_max    - [n x 1] Zero-torque speeds (> 0 enables the law)
%    logf_h    - (Optional) logger handle for diagnostics
%
%  Outputs:
%    tau - [n x 1] Self-consistent effective motor torques
%
%  References:
%    [1] Berg, H. C. & Turner, L. (1993). Biophys. J. 65, 2201-2216.
%    [2] Berg, H. C. (2003). Annu. Rev. Biochem. 72, 19-54.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-13
%  Updated: 2026-09-04 — signed asymmetric linear law; exact linear solve

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
if any(~isfinite(a)) || any(~isfinite(B(:))) || any(~isfinite(tau_motor))
    error('solve_motor_self_consistent:NonFiniteInput', ...
        'Motor self-consistency inputs must be finite.');
end

% Slope of the torque-speed law per motor: d(tau)/d(omega) = -|tau_motor|/om_max
% (in the signed omega variable: tau = tau_motor - s*tau_motor/om_max * omega,
%  with s = sign(tau_motor), so s*tau_motor = |tau_motor|).
% Motors with om_max <= 0 or non-finite om_max are unlimited (slope 0).
limited = isfinite(om_max) & (om_max > 0);
slope = zeros(n, 1);
slope(limited) = abs(tau_motor(limited)) ./ om_max(limited);

if ~any(limited)
    tau = tau_motor;
    return;
end

%  tau = tau_motor - slope .* (a + B*tau)
%  =>  (I + diag(slope)*B) tau = tau_motor - slope .* a
A_sys = eye(n) + diag(slope) * B;
rhs   = tau_motor - slope .* a;
tau   = A_sys \ rhs;

if any(~isfinite(tau)) || rcond(A_sys) < 1e-14
    if ~isempty(logf_h)
        logf_h(['[ODE] ERROR: motor torque-speed solve is singular or ', ...
                'non-finite (n=%d, rcond=%.3e).\n'], n, rcond(A_sys));
    end
    error('solve_motor_self_consistent:SingularSystem', ...
        ['Motor torque-speed solve is singular or non-finite ', ...
         '(n=%d, rcond=%.3e).'], n, rcond(A_sys));
end

end
