function q = parse_initial_orientation(spec)
%PARSE_INITIAL_ORIENTATION  Convert user-specified initial orientation to quaternion
%
%  Description:
%    Translates a user-friendly orientation specification into the
%    scalar-first unit quaternion q = [q0; q1; q2; q3] used internally
%    by the simulator to represent body orientation.
%
%    Accepted input formats (dispatch by presence of field):
%      struct('quaternion',  [w, x, y, z])
%        - Scalar-first unit quaternion (will be normalised).
%      struct('axis_angle',  struct('axis', [x, y, z], 'angle_rad', theta))
%        - Rotation of angle theta around axis (axis normalised internally).
%      struct('euler_zyx',   [yaw, pitch, roll])    % Tait-Bryan (rad)
%        - Intrinsic rotations: yaw about Z, then pitch about new Y,
%          then roll about new X.  Equivalent to:
%              R = Rz(yaw) * Ry(pitch) * Rx(roll)
%
%  Syntax:
%    q = parse_initial_orientation(spec)
%
%  Input:
%    spec - Struct containing exactly ONE of the above fields, OR empty
%           ([] / empty struct) to signal "unspecified".
%
%  Output:
%    q - [4 x 1] Scalar-first unit quaternion [q0; q1; q2; q3], or []
%        if spec is empty / unspecified.
%
%  Conventions:
%    - Active rotation: v_lab = quat_to_rotmat(q) * v_body.
%    - Consistent with quat_multiply / quat_to_rotmat used elsewhere.
%
%  See also: quat_to_rotmat, quat_multiply
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-05-04

if nargin < 1 || isempty(spec) || (isstruct(spec) && isempty(fieldnames(spec)))
    q = [];
    return;
end

if ~isstruct(spec)
    error('parse_initial_orientation:BadInput', ...
        'initial_orientation must be a struct, got %s.', class(spec));
end

% Count which recognised fields are present
has_quat   = isfield(spec, 'quaternion') && ~isempty(spec.quaternion);
has_axis   = isfield(spec, 'axis_angle') && ~isempty(spec.axis_angle);
has_euler  = isfield(spec, 'euler_zyx')  && ~isempty(spec.euler_zyx);
n_formats  = has_quat + has_axis + has_euler;

if n_formats == 0
    error('parse_initial_orientation:NoFormat', ...
        ['initial_orientation struct has none of the recognised fields ', ...
         '{quaternion, axis_angle, euler_zyx}. Please specify exactly one.']);
end
if n_formats > 1
    error('parse_initial_orientation:AmbiguousFormat', ...
        ['initial_orientation struct has multiple format fields; ', ...
         'please specify exactly ONE of {quaternion, axis_angle, euler_zyx}.']);
end

%% Dispatch
if has_quat
    q = parse_quaternion(spec.quaternion);
elseif has_axis
    q = parse_axis_angle(spec.axis_angle);
else
    q = parse_euler_zyx(spec.euler_zyx);
end

%% Final normalisation (safety)
q = q(:);
q_norm = norm(q);
if q_norm < 1e-30 || ~isfinite(q_norm)
    error('parse_initial_orientation:DegenerateQuaternion', ...
        'Parsed quaternion has zero or non-finite norm (%.3e).', q_norm);
end
q = q / q_norm;

end


%% ==================== Local helpers ====================

function q = parse_quaternion(v)
    v = v(:);
    if numel(v) ~= 4 || ~isnumeric(v) || any(~isfinite(v))
        error('parse_initial_orientation:BadQuaternion', ...
            'quaternion must be a 4-element finite numeric vector [w x y z].');
    end
    q = v;
end


function q = parse_axis_angle(aa)
    if ~isstruct(aa)
        error('parse_initial_orientation:BadAxisAngle', ...
            'axis_angle must be a struct with fields .axis and .angle_rad.');
    end
    if ~isfield(aa, 'axis') || ~isfield(aa, 'angle_rad')
        error('parse_initial_orientation:BadAxisAngle', ...
            'axis_angle requires fields .axis (3-vector) and .angle_rad (scalar).');
    end
    axis_v = aa.axis(:);
    if numel(axis_v) ~= 3 || ~isnumeric(axis_v) || any(~isfinite(axis_v))
        error('parse_initial_orientation:BadAxis', ...
            'axis_angle.axis must be a 3-element finite numeric vector.');
    end
    theta = aa.angle_rad;
    if ~isnumeric(theta) || ~isscalar(theta) || ~isfinite(theta)
        error('parse_initial_orientation:BadAngle', ...
            'axis_angle.angle_rad must be a finite scalar.');
    end

    n_axis = norm(axis_v);
    if n_axis < 1e-30
        error('parse_initial_orientation:ZeroAxis', ...
            'axis_angle.axis has zero magnitude.');
    end
    axis_u = axis_v / n_axis;

    half = theta / 2;
    q = [cos(half); sin(half) * axis_u];
end


function q = parse_euler_zyx(ang)
    ang = ang(:);
    if numel(ang) ~= 3 || ~isnumeric(ang) || any(~isfinite(ang))
        error('parse_initial_orientation:BadEulerZYX', ...
            'euler_zyx must be a 3-element finite numeric vector [yaw pitch roll] in radians.');
    end
    yaw   = ang(1);
    pitch = ang(2);
    roll  = ang(3);

    % Quaternions for each intrinsic rotation (scalar-first)
    qz = [cos(yaw/2);   0; 0; sin(yaw/2)];
    qy = [cos(pitch/2); 0; sin(pitch/2); 0];
    qx = [cos(roll/2);  sin(roll/2); 0; 0];

    % Intrinsic sequence Z-Y-X -> q = qz * qy * qx
    q = quat_multiply(quat_multiply(qz, qy), qx);
end
