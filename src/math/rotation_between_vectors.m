function R = rotation_between_vectors(u, v)
%ROTATION_BETWEEN_VECTORS  Rotate one direction onto another
%
%  Description:
%    Computes the shortest-arc proper rotation R satisfying
%      R*u_hat = v_hat
%    for finite, nonzero three-vectors u and v.  Inputs are reshaped to
%    columns and divided by max(norm(vector), 1e-30).  The 1e-30 floor
%    prevents division by zero but does not define a meaningful direction:
%    callers must not use zero or numerically vanishing vectors.
%
%    For the general case, the routine applies Rodrigues' formula about
%    cross(u_hat,v_hat).  Dot products above 1 - 1e-12 are treated as
%    parallel and return the identity.  Dot products below -1 + 1e-12 are
%    treated as anti-parallel and return a pi rotation about a deterministic
%    perpendicular axis.  The anti-parallel rotation is necessarily
%    non-unique.  These 1e-12 angular tests avoid normalising a vanishing
%    cross product near either limiting case.
%
%  Syntax:
%    R = rotation_between_vectors(u, v)
%
%  Inputs:
%    u - Three-element finite, nonzero source direction
%    v - Three-element finite, nonzero target direction
%
%  Outputs:
%    R - [3 x 3] Proper rotation mapping u_hat to v_hat, within round-off
%
%  References:
%    [1] Rodrigues' rotation formula.
%        See: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-28
%  Updated: 2026-07-12

%% Normalise input vectors
u = u(:) / max(norm(u), 1e-30);
v = v(:) / max(norm(v), 1e-30);

%% Compute cosine of the angle between u and v
c = dot(u, v);

%% Handle special cases
if c > 1.0 - 1e-12
    % u and v are nearly parallel: identity rotation
    R = eye(3);
    return;
end

if c < -1.0 + 1e-12
    % u and v are nearly anti-parallel: 180-degree rotation
    % Find a vector perpendicular to u to use as rotation axis
    if abs(u(1)) < 0.9
        axis = cross(u, [1; 0; 0]);
    else
        axis = cross(u, [0; 1; 0]);
    end
    axis = axis / norm(axis);

    % 180-degree Rodrigues formula: R = 2 * k * k' - I
    R = 2 * (axis * axis') - eye(3);
    return;
end

%% General case: Rodrigues formula
%  Rotation axis: k = u x v (normalised)
k = cross(u, v);
k = k / norm(k);

% Skew-symmetric matrix of k
K = [0,    -k(3),  k(2);
     k(3),  0,    -k(1);
    -k(2),  k(1),  0   ];

% Rodrigues formula: R = I + sin(theta) * K + (1 - cos(theta)) * K^2
s = sqrt(1 - c^2);     % sin(theta)
R = eye(3) + s * K + (1 - c) * (K * K);

end
