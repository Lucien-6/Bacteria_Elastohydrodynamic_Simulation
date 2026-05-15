function R = rotation_between_vectors(u, v)
%ROTATION_BETWEEN_VECTORS  Compute rotation matrix aligning vector u to v
%
%  Description:
%    Computes the minimal (shortest-arc) rotation matrix R that maps
%    unit vector u to unit vector v using the Rodrigues formula.
%    Both inputs are normalised internally.  Handles parallel and
%    anti-parallel vectors as special cases.
%
%  Syntax:
%    R = rotation_between_vectors(u, v)
%
%  Inputs:
%    u - [3 x 1] Source direction vector (will be normalised)
%    v - [3 x 1] Target direction vector (will be normalised)
%
%  Output:
%    R - [3 x 3] Rotation matrix such that R * u_hat = v_hat
%        (up to numerical precision)
%
%  References:
%    [1] Rodrigues' rotation formula.
%        See: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-28

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
