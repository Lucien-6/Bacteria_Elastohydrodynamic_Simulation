function kb = compute_curvature_binormal(x)
%COMPUTE_CURVATURE_BINORMAL  Discrete curvature binormal vectors for an elastic rod
%
%  Description:
%    Computes the discrete curvature binormal vectors for an inextensible
%    rod using the formula from Bergou et al. (2008) "Discrete Elastic
%    Rods".  The curvature binormal is defined at interior vertices and
%    captures both the magnitude and direction of curvature.
%
%  Syntax:
%    kb = compute_curvature_binormal(x)
%
%  Inputs:
%    x  - [3 x N] Node positions [m], where N is the number of vertices
%
%  Outputs:
%    kb - [3 x (N-2)] Curvature binormal vectors at interior vertices.
%         kb(:,i) corresponds to vertex i+1 (1-indexed).
%
%  Formula:
%    (kappa * b)_i = 2 * (e_{i-1} x e_i)
%                    / (|e_{i-1}| * |e_i| + e_{i-1} . e_i)
%
%  Properties:
%    - |kb| = 2 * tan(theta/2), where theta is the turning angle
%    - kb is perpendicular to the osculating plane
%    - kb = 0 for collinear edges (straight rod)
%
%  References:
%    [1] Bergou, M. et al. (2008). Discrete Elastic Rods. ACM TOG 27(3):63.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-01-23
%  Updated: 2026-02-11 — Renamed to snake_case; standardized formatting

%% Input validation
if size(x, 1) ~= 3
    error('compute_curvature_binormal:InvalidInput', ...
          'Input x must be a 3 x N matrix.');
end

N = size(x, 2);  % Number of vertices

if N < 3
    kb = zeros(3, 0);
    return;
end

%% Compute edge vectors and lengths
e   = x(:, 2:N) - x(:, 1:N-1);   % [3 x (N-1)]
len = vecnorm(e);                  % [1 x (N-1)]

% Numerical safety: prevent division by zero
min_length = 1e-20;
len(len < min_length) = min_length;

%% Compute curvature binormal at each interior vertex
num_interior = N - 2;
kb = zeros(3, num_interior);

for i = 1:num_interior
    e_prev   = e(:, i);       % Edge before interior vertex i+1
    e_curr   = e(:, i + 1);   % Edge after  interior vertex i+1
    len_prev = len(i);
    len_curr = len(i + 1);

    cross_e = cross(e_prev, e_curr);
    dot_e   = dot(e_prev, e_curr);
    denom   = len_prev * len_curr + dot_e;

    % Protect against nearly anti-parallel edges (theta -> pi)
    if abs(denom) < 1e-12
        denom = 1e-12 * sign(denom + 1e-20);
    end

    kb(:, i) = 2 * cross_e / denom;
end

end
