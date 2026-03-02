function X = pili_spatial_coords(x_root, theta, phi, rot_k, ell_k, N_k)
%PILI_SPATIAL_COORDS  Pili node positions in lab frame from Euler angles
%
%  Description:
%    Reconstructs spatial positions of pili nodes from the root position,
%    Euler angles, and the pili's rotated basis.  Adapted from Walker
%    (2020) spatial_coords.m.
%
%  Syntax:
%    X = pili_spatial_coords(x_root, theta, phi, rot_k, ell_k, N_k)
%
%  Inputs:
%    x_root - [3 x 1] Root position in lab frame (dimensionless)
%    theta  - [N_k x 1] Polar Euler angles in rot_k basis
%    phi    - [N_k x 1] Azimuthal Euler angles in rot_k basis
%    rot_k  - [3 x 3] Rotation matrix (lab -> basis)
%    ell_k  - Dimensionless pili length (L_pili / L_0)
%    N_k    - Number of segments
%
%  Outputs:
%    X - [(N_k+1) x 3] Node positions in lab frame
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%         (spatial_coords.m; cumulative Euler-angle-to-position integration)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Standardised formatting

ds = ell_k / N_k;

% Tangent vectors in rot_k basis, then transform to lab frame
[~, ~, d3_basis] = pili_directors(theta, phi, zeros(N_k, 1));
d3_lab = (rot_k' * d3_basis')';    % [N_k x 3]

% Cumulative integration from root
X = zeros(N_k + 1, 3);
X(1, :) = x_root(:)';
for j = 1:N_k
    X(j+1, :) = X(j, :) + ds * d3_lab(j, :);
end

end
