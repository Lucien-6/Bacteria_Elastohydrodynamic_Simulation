function X = pili_spatial_coords(x_root, theta, phi, rot_k, ell_k, N_k)
%PILI_SPATIAL_COORDS  Integrate pilus tangents into lab-frame node positions
%
%  Description:
%    Reconstructs N_k+1 node positions from the root and N_k piecewise
%    constant tangent directions.  theta and phi are in radians and are
%    expressed in the pilus rotated basis.  rot_k maps lab coordinates to
%    basis coordinates, so rot_k' maps each computed tangent back to the
%    lab frame before cumulative integration with ds = ell_k/N_k.
%
%    Position depends only on d3(theta, phi).  The material twist angle psi
%    rotates d1 and d2 about d3 and therefore cannot affect the centreline;
%    pili_directors is consequently evaluated with psi = 0 here.
%
%  Syntax:
%    X = pili_spatial_coords(x_root, theta, phi, rot_k, ell_k, N_k)
%
%  Inputs:
%    x_root - [3 x 1] Dimensionless root position in the lab frame
%    theta  - [N_k x 1] Basis-frame polar angles, in radians
%    phi    - [N_k x 1] Basis-frame azimuthal angles, in radians
%    rot_k  - [3 x 3] Orthogonal map from lab frame to pilus basis
%    ell_k  - Dimensionless pilus length, L_pili/L_0
%    N_k    - Positive integer number of equal-length segments
%
%  Outputs:
%    X - [(N_k+1) x 3] Dimensionless node positions in the lab frame
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%        (spatial_coords.m; cumulative Euler-angle-to-position integration)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-07-12

ds = ell_k / N_k;

% Tangent vectors in rot_k basis, then transform to lab frame
[~, ~, d3_basis] = pili_directors(theta, phi, zeros(N_k, 1));
d3_lab = (rot_k' * d3_basis')';    % [N_k x 3]

% Cumulative integration from root (vectorised)
X = [x_root(:)'; x_root(:)' + ds * cumsum(d3_lab, 1)];

end
