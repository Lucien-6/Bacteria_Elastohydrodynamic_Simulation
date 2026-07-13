function Z = enforce_clamped_geometry(Z, params, force_full_projection)
%ENFORCE_CLAMPED_GEOMETRY  Correct root-only director drift
%
%  Description:
%    Corrects geometric drift only at each pilus-body junction.  Interior
%    segment angles are never modified.  The state is
%    laid out as
%      Z = [x_body(3); q_body(4); theta_1; phi_1; psi_1; ...],
%    where each pilus contributes consecutive theta, phi, and psi blocks,
%    each of length N_segments(k).  All Euler angles are in radians.
%
%    The required first-segment frame is formed in the lab frame as
%      D_lab = R_body*params.D_local{k}.
%    params.rot{k} maps lab-frame vectors into pilus-k basis coordinates,
%    so D_basis = params.rot{k}*D_lab is converted to Euler angles before
%    the first element of each relevant angle block is overwritten.
%
%    theta_1 and phi_1 are always projected.  psi_1 is also overwritten
%    for a clamped root, or whenever force_full_projection is true.  For
%    a free_torsion root under the default call, psi_1 is left unchanged
%    so that torsional dynamics remain free.  The full projection option
%    is intended for consistent initialisation.
%
%    Call after quaternion renormalisation and after each ODE integration
%    interval.  Only first-segment angles are modified.
%
%  Syntax:
%    Z = enforce_clamped_geometry(Z, params)
%    Z = enforce_clamped_geometry(Z, params, force_full_projection)
%
%  Inputs:
%    Z                     - [n x 1] Global state vector in the layout above
%    params                - Parameter structure containing body-to-root
%                            frames, basis rotations, and segment counts
%    force_full_projection - Optional logical scalar; default false
%
%  Outputs:
%    Z - [n x 1] State with required first-segment angles overwritten
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%        (Clamped-end boundary condition; geometric projection step)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-07
%  Updated: 2026-07-13

if nargin < 3 || isempty(force_full_projection)
    force_full_projection = false;
end

M = params.M_pili;
if M == 0
    return;
end

%% Unpack body orientation
q      = Z(4:7);
R_body = quat_to_rotmat(q);

%% Correct each pilus's first-segment Euler angles
N_segs = params.N_segments;
idx    = 8;

for k = 1:M
    N_k   = N_segs(k);
    rot_k = params.rot{k};

    % Expected director frame in lab frame
    D_lab   = R_body * params.D_local{k};
    D_basis = rot_k * D_lab;

    d1_basis = D_basis(:, 1)';
    d2_basis = D_basis(:, 2)';
    d3_basis = D_basis(:, 3)';

    [theta_c, phi_c, psi_c] = pili_euler_angles(d1_basis, d2_basis, d3_basis);

    % Overwrite first-segment angles based on root coupling type.
    rc_k = 'clamped';
    if isfield(params, 'root_coupling') && numel(params.root_coupling) >= k
        rc_k = params.root_coupling{k};
    end

    Z(idx)       = theta_c;     % theta_1 always projected
    Z(idx + N_k) = phi_c;       % phi_1 always projected
    if strcmp(rc_k, 'clamped') || force_full_projection
        Z(idx + 2*N_k) = psi_c;  % psi_1 projected for clamped, or on demand
    end
    % free_torsion (default call): psi_1 is NOT projected — it evolves freely.

    idx = idx + 3 * N_k;
end

end
