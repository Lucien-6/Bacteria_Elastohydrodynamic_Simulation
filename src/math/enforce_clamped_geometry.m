function Z = enforce_clamped_geometry(Z, params, force_full_projection)
%ENFORCE_CLAMPED_GEOMETRY  Project first-segment Euler angles to match body
%
%  Description:
%    Corrects constraint drift at the pili-body junction by enforcing:
%
%        D_lab^{(1)} = R_body * D_local{k}
%
%    For clamped pili, all three first-segment Euler angles (theta_1,
%    phi_1, psi_1) are projected.  For free_torsion pili, only theta_1
%    and phi_1 are projected by default; psi_1 is left free to evolve
%    under the torsional dynamics (elastic moment + motor torque).
%
%    The optional `force_full_projection` flag lets callers override the
%    free_torsion behaviour and project all three angles (including
%    psi_1) regardless of root coupling type.  This is typically used at
%    initialisation (`init_coupled_state`) so that the starting state
%    already satisfies the full clamped geometry — the first ODE step
%    can then evolve psi_1 freely from a physically self-consistent
%    initial condition.
%
%    Call after each ODE integration interval and after quaternion
%    renormalisation, with `force_full_projection` left at its default
%    value (false) so that free_torsion ψ_1 is preserved during the run.
%
%  Syntax:
%    Z = enforce_clamped_geometry(Z, params)
%    Z = enforce_clamped_geometry(Z, params, force_full_projection)
%
%  Inputs:
%    Z                       - [n x 1] Global state vector
%    params                  - Parameter structure (from init_coupled_params)
%    force_full_projection   - (Optional, default false) If true, project
%                              all three first-segment Euler angles
%                              (theta_1, phi_1, psi_1) regardless of
%                              root_coupling.  Intended for initialisation.
%
%  Outputs:
%    Z - [n x 1] State vector with corrected first-segment angles
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%         (Clamped-end boundary condition; geometric projection step)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-07

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
