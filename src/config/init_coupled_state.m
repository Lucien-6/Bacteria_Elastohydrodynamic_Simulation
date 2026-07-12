function Z0 = init_coupled_state(params, logf_h)
%INIT_COUPLED_STATE  Initialise the global state vector Z0
%
%  Description:
%    Creates the initial state vector for the coupled body-pili system.
%    The body starts at the origin with the pre-resolved orientation in
%    params.q_0_init.  Each pilus starts in its reference centreline
%    configuration, rotated consistently with the body, and its first
%    director frame is projected onto the root boundary condition.
%
%  Syntax:
%    Z0 = init_coupled_state(params)
%    Z0 = init_coupled_state(params, logf_h)
%
%  Inputs:
%    params  - Parameter structure (from init_coupled_params)
%              Must include q_0_init, M_pili, N_segments, rot, pili_x0,
%              and D_local, together with fields required by
%              enforce_clamped_geometry and reconstruct_geometry.
%              q_0_init is a scalar-first quaternion [w;x;y;z] and is
%              normalised before use.  When an explicit initial_orientation
%              was resolved, params.alpha_0 may be NaN; q_0_init remains the
%              sole authoritative orientation input.
%    logf_h  - (Optional) Logger function handle. If omitted or empty,
%              falls back to fprintf.
%
%  Output:
%    Z0 - [7 + 3*sum(params.N_segments) x 1] state vector:
%           Z0(1:3)   = X_cm  (body position, dimensionless)
%           Z0(4:7)   = q     (body quaternion)
%           Z0(8:end)         = per-pilus [theta; phi; psi] segment angles
%
%  Warnings and Side Effects:
%    Logs and emits init_coupled_state:GeometryMismatch when reconstructed
%    pili coordinates differ from their rotated references by more than
%    1e-6.  Logs and emits init_coupled_state:DirectorMismatch when a root
%    director frame differs from its boundary-condition target by more than
%    1e-10.
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%         (Quaternion-based body orientation; pili initial configuration)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-07-12

if nargin < 2 || isempty(logf_h)
    logf_h = @(varargin) fprintf(varargin{:});
end

%% Body initial state
X_cm_0 = [0; 0; 0];

% Initial quaternion is pre-resolved in init_coupled_params:
%   - Parametric mode        : derived from alpha_0 (Y-axis tilt)
%   - Custom mode, alpha_0   : same Y-axis tilt convention
%   - Custom mode, explicit  : user-supplied initial_orientation
%                              (quaternion / axis_angle / euler_zyx);
%                              params.alpha_0 is NaN on this path
% See parse_initial_orientation.m for the custom-mode spec.
q_0 = params.q_0_init(:);
q_0 = q_0 / max(norm(q_0), eps);          % Ensure unit quaternion

%% Build state vector
M      = params.M_pili;
N_segs = params.N_segments;

n_total = 7 + 3 * sum(N_segs);
Z0 = zeros(n_total, 1);

Z0(1:3) = X_cm_0;
Z0(4:7) = q_0;

%% Pili initial Euler angles
idx    = 8;
R_body = quat_to_rotmat(q_0);

for k = 1:M
    N_k   = N_segs(k);
    rot_k = params.rot{k};

    % Pili coords in params.pili_x0{k} are in the body-local frame
    % (parametric: body-X is the long axis; custom: user-defined axes).
    % Rotate them to the lab frame matching the initial body orientation.
    coords_k_orig = params.pili_x0{k};              % [(N_k+1) x 3]
    coords_k = (R_body * coords_k_orig')';          % [(N_k+1) x 3] in lab frame

    tangents = diff(coords_k, 1, 1);
    seg_lens = sqrt(sum(tangents.^2, 2));
    tangents = tangents ./ seg_lens;

    % Transform tangents to pili basis
    tangents_basis = (rot_k * tangents')';

    theta_k = acos(max(-1, min(1, tangents_basis(:, 3))));
    phi_k   = atan2(tangents_basis(:, 2), tangents_basis(:, 1));
    psi_k   = zeros(N_k, 1);

    Z0(idx        : idx+N_k-1)   = theta_k;
    Z0(idx+N_k    : idx+2*N_k-1) = phi_k;
    Z0(idx+2*N_k  : idx+3*N_k-1) = psi_k;
    idx = idx + 3*N_k;
end

%% Enforce clamped geometry: correct first-segment Euler angles
%  The first-segment directors must satisfy D_lab = R_body * D_local{k}.
%  Setting psi = 0 above uses the Walker (2020) convention default, which
%  does NOT account for the body rotation when alpha_0 != 0.  The enforce
%  step overwrites the first-segment (theta, phi, psi) to be exactly
%  consistent with the clamped boundary condition at the position level.
%
%  We pass force_full_projection = true so that ψ_1 is projected even for
%  free_torsion pili.  Otherwise, with alpha_0 ≠ 0 (or any non-identity
%  initial body orientation), the first-segment director frame would not
%  match the BC requirement at t = 0, triggering a spurious torsional
%  elastic energy that "releases" during the first ODE step.  At runtime
%  (run_simulation main loop, post-quaternion-renorm, post-Brownian-kick)
%  the function is called WITHOUT this flag so free_torsion ψ_1 keeps
%  evolving freely as designed.
Z0 = enforce_clamped_geometry(Z0, params, true);

%% Verify reconstruction consistency
geom = reconstruct_geometry(Z0, params);

for k = 1:M
    X_recon = geom.X_pili{k};
    % Compare against rotated original coords (matching body orientation)
    X_orig_lab = (R_body * params.pili_x0{k}')';
    max_err = max(abs(X_recon - X_orig_lab), [], 'all');
    if max_err > 1e-6
        logf_h('[INIT] WARNING: Pili %d reconstruction error: %.2e (should be < 1e-6)\n', k, max_err);
        warning('init_coupled_state:GeometryMismatch', ...
            'Pili %d reconstruction error: %.2e (should be < 1e-6)', k, max_err);
    end
end

% Verify director consistency with clamped BC
for k = 1:M
    rot_k = params.rot{k};
    D_lab_expected = R_body * params.D_local{k};
    D_basis_expected = rot_k * D_lab_expected;

    d1b = geom.d1_basis{k}(1,:)';
    d2b = geom.d2_basis{k}(1,:)';
    d3b = geom.d3_basis{k}(1,:)';
    D_basis_actual = [d1b, d2b, d3b];

    dir_err = max(abs(D_basis_expected - D_basis_actual), [], 'all');
    if dir_err > 1e-10
        logf_h('[INIT] WARNING: Pili %d first-segment director error: %.2e (should be < 1e-10)\n', k, dir_err);
        warning('init_coupled_state:DirectorMismatch', ...
            'Pili %d first-segment director error: %.2e (should be < 1e-10)', k, dir_err);
    end
end

end
