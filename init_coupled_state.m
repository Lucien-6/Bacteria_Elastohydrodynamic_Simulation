function Z0 = init_coupled_state(params)
%INIT_COUPLED_STATE  Initialise the global state vector Z0
%
%  Description:
%    Creates the initial state vector for the coupled body-pili system.
%    The body starts at the origin with orientation determined by the
%    initial tilt angle alpha_0 (angle of body long axis from +X in the
%    XZ plane).  Each pilus starts in its rest (intrinsic) configuration,
%    rotated consistently with the body.
%
%  Syntax:
%    Z0 = init_coupled_state(params)
%
%  Inputs:
%    params - Parameter structure (from init_coupled_params)
%             Must include params.alpha_0 (tilt angle in rad)
%
%  Outputs:
%    Z0 - [n x 1] Initial global state vector
%           Z0(1:3)   = X_cm  (body position, dimensionless)
%           Z0(4:7)   = q     (body quaternion)
%           Z0(8:end)         = pili Euler angles
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-12 — Added alpha_0 initial tilt angle support

%% Body initial state
X_cm_0 = [0; 0; 0];

% Compute initial quaternion from tilt angle alpha_0
%   alpha_0 : angle between body long axis and +X axis in XZ plane
%   Body long axis is X in body frame (capsule_body_building convention).
%   To map body-X -> (cos(alpha_0), 0, sin(alpha_0)) in lab frame,
%   rotate about Y by theta_Y = -alpha_0.
%
%   Proof:  R_Y(theta_Y) * [1;0;0] = [cos(theta_Y); 0; -sin(theta_Y)]
%           Setting this equal to [cos(alpha_0); 0; sin(alpha_0)] gives
%           theta_Y = -alpha_0.
%
%   Quaternion: q = [cos(theta_Y/2); 0; sin(theta_Y/2); 0]
%
%   Special cases:
%     alpha_0 =  pi/2  -> theta_Y = -pi/2 -> vertical (+Z)
%     alpha_0 =  0     -> theta_Y = 0      -> horizontal (+X)
%     alpha_0 = -pi/2  -> theta_Y = pi/2   -> inverted (-Z)
alpha_0   = params.alpha_0;
theta_Y   = -alpha_0;
q_0       = [cos(theta_Y/2); 0; sin(theta_Y/2); 0];
q_0       = q_0 / norm(q_0);              % Ensure unit quaternion

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

    % Original pili coords are in the default body frame (long axis along Z).
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

%% Verify reconstruction consistency
geom = reconstruct_geometry(Z0, params);

for k = 1:M
    X_recon = geom.X_pili{k};
    % Compare against rotated original coords (matching body orientation)
    X_orig_lab = (R_body * params.pili_x0{k}')';
    max_err = max(abs(X_recon - X_orig_lab), [], 'all');
    if max_err > 1e-6
        warning('init_coupled_state:GeometryMismatch', ...
            'Pili %d reconstruction error: %.2e (should be < 1e-6)', k, max_err);
    end
end

fprintf('Initial state vector created: %d DOFs\n', n_total);
fprintf('  Body tilt: alpha_0 = %.4f rad (%.1f deg), q_0 = [%.4f; %.4f; %.4f; %.4f]\n', ...
    alpha_0, rad2deg(alpha_0), q_0(1), q_0(2), q_0(3), q_0(4));

end
