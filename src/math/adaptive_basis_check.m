function [Z_new, params_new, rotated_flags] = adaptive_basis_check(Z, params)
%ADAPTIVE_BASIS_CHECK  Check and update Euler angle bases for all pili
%
%  Description:
%    For each pilus, checks whether any theta is near the singularity at
%    0 or pi.  If so, selects a new optimal basis that moves the
%    singularity away from the current configuration and transforms the
%    Euler angles accordingly.
%
%    Implements per-pili independent adaptive basis rotation from
%    Walker (2020), main_physical.m lines 153-174.
%
%  Syntax:
%    [Z_new, params_new, rotated_flags] = adaptive_basis_check(Z, params)
%
%  Inputs:
%    Z      - [n x 1] Current state vector
%    params - Parameter structure
%
%  Outputs:
%    Z_new         - [n x 1] Updated state vector
%    params_new    - Updated params (new rot{k} for rotated pili)
%    rotated_flags - [M x 1] Logical, true if pilus k's basis was rotated
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06

%% Initialise
delta  = params.delta;
M      = params.M_pili;
N_segs = params.N_segments;

Z_new         = Z;
params_new    = params;
rotated_flags = false(M, 1);

if M == 0
    return;
end

%% Pre-compute test points on unit sphere
phis_test   = linspace(0, 2*pi, 100);
thetas_test = linspace(0, pi,   100);
[thetas_grid, phis_grid] = meshgrid(thetas_test, phis_test);
thetas_flat = thetas_grid(:);
phis_flat   = phis_grid(:);
v1s = sin(thetas_flat) .* cos(phis_flat);
v2s = sin(thetas_flat) .* sin(phis_flat);
v3s = cos(thetas_flat);

%% Check each pilus
[~, ~, pili_angles] = unpack_state(Z, N_segs);
idx = 8;

for k = 1:M
    N_k     = N_segs(k);
    theta_k = pili_angles{k}.theta;

    % Check proximity to singularity
    if min(pi - theta_k) < delta || min(theta_k) < delta
        rotated_flags(k) = true;

        phi_k = pili_angles{k}.phi;
        psi_k = pili_angles{k}.psi;
        [d1_k, d2_k, d3_k] = pili_directors(theta_k, phi_k, psi_k);

        % Find optimal new z-axis direction
        x_test = [d3_k(:,1); -d3_k(:,1)];
        y_test = [d3_k(:,2); -d3_k(:,2)];
        z_test = [d3_k(:,3); -d3_k(:,3)];

        d_dist = (x_test - v1s').^2 + (y_test - v2s').^2 + (z_test - v3s').^2;
        [~, ind] = max(min(d_dist, [], 1));
        v_opt = [v1s(ind), v2s(ind), v3s(ind)];

        % Average with z-axis and normalise (Walker's smoothing)
        v_opt = mean([v_opt; 0, 0, 1]);
        v_opt_norm = norm(v_opt);
        if v_opt_norm < 1e-12
            % Degenerate case: v_opt was nearly anti-parallel to z-axis.
            % Fall back to x-axis as rotation axis (guaranteed safe since
            % the singularity axis is z).
            v_opt = [1, 0, 0];
        else
            v_opt = v_opt / v_opt_norm;
        end

        % Basis rotation by pi about the axis v_opt:
        %   K = skew(v_opt),  rot_local = I + 2*K^2 = 2*v_opt*v_opt' - I
        %   For a UNIT vector v_opt this is a PROPER rotation of angle pi
        %   about v_opt (eigenvalues {+1, -1, -1}, so det = +1) -- NOT a
        %   reflection.  It maps the basis z-axis (the Euler-angle
        %   singularity location) to its image mirrored across v_opt,
        %   relocating the singularity away from the current tangents while
        %   PRESERVING frame handedness (the right-handed triad d1 x d2 = d3
        %   is retained, so curvature/twist signs are unaffected).  Walker
        %   (2020) uses this form for its simplicity and numerical stability.
        K = [0, -v_opt(3), v_opt(2); ...
             v_opt(3), 0, -v_opt(1); ...
            -v_opt(2), v_opt(1), 0];
        rot_local = eye(3) + 2 * K^2;

        rot_k_old = params.rot{k};
        rot_k_new = rot_local * rot_k_old;

        % Transform directors to new basis
        d1_new = (rot_local * d1_k')';
        d2_new = (rot_local * d2_k')';
        d3_new = (rot_local * d3_k')';

        [theta_new, phi_new, psi_new] = pili_euler_angles(d1_new, d2_new, d3_new);

        % Post-rotation safety check: verify new thetas are away from
        % singularity.  If not, issue a warning (the next ODE segment
        % will re-trigger the event and retry).
        new_margin = min([theta_new; pi - theta_new]);
        if new_margin < delta / 2
            logf_h = get_basis_logf(params);
            logf_h(sprintf(['[ODE] WARNING: Pili %d basis rotation insufficient. ' ...
                'Post-rotation min(theta, pi-theta) = %.4f < delta/2 = %.4f. ' ...
                'ODE event will re-trigger.\n'], k, new_margin, delta/2));
            warning_msg = sprintf( ...
                'Pili %d: basis rotation was insufficient. Post-rotation min(theta, pi-theta) = %.4f < delta/2 = %.4f. ODE event will re-trigger.', ...
                k, new_margin, delta/2);
            warning('adaptive_basis_check:InsufficientRotation', warning_msg);
        end

        % Update state vector
        Z_new(idx        : idx+N_k-1)   = theta_new;
        Z_new(idx+N_k    : idx+2*N_k-1) = phi_new;
        Z_new(idx+2*N_k  : idx+3*N_k-1) = psi_new;

        % Update rotation matrix (re-orthogonalise to prevent accumulated
        % floating-point drift from repeated matrix multiplications).
        [U_rot, ~, V_rot] = svd(rot_k_new);
        rot_k_new = U_rot * V_rot';
        params_new.rot{k} = rot_k_new;

        % Intrinsic curvature is basis-invariant (dot products preserved)
        params_new.kb0{k} = params.kb0{k};
    end

    idx = idx + 3 * N_k;
end

end


%% ==================== Local Helper Function ====================

function logf_h = get_basis_logf(params)
%GET_BASIS_LOGF  Retrieve logf handle from params, or fallback to fprintf
    if isstruct(params) && isfield(params, 'logf') && ~isempty(params.logf)
        logf_h = params.logf;
    else
        logf_h = @(varargin) fprintf(varargin{:});
    end
end
