function M_active = coupled_active_moment(t, k, N_k, d1, d2, d3, d1N, d2N, d3N, am_cfg, root_coupling_k)
%COUPLED_ACTIVE_MOMENT  Per-pili time-varying active internal moments
%
%  Description:
%    Computes the integrated active internal moments for pilus k at
%    dimensionless time t.  The result is expressed in the local
%    director basis (d1N, d2N, d3N) at each node, matching the layout
%    of elastic curvature entries in R_global.
%
%    The algorithm follows 3D-Flexible-Filaments/integrated_internal_moments.m:
%      1. Specify moment density as antiderivatives f(s) on s in [0,1].
%      2. Evaluate segment integrals: f(s_{j+1}) - f(s_j).
%      3. Reverse cumulative sum to get integral from s_i to 1.
%      4. Project onto local director basis at each node.
%
%    Supported modes:
%      'off'             - Zero active moments (default)
%      'traveling_wave'  - Sinusoidal traveling wave along the pilus
%      'custom'          - User-provided function handle
%
%  Syntax:
%    M_active = coupled_active_moment(t, k, N_k, d1, d2, d3, ...
%                                     d1N, d2N, d3N, am_cfg)
%    M_active = coupled_active_moment(t, k, N_k, d1, d2, d3, ...
%                                     d1N, d2N, d3N, am_cfg, root_coupling_k)
%
%  Inputs:
%    t      - Current dimensionless time
%    k      - Pilus index (for custom mode identification)
%    N_k    - Number of segments for this pilus
%    d1,d2,d3   - [N_k x 3] Segment directors (rot_k basis)
%    d1N,d2N,d3N - [(N_k+1) x 3] Nodal directors (rot_k basis)
%    am_cfg - Active moment configuration struct with fields:
%               .type  - 'off' | 'traveling_wave' | 'custom'
%             For 'traveling_wave':
%               .amplitude_dimless - Dimensionless amplitude (pre-converted)
%               .omega_dimless     - Dimensionless angular frequency
%               .wavenumber        - Spatial wavenumber (on s in [0,1])
%               .phase_rad         - Optional phase offset [rad] in cos(kw*s - w*t + phase_rad)
%                                    (default 0); shifts temporal phase of the wave
%               .d1_amp  - Weight for d1 (normal) bending component
%               .d2_amp  - Weight for d2 (binormal) bending component
%               .d3_amp  - Weight for d3 (tangent) torsional component
%                          Non-zero d3_amp drives twist rather than bending;
%                          typically zero for standard flagellar/pili beating.
%             For 'custom':
%               .fun - Function handle: fun(t, k, N_k, d1, d2, d3, d1N, d2N, d3N)
%    root_coupling_k - (Optional) 'clamped' (default) | 'free_torsion'
%                      Root boundary condition type for this pilus.
%
%  Outputs:
%    M_active - [3 x N_k] Integrated active moments in local director basis.
%               Row 1: d1 component, Row 2: d2 component, Row 3: d3 component.
%               For clamped coupling, first segment is zeroed (overwritten by BC).
%               For free_torsion, only bending components (rows 1-2) of the
%               first segment are zeroed; torsional component (row 3) is preserved.
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] 3D-Flexible-Filaments/integrated_internal_moments.m
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-03

if nargin < 11 || isempty(root_coupling_k)
    root_coupling_k = 'clamped';
end

if strcmp(am_cfg.type, 'off')
    M_active = zeros(3, N_k);
    return
end

switch am_cfg.type
    case 'traveling_wave'
        A    = am_cfg.amplitude_dimless;
        w    = am_cfg.omega_dimless;
        kw   = am_cfg.wavenumber;
        a1   = am_cfg.d1_amp;
        a2   = am_cfg.d2_amp;
        a3   = am_cfg.d3_amp;
        if isfield(am_cfg, 'phase_rad') && ~isempty(am_cfg.phase_rad)
            phi_rad = am_cfg.phase_rad;
        else
            phi_rad = 0;
        end

        s = (0 : 1/N_k : 1)';

        % Antiderivatives of A * a_i * cos(kw*s - w*t + phi_rad):
        %   integral = (A * a_i / kw) * sin(kw*s - w*t + phi_rad)
        phase = @(s_val) kw * s_val - w * t + phi_rad;

        if abs(kw) > 1e-15
            inv_kw = 1.0 / kw;
            d1_anti = @(s_val) A * a1 * inv_kw * sin(phase(s_val));
            d2_anti = @(s_val) A * a2 * inv_kw * sin(phase(s_val));
            d3_anti = @(s_val) A * a3 * inv_kw * sin(phase(s_val));
        else
            d1_anti = @(s_val) A * a1 * s_val .* cos(phi_rad - w * t);
            d2_anti = @(s_val) A * a2 * s_val .* cos(phi_rad - w * t);
            d3_anti = @(s_val) A * a3 * s_val .* cos(phi_rad - w * t);
        end

        seg_d1 = d1_anti(s(2:end)) - d1_anti(s(1:end-1));
        seg_d2 = d2_anti(s(2:end)) - d2_anti(s(1:end-1));
        seg_d3 = d3_anti(s(2:end)) - d3_anti(s(1:end-1));

        seg_moments = seg_d1 .* d1 + seg_d2 .* d2 + seg_d3 .* d3;

        integrated_moments = cumsum(seg_moments, 1, 'reverse');

        M_active = [dot(d1N(1:end-1,:), integrated_moments, 2), ...
                    dot(d2N(1:end-1,:), integrated_moments, 2), ...
                    dot(d3N(1:end-1,:), integrated_moments, 2)]';

    case 'custom'
        M_active = am_cfg.fun(t, k, N_k, d1, d2, d3, d1N, d2N, d3N);

    otherwise
        error('coupled_active_moment:UnknownType', ...
            'Unknown active_moment type: ''%s'' for pili %d.', am_cfg.type, k);
end

if strcmp(root_coupling_k, 'clamped')
    M_active(:, 1) = 0;     % All components zeroed (overwritten by clamped BC)
else
    M_active(1:2, 1) = 0;   % Bending components zeroed (theta_1, phi_1 constrained)
    % M_active(3, 1) preserved — torsional component drives free psi_1 dynamics
end

end
