function params = init_coupled_params(varargin)
%INIT_COUPLED_PARAMS  Initialise parameter structure for coupled simulation
%
%  Description:
%    Sets up the parameter structure containing all physical and numerical
%    parameters for the coupled body-pili simulation.  Uses the unified
%    nondimensionalisation scheme from the integration plan.
%
%  Syntax:
%    params = init_coupled_params()
%    params = init_coupled_params('Name', Value, ...)
%
%  Name-Value Pairs:
%    'geometry_mode' Geometry mode: 'parametric' (default) or 'custom'
%    'custom_geometry' Struct with custom geometry config (required if 'custom')
%    'N_body'      Number of body surface points      (default: 300)
%                  [parametric mode only]
%    'minor_axis'  Semi-minor axis [m]                 (default: 0.4e-6)
%                  [parametric mode only]
%    'major_axis'  Semi-major axis [m]                 (default: 1.25e-6)
%                  [parametric mode only; L_0 source]
%    'mu'          Dynamic viscosity [Pa*s]            (default: 1e-3)
%    'E'           Young's modulus [Pa]                (default: 3.279e8)
%    'r_pili'      Pili radius [m]                     (default: 3e-9)
%    'sigma'       Poisson's ratio                     (default: 0.2)
%    'rho_bac'     Bacteria density [kg/m^3]           (default: 1180)
%    'rho_fluid'   Fluid density [kg/m^3]              (default: 998.232)
%    'N_segments'  Segment counts [M x 1]              (default: 300)
%                  [parametric mode only]
%    'Pili_Matrix' Pili parameter matrix [8 x M]       (default: single pole)
%                  [parametric mode only]
%    'seg_per_um'  Segments per micron                 (default: 10)
%    'delta'       Euler-angle singularity threshold   (default: pi/20)
%    'alpha_0'     Initial body tilt angle [rad]       (default: pi/2)
%                  Angle of body long axis (body X) from +X in XZ plane
%                  pi/2 = vertical (+Z), 0 = horizontal (+X),
%                  NOTE: the inline demo in main_coupled_simulation.m sets
%                  alpha_0 = 0 (horizontal), which is intentionally non-default
%                  to isolate motor-driven rotation; the README sinking
%                  workflow uses the default pi/2 (vertical).
%                  -pi/2 = inverted (-Z).  Range: [-pi/2, pi/2]
%    'Temper'      Absolute temperature [K]            (default: 293.15)
%    'enable_brown' Enable body Brownian motion        (default: false)
%    'eps_body'      Optional RSM regularisation length [L0], dimensionless.
%                    If empty (default), uses eps_body = h_mean/2.6667 from the
%                    surface mesh.  If set to a positive scalar, overrides
%                    that default (for parameter studies).
%
%  Outputs:
%    params - Structure with all simulation parameters
%
%  Note:
%    'dis' (gap between body and pili root) is computed automatically
%    as dis = 1.1 * (eps_body + eps_pili) to ensure proper separation.
%
%    Pili use the RSS (Regularized Stokeslet Segments) method
%    (Cortez 2018) with exact segment integrals + rotlet coupling.
%    Body uses RSM (Cortez 2005) point forces on surface nodes.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06

%% Parse inputs
p = inputParser;
addParameter(p, 'geometry_mode',  'parametric', @(x) any(validatestring(x, {'parametric', 'custom'})));
addParameter(p, 'custom_geometry', struct(),    @isstruct);
addParameter(p, 'N_body',      300,    @isnumeric);
addParameter(p, 'minor_axis',  0.4e-6, @isnumeric);
addParameter(p, 'major_axis',  1.25e-6, @isnumeric);
addParameter(p, 'mu',          1e-3,   @isnumeric);
addParameter(p, 'E',           3.279e8, @isnumeric);
addParameter(p, 'r_pili',      3e-9,   @isnumeric);
addParameter(p, 'sigma',       0.2,    @isnumeric);
addParameter(p, 'rho_bac',     1180.0, @isnumeric);
addParameter(p, 'rho_fluid',   998.232, @isnumeric);
addParameter(p, 'N_segments',  300,    @isnumeric);
addParameter(p, 'Pili_Matrix', [1; 3; 0; 0; NaN; NaN; NaN; 1], @isnumeric);
addParameter(p, 'seg_per_um',  10,     @isnumeric);
addParameter(p, 'delta',       pi/20,  @isnumeric);
addParameter(p, 'alpha_0',    pi/2,   @isnumeric);
addParameter(p, 'initial_orientation', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'Temper',      293.15,  @isnumeric);
addParameter(p, 'enable_brown', false,  @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'intrinsic_curvature', {}, @iscell);
addParameter(p, 'active_moment', {}, @iscell);
addParameter(p, 'root_coupling', {}, @iscell);
addParameter(p, 'motor_torque', {}, @iscell);
addParameter(p, 'background_flow', struct('type','off'), @isstruct);
addParameter(p, 'eps_body', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
addParameter(p, 'logf', @(varargin) fprintf(varargin{:}), @(x) isa(x, 'function_handle'));
addParameter(p, 'logf_file', @(varargin) [], @(x) isa(x, 'function_handle'));
parse(p, varargin{:});

geometry_mode = p.Results.geometry_mode;

%% Store dimensional parameters
params.dim.minor_axis = p.Results.minor_axis;
params.dim.major_axis = p.Results.major_axis;
params.dim.mu         = p.Results.mu;
params.dim.E          = p.Results.E;
params.dim.r_pili     = p.Results.r_pili;
params.dim.rho_bac    = p.Results.rho_bac;
params.dim.rho_fluid  = p.Results.rho_fluid;
params.dim.Temper     = p.Results.Temper;
params.brown_enabled  = p.Results.enable_brown;
params.geometry_mode  = geometry_mode;

% Store logger handles (injected by run_simulation; fallback to fprintf)
params.logf      = p.Results.logf;
params.logf_file = p.Results.logf_file;
logf_h = params.logf;   % shorthand for frequent use below

%% Characteristic scales for nondimensionalisation
switch geometry_mode
    case 'parametric'
        L_0 = p.Results.major_axis;                       % Length scale [m]
    case 'custom'
        custom_cfg = p.Results.custom_geometry;
        if isfield(custom_cfg, 'L_0_override') && ~isempty(custom_cfg.L_0_override)
            L_0 = custom_cfg.L_0_override;
        else
            % Will be set after reading body points
            L_0 = NaN;
        end
end
I_pili = pi * p.Results.r_pili^4 / 4;             % Second moment of area [m^4]
EI  = p.Results.E * I_pili;                        % Bending stiffness [N*m^2]

%% Build body geometry (dispatch by geometry_mode)
switch geometry_mode
    case 'parametric'
        % --- Parametric capsule body ---
        N_body_requested = p.Results.N_body;
        [Capsule_str, N_body_actual] = capsule_body_building( ...
            N_body_requested, p.Results.minor_axis, p.Results.major_axis);

        % Store body surface points in LOCAL frame (nondimensionalised)
        params.N_body    = N_body_actual;
        body_coords_m    = [Capsule_str.rxB'; Capsule_str.ryB'; Capsule_str.rzB'];

        % Capsule (spherocylinder) volume
        a = p.Results.major_axis;
        b = p.Results.minor_axis;
        V_body = pi * b^2 * (4/3*b + 2*(a - b));

        params.body_type   = 'capsule';
        params.body_faces  = [];       % No triangulation needed; surf() formula used
        params.body_verts  = [];
        params.dim.volume  = V_body;

        % L_0 is major_axis for parametric mode
        L_0 = p.Results.major_axis;

    case 'custom'
        % --- Custom body from point cloud ---
        [Body_str, body_meta] = custom_body_building(custom_cfg);

        N_body_actual = body_meta.N_body;
        params.N_body = N_body_actual;
        body_coords_m = [Body_str.rxB'; Body_str.ryB'; Body_str.rzB'];

        V_body = body_meta.volume;

        % Guard: surface triangulation failure yields NaN volume,
        % which would silently corrupt F_sink and Eg downstream.
        if isnan(V_body) || V_body <= 0
            error('init_coupled_params:InvalidVolume', ...
                ['Custom body triangulation produced an invalid volume (%.4e m^3). ', ...
                 'Surface reconstruction (alphaShape / convhulln) failed or the point ', ...
                 'cloud does not define a closed surface.  Verify the body point cloud file.'], ...
                V_body);
        end

        params.body_type   = 'custom';
        params.body_faces  = body_meta.faces;
        params.body_verts  = body_meta.vertices;         % [m] body-local frame
        params.custom_attachment_sites = body_meta.attachment_sites;
        params.dim.volume  = V_body;

        % Persist the alphaShape for downstream visualisation (flow-field
        % SDF masking in custom mode).  Empty if convhulln fallback was used.
        if isfield(body_meta, 'alpha_shape')
            params.body_alpha_shape = body_meta.alpha_shape;
        else
            params.body_alpha_shape = [];
        end

        % L_0: use override if provided, else use suggested
        if isnan(L_0)
            L_0 = body_meta.L_0_suggested;
        end

        % Store custom-specific spacing info (before nondim)
        params.dim.custom_h_mean = body_meta.h_mean;
        params.dim.custom_h_min  = body_meta.h_min;
end

% Nondimensionalise body coordinates
params.body_local = body_coords_m / L_0;

% Nondimensionalised triangulation vertices for visualisation modules
% (compute_plot_flowfield uses these for SDF masking in custom mode).
if strcmp(geometry_mode, 'custom') && ~isempty(params.body_verts)
    params.body_verts_local = params.body_verts / L_0;   % [L_0]
else
    params.body_verts_local = [];
end

% Store L_0 and derived scales
T_0 = 8 * pi * p.Results.mu * L_0^4 / EI;         % Time scale [s]
U_0 = L_0 / T_0;                                   % Velocity scale [m/s]
F_0 = EI / L_0^2;                                  % Force scale [N]

params.scales.L_0 = L_0;
params.scales.T_0 = T_0;
params.scales.U_0 = U_0;
params.scales.F_0 = F_0;
params.scales.EI  = EI;

%% Compute sinking force and dimensionless parameters
Delta_rho = p.Results.rho_bac - p.Results.rho_fluid;
g         = 9.81;
F_sink    = Delta_rho * V_body * g;             % Sinking force [N]

params.Eg      = F_sink * L_0^2 / EI;             % Elastogravitational number
params.epsilon = p.Results.r_pili / L_0;           % Global slenderness
params.sigma   = p.Results.sigma;
params.g_hat   = [0; 0; -1];                       % Gravity direction
params.delta   = p.Results.delta;

%% Resolve initial body orientation
%  Priority: initial_orientation (if non-empty) overrides alpha_0.
%  In parametric mode, initial_orientation is always empty, so the
%  alpha_0-derived quaternion is used (legacy behaviour).
alpha_0_input         = p.Results.alpha_0;
init_orient_input     = p.Results.initial_orientation;
orientation_source    = 'alpha_0';

if ~isempty(init_orient_input)
    % Only meaningful in custom mode
    if ~strcmp(geometry_mode, 'custom')
        error('init_coupled_params:OrientationModeMismatch', ...
            'initial_orientation is only supported in geometry_mode = ''custom''.');
    end
    % Warn if alpha_0 was also explicitly non-default
    if abs(alpha_0_input) > 1e-12
        warning('init_coupled_params:ConflictingOrientation', ...
            ['Both initial_orientation and alpha_0 (%.4f rad) are set; ', ...
             'initial_orientation takes precedence and alpha_0 is ignored.'], ...
            alpha_0_input);
    end
    q_0_init = parse_initial_orientation(init_orient_input);
    orientation_source = 'initial_orientation';
    params.alpha_0 = NaN;    % not applicable under explicit orientation
else
    % Legacy alpha_0 path (validated to [-pi/2, pi/2])
    if alpha_0_input < -pi/2 || alpha_0_input > pi/2
        error('init_coupled_params:InvalidAlpha0', ...
            'alpha_0 = %.4f rad is outside valid range [-pi/2, pi/2].', alpha_0_input);
    end
    params.alpha_0 = alpha_0_input;

    % Custom-mode warning: alpha_0 inherits the parametric convention
    % (rotate body frame about +Y so that body-X tilts by alpha_0 toward
    % +Z).  This only corresponds to physical long-axis tilt if the user
    % point cloud has its longest extent along body-X.  Otherwise,
    % initial_orientation should be used instead (m16).
    if strcmp(geometry_mode, 'custom') && abs(alpha_0_input) > 1e-12
        warning('init_coupled_params:CustomAlpha0Assumption', ...
            ['Custom mode: alpha_0 = %.4f rad rotates the body about +Y ', ...
             'assuming the long axis is body-X. For arbitrary custom ', ...
             'bodies, prefer specifying custom_geometry.initial_orientation ', ...
             '(quaternion / axis_angle / euler_zyx).'], alpha_0_input);
    end

    % alpha_0: rotation about +Y by theta_Y = -alpha_0
    theta_Y  = -alpha_0_input;
    q_0_init = [cos(theta_Y/2); 0; sin(theta_Y/2); 0];
    q_0_init = q_0_init / norm(q_0_init);
end

params.q_0_init            = q_0_init;
params.orientation_source  = orientation_source;

%% Compute body surface characteristic spacing for RSM regularisation
X_body_local = params.body_local';                  % [N_body x 3]
D2_body = sum((reshape(X_body_local, [N_body_actual, 1, 3]) - ...
               reshape(X_body_local, [1, N_body_actual, 3])).^2, 3);
D2_body(1:N_body_actual+1:end) = Inf;              % Exclude self
min_dists_body = sqrt(min(D2_body, [], 2));

h_min  = min(min_dists_body);
h_mean = mean(min_dists_body);

if isempty(p.Results.eps_body)
    % Default from batch_sweep_N_body_eps_RSM (Stokes sphere): eps/ds_avg ~ 0.375
    params.eps_body = h_mean / 2.6667;
else
    params.eps_body = p.Results.eps_body;
end
params.body_h_min  = h_min;
params.body_h_mean = h_mean;

%% Compute adaptive gap distance
eps_pili_dimless = params.epsilon;
eps_body_dim     = params.eps_body * L_0;
eps_pili_dim     = eps_pili_dimless * L_0;
dis              = 1.1 * (eps_body_dim + eps_pili_dim);
params.dim.dis   = dis;

%% Build pili geometry (dispatch by geometry_mode)
switch geometry_mode
    case 'parametric'
        Pili_Matrix = p.Results.Pili_Matrix;
    case 'custom'
        if isfield(custom_cfg, 'pili') && ~isempty(custom_cfg.pili)
            Pili_Matrix = zeros(8, 0);   % Placeholder; custom pili built separately
        else
            Pili_Matrix = zeros(8, 0);   % No pili defined
        end
end

M_pili      = size(Pili_Matrix, 2);
params.M_pili = M_pili;

if M_pili > 0 && strcmp(geometry_mode, 'parametric')
    % Validate: pili lengths must be strictly positive
    pili_lengths_check = Pili_Matrix(2, :);
    zero_mask = (pili_lengths_check <= 0);
    if any(zero_mask)
        zero_idx = find(zero_mask);
        error('init_coupled_params:ZeroPiliLength', ...
            ['Pili [%s] has zero or negative length (%.4g um).\n' ...
             'To simulate a body with no pili, use an empty Pili_Matrix: ' ...
             'Pili_Matrix = zeros(8, 0).'], ...
            num2str(zero_idx), pili_lengths_check(zero_idx(1)));
    end

    % Validate: check for duplicate attachment sites
    attach_sites = Pili_Matrix(8, :);
    [unique_sites, ~, ic] = unique(attach_sites);
    if length(unique_sites) < M_pili
        counts = accumarray(ic, 1);
        dup_mask = counts > 1;
        dup_sites = unique_sites(dup_mask);
        for ds = 1:length(dup_sites)
            dup_idx = find(attach_sites == dup_sites(ds));
            warning('init_coupled_params:DuplicateAttachSite', ...
                'Pili [%s] share the same attachment site %d. This may be physically unrealistic.', ...
                num2str(dup_idx), dup_sites(ds));
        end
    end

    % Generate initial pili geometry (to be resampled below)
    Pili_struct = pili_building(Pili_Matrix, p.Results.seg_per_um, ...
        p.Results.major_axis, p.Results.minor_axis, dis, ...
        p.Results.E, p.Results.r_pili);

    % Handle N_segments (scalar or vector)
    N_segments_input = p.Results.N_segments;
    if isscalar(N_segments_input)
        N_segments = repmat(N_segments_input, M_pili, 1);
    else
        N_segments = N_segments_input(:);
    end
    params.N_segments = N_segments;

    % Allocate pili storage
    params.r_local  = cell(M_pili, 1);
    params.d_local  = cell(M_pili, 1);
    params.D_local  = cell(M_pili, 1);
    params.ell      = zeros(M_pili, 1);
    params.rot      = cell(M_pili, 1);
    params.kb0      = cell(M_pili, 1);
    params.pili_x0  = cell(M_pili, 1);

    for k = 1:M_pili
        % Root position and direction (nondimensionalised)
        params.r_local{k} = Pili_struct(k).root_local / L_0;
        params.d_local{k} = Pili_struct(k).root_dir;

        % Dimensionless pili length
        params.ell(k) = Pili_struct(k).L / L_0;

        % Initial adaptive basis: identity
        params.rot{k} = eye(3);

        % Initial pili node coordinates (nondimensionalised, lab frame)
        x_pili_dim     = Pili_struct(k).x;
        x_pili_dimless = x_pili_dim / L_0;

        % Resample to N_segments(k)+1 nodes if needed
        N_orig   = size(x_pili_dimless, 2);
        N_target = N_segments(k) + 1;

        if N_orig ~= N_target
            s_orig   = [0, cumsum(sqrt(sum(diff(x_pili_dimless, 1, 2).^2, 1)))];
            s_target = linspace(0, s_orig(end), N_target);
            x_resampled = zeros(3, N_target);
            for dim = 1:3
                x_resampled(dim,:) = interp1(s_orig, x_pili_dimless(dim,:), ...
                                              s_target, 'pchip');
            end
            params.pili_x0{k} = x_resampled';
        else
            params.pili_x0{k} = x_pili_dimless';
        end

        % Compute intrinsic curvature in Euler-angle representation
        coords_k    = params.pili_x0{k};
        tangents_k  = diff(coords_k, 1, 1);
        seg_lens_k  = sqrt(sum(tangents_k.^2, 2));
        tangents_k  = tangents_k ./ seg_lens_k;

        cos_theta0 = max(-1, min(1, tangents_k(:, 3)));
        theta0 = acos(cos_theta0);
        phi0   = atan2(tangents_k(:, 2), tangents_k(:, 1));
        psi0   = zeros(N_segments(k), 1);

        [d1_0, d2_0, d3_0] = pili_directors(theta0, phi0, psi0);

        % Store first-segment director frame (body-local)
        params.D_local{k} = [d1_0(1,:)', d2_0(1,:)', d3_0(1,:)'];

        % Compute intrinsic curvature at nodes.
        % Use the SAME node-director operator (compute_node_directors:
        % SVD-projected segment averaging) and the SAME finite-difference
        % stencil as the run-time curvature in compute_pili_elastic_moments.m.
        % This guarantees kappa_live - kb0 == 0 (bit-exact) in the reference
        % configuration, so a 'frozen' pilus starts as an exact zero-moment
        % equilibrium instead of carrying a spurious restoring moment.
        [d1NC0, d2NC0, d3NC0] = compute_node_directors(d1_0, d2_0, d3_0, N_segments(k));

        d1NCds0 = [d1NC0(2,:)-d1NC0(1,:); ...
                   d1NC0(2:end,:)-d1NC0(1:end-1,:)] * N_segments(k);
        d2NCds0 = [d2NC0(2,:)-d2NC0(1,:); ...
                   d2NC0(2:end,:)-d2NC0(1:end-1,:)] * N_segments(k);
        d3NCds0 = [d3NC0(2,:)-d3NC0(1,:); ...
                   d3NC0(2:end,:)-d3NC0(1:end-1,:)] * N_segments(k);

        kb0_k = zeros(N_segments(k)+1, 3);
        kb0_k(:, 1) = dot(d3NC0, d2NCds0, 2);
        kb0_k(:, 2) = dot(d1NC0, d3NCds0, 2);
        kb0_k(:, 3) = dot(d2NC0, d1NCds0, 2);
        params.kb0{k} = kb0_k;
    end

    %% Apply intrinsic curvature overrides (per-pili, static)
    ic_input = p.Results.intrinsic_curvature;
    for k = 1:M_pili
        if k <= numel(ic_input) && ~isempty(ic_input{k})
            ic_k = ic_input{k};
            if ischar(ic_k) || isstring(ic_k)
                ic_k = struct('type', char(ic_k));
            end
        else
            ic_k = struct('type', 'frozen');
        end

        switch ic_k.type
            case 'frozen'
                % Keep kb0{k} computed from initial geometry (no change)
            case 'zero'
                params.kb0{k} = zeros(N_segments(k)+1, 3);
            case 'constant'
                kv = ic_k.kappa(:)';
                params.kb0{k} = repmat(kv, N_segments(k)+1, 1);
            case 'helix'
                s_nodes = linspace(0, 1, N_segments(k)+1)';
                hk1 = 0 * s_nodes;
                hk2 = 0 * s_nodes + ic_k.kappa_d2;
                hk3 = 0 * s_nodes + ic_k.kappa_d3;
                params.kb0{k} = [hk1, hk2, hk3];
            otherwise
                error('init_coupled_params:BadIC', ...
                    'Unknown intrinsic_curvature type ''%s'' for pili %d.', ic_k.type, k);
        end
    end

    %% Configure active moments (per-pili, time-varying)
    am_input = p.Results.active_moment;
    params.am_cfg = cell(M_pili, 1);
    for k = 1:M_pili
        if k <= numel(am_input) && ~isempty(am_input{k})
            am_k = am_input{k};
            if ischar(am_k) || isstring(am_k)
                am_k = struct('type', char(am_k));
            end
        else
            am_k = struct('type', 'off');
        end

        switch am_k.type
            case 'off'
                params.am_cfg{k} = am_k;

            case 'traveling_wave'
                L_pili_phys = params.ell(k) * L_0;
                am_k.amplitude_dimless = am_k.amplitude_Nm * L_pili_phys / EI;
                am_k.omega_dimless     = am_k.omega_Hz * 2 * pi * T_0;
                if ~isfield(am_k, 'wavenumber'), am_k.wavenumber = 2*pi; end
                if ~isfield(am_k, 'd1_amp'), am_k.d1_amp = 1; end
                if ~isfield(am_k, 'd2_amp'), am_k.d2_amp = 0; end
                if ~isfield(am_k, 'd3_amp'), am_k.d3_amp = 0; end
                if ~isfield(am_k, 'phase_rad'), am_k.phase_rad = 0; end
                params.am_cfg{k} = am_k;

            case 'custom'
                params.am_cfg{k} = am_k;

            otherwise
                error('init_coupled_params:BadAM', ...
                    'Unknown active_moment type ''%s'' for pili %d.', am_k.type, k);
        end
    end

    %% Configure root coupling and motor torque (per-pili)
    params.root_coupling = cell(M_pili, 1);
    params.motor_cfg    = cell(M_pili, 1);
    rc_input  = p.Results.root_coupling;
    mt_input  = p.Results.motor_torque;
    for k = 1:M_pili
        % Root coupling type
        if k <= numel(rc_input) && ~isempty(rc_input{k})
            params.root_coupling{k} = char(rc_input{k});
        else
            params.root_coupling{k} = 'clamped';
        end

        % Motor torque configuration
        mt_k = struct('type', 'off');
        if k <= numel(mt_input) && ~isempty(mt_input{k})
            mt_k = mt_input{k};
        end

        switch mt_k.type
            case 'off'
                params.motor_cfg{k} = mt_k;

            case 'constant'
                % Nondimensionalise the motor torque.  The motor torque is
                % injected into the segment-1 d3 (psi_1) twist moment-balance
                % row in ode_coupled_rhs.m, i.e. the SAME row type that carries
                % the elastic / active internal moments.  Because curvature and
                % twist are discretised on the per-pilus normalised arc length
                % s' in [0,1] (ds = 1/N_k), that row's moment scale is
                % EI / L_pili (NOT EI / L_0).  The torque must therefore use the
                % pilus length L_pili = ell * L_0, consistent with the active-
                % moment nondimensionalisation above.
                %   torque_pN_nm * 1e-21 gives N*m
                L_pili_phys = params.ell(k) * L_0;
                mt_k.torque_dimless = mt_k.torque_pN_nm * 1e-21 / (EI / L_pili_phys);
                % Nondimensionalise: max_omega_dimless = max_omega [Hz] * 2*pi * T_0
                mt_k.max_omega_dimless = mt_k.max_omega_Hz * 2 * pi * T_0;
                params.motor_cfg{k} = mt_k;

            otherwise
                error('init_coupled_params:BadMT', ...
                    'Unknown motor_torque type ''%s'' for pili %d.', mt_k.type, k);
        end
    end

else
    % --- Custom pili geometry ---
    if strcmp(geometry_mode, 'custom') && isfield(custom_cfg, 'pili') && ~isempty(custom_cfg.pili)
        custom_pili_cfg = custom_cfg.pili;
        M_pili = numel(custom_pili_cfg);
        params.M_pili = M_pili;

        % Build attachment sites from custom_body_building output
        attachment_sites = params.custom_attachment_sites;

        % Generate initial custom pili geometry
        Pili_struct = custom_pili_building(custom_pili_cfg, attachment_sites, ...
            p.Results.seg_per_um, L_0, dis, ...
            p.Results.E, p.Results.r_pili);

        % Compute N_segments from actual pili lengths
        N_segments = zeros(M_pili, 1);
        for k = 1:M_pili
            N_segments(k) = Pili_struct(k).N - 1;
        end
        params.N_segments = N_segments;

        % Allocate pili storage
        params.r_local  = cell(M_pili, 1);
        params.d_local  = cell(M_pili, 1);
        params.D_local  = cell(M_pili, 1);
        params.ell      = zeros(M_pili, 1);
        params.rot      = cell(M_pili, 1);
        params.kb0      = cell(M_pili, 1);
        params.pili_x0  = cell(M_pili, 1);

        for k = 1:M_pili
            % Root position and direction (nondimensionalised)
            params.r_local{k} = Pili_struct(k).root_local / L_0;
            params.d_local{k} = Pili_struct(k).root_dir;

            % Dimensionless pili length
            params.ell(k) = Pili_struct(k).L / L_0;

            % Initial adaptive basis: identity
            params.rot{k} = eye(3);

            % Pili node coordinates (nondimensionalised, lab frame)
            x_pili_dimless = Pili_struct(k).x / L_0;
            params.pili_x0{k} = x_pili_dimless';

            % Compute intrinsic curvature in Euler-angle representation
            coords_k    = params.pili_x0{k};
            tangents_k  = diff(coords_k, 1, 1);
            seg_lens_k  = sqrt(sum(tangents_k.^2, 2));
            tangents_k  = tangents_k ./ seg_lens_k;

            cos_theta0 = max(-1, min(1, tangents_k(:, 3)));
            theta0 = acos(cos_theta0);
            phi0   = atan2(tangents_k(:, 2), tangents_k(:, 1));
            % Initial twist psi0 = 0 for custom pili.
            % Rationale: the user-supplied centerline defines the tangent
            % field t(s), but there is no unique natural binormal for an
            % arbitrary 3D curve (Frenet frames fail for straight segments,
            % and parallel-transport frames require a seed).  We therefore
            % adopt psi0 = 0 in the spherical-angle (theta, phi, psi)
            % parameterisation as a canonical starting point, consistent
            % with the parametric branch.  If the user needs an intrinsic
            % twist, they should use intrinsic_curvature = 'constant' /
            % 'helix' (see parse_pili_config and USER_GUIDE).
            psi0   = zeros(N_segments(k), 1);

            [d1_0, d2_0, d3_0] = pili_directors(theta0, phi0, psi0);

            % Store first-segment director frame (body-local)
            params.D_local{k} = [d1_0(1,:)', d2_0(1,:)', d3_0(1,:)'];

            % Compute intrinsic curvature at nodes.
            % Use the SAME node-director operator (compute_node_directors:
            % SVD-projected segment averaging) and the SAME finite-difference
            % stencil as the run-time curvature in compute_pili_elastic_moments.m.
            % This guarantees kappa_live - kb0 == 0 (bit-exact) in the reference
            % configuration, so a 'frozen' pilus starts as an exact zero-moment
            % equilibrium instead of carrying a spurious restoring moment.
            [d1NC0, d2NC0, d3NC0] = compute_node_directors(d1_0, d2_0, d3_0, N_segments(k));

            d1NCds0 = [d1NC0(2,:)-d1NC0(1,:); ...
                       d1NC0(2:end,:)-d1NC0(1:end-1,:)] * N_segments(k);
            d2NCds0 = [d2NC0(2,:)-d2NC0(1,:); ...
                       d2NC0(2:end,:)-d2NC0(1:end-1,:)] * N_segments(k);
            d3NCds0 = [d3NC0(2,:)-d3NC0(1,:); ...
                       d3NC0(2:end,:)-d3NC0(1:end-1,:)] * N_segments(k);

            kb0_k = zeros(N_segments(k)+1, 3);
            kb0_k(:, 1) = dot(d3NC0, d2NCds0, 2);
            kb0_k(:, 2) = dot(d1NC0, d3NCds0, 2);
            kb0_k(:, 3) = dot(d2NC0, d1NCds0, 2);
            params.kb0{k} = kb0_k;
        end

        %% Apply intrinsic curvature overrides (per-pili, static) — custom mode
        ic_input = p.Results.intrinsic_curvature;
        for k = 1:M_pili
            if k <= numel(ic_input) && ~isempty(ic_input{k})
                ic_k = ic_input{k};
                if ischar(ic_k) || isstring(ic_k)
                    ic_k = struct('type', char(ic_k));
                end
            else
                ic_k = struct('type', 'frozen');
            end

            switch ic_k.type
                case 'frozen'
                    % Keep kb0{k} computed from initial geometry
                case 'zero'
                    params.kb0{k} = zeros(N_segments(k)+1, 3);
                case 'constant'
                    kv = ic_k.kappa(:)';
                    params.kb0{k} = repmat(kv, N_segments(k)+1, 1);
                case 'helix'
                    s_nodes = linspace(0, 1, N_segments(k)+1)';
                    hk1 = 0 * s_nodes;
                    hk2 = 0 * s_nodes + ic_k.kappa_d2;
                    hk3 = 0 * s_nodes + ic_k.kappa_d3;
                    params.kb0{k} = [hk1, hk2, hk3];
                otherwise
                    error('init_coupled_params:BadIC', ...
                        'Unknown intrinsic_curvature type ''%s'' for pili %d.', ic_k.type, k);
            end
        end

        %% Configure active moments (per-pili, time-varying) — custom mode
        am_input = p.Results.active_moment;
        params.am_cfg = cell(M_pili, 1);
        for k = 1:M_pili
            if k <= numel(am_input) && ~isempty(am_input{k})
                am_k = am_input{k};
                if ischar(am_k) || isstring(am_k)
                    am_k = struct('type', char(am_k));
                end
            else
                am_k = struct('type', 'off');
            end

            switch am_k.type
                case 'off'
                    params.am_cfg{k} = am_k;

                case 'traveling_wave'
                    L_pili_phys = params.ell(k) * L_0;
                    am_k.amplitude_dimless = am_k.amplitude_Nm * L_pili_phys / EI;
                    am_k.omega_dimless     = am_k.omega_Hz * 2 * pi * T_0;
                    if ~isfield(am_k, 'wavenumber'), am_k.wavenumber = 2*pi; end
                    if ~isfield(am_k, 'd1_amp'), am_k.d1_amp = 1; end
                    if ~isfield(am_k, 'd2_amp'), am_k.d2_amp = 0; end
                    if ~isfield(am_k, 'd3_amp'), am_k.d3_amp = 0; end
                    if ~isfield(am_k, 'phase_rad'), am_k.phase_rad = 0; end
                    params.am_cfg{k} = am_k;

                case 'custom'
                    params.am_cfg{k} = am_k;

                otherwise
                    error('init_coupled_params:BadAM', ...
                        'Unknown active_moment type ''%s'' for pili %d.', am_k.type, k);
            end
        end

        %% Configure root coupling and motor torque — custom mode
        params.root_coupling = cell(M_pili, 1);
        params.motor_cfg    = cell(M_pili, 1);
        rc_input  = p.Results.root_coupling;
        mt_input  = p.Results.motor_torque;
        for k = 1:M_pili
            % Root coupling type
            if k <= numel(rc_input) && ~isempty(rc_input{k})
                params.root_coupling{k} = char(rc_input{k});
            else
                params.root_coupling{k} = 'clamped';
            end

            % Motor torque configuration
            mt_k = struct('type', 'off');
            if k <= numel(mt_input) && ~isempty(mt_input{k})
                mt_k = mt_input{k};
            end

            switch mt_k.type
                case 'off'
                    params.motor_cfg{k} = mt_k;

                case 'constant'
                    % Motor torque enters the segment-1 psi_1 twist moment-
                    % balance row, whose scale is EI / L_pili (per-pilus s' in
                    % [0,1] arc-length normalisation), matching the active-
                    % moment convention.  Use L_pili = ell * L_0, not L_0.
                    L_pili_phys = params.ell(k) * L_0;
                    mt_k.torque_dimless = mt_k.torque_pN_nm * 1e-21 / (EI / L_pili_phys);
                    mt_k.max_omega_dimless = mt_k.max_omega_Hz * 2 * pi * T_0;
                    params.motor_cfg{k} = mt_k;

                otherwise
                    error('init_coupled_params:BadMT', ...
                        'Unknown motor_torque type ''%s'' for pili %d.', mt_k.type, k);
            end
        end

    else
        % No pili (either mode)
        params.N_segments = [];
        params.r_local    = {};
        params.d_local    = {};
        params.D_local    = {};
        params.ell        = [];
        params.rot        = {};
        params.kb0        = {};
        params.pili_x0    = {};
        params.am_cfg     = {};
        params.root_coupling = {};
        params.motor_cfg  = {};
        ic_input          = {};
    end
end

%% Precompute trapezoidal quadrature weights for the global point cloud
%  quad_weights(i) = 1 for body points, and trapezoidal weight ds/2 or ds
%  for pili node points.  Used by build_coupled_lin_sys to vectorise the
%  angular-moment rows of B_global (replaces get_trapezoidal_weight).
N_total_pts = params.N_body + sum(params.N_segments) + params.M_pili;
params.quad_weights = ones(N_total_pts, 1);
off_qw = params.N_body;
for k = 1:params.M_pili
    N_k  = params.N_segments(k);
    ds_k = params.ell(k) / N_k;
    w_k  = ds_k * ones(N_k + 1, 1);
    w_k(1)   = ds_k / 2;
    w_k(end) = ds_k / 2;
    params.quad_weights(off_qw + 1 : off_qw + N_k + 1) = w_k;
    off_qw = off_qw + N_k + 1;
end

%% Configure background flow (global, steady)
bg_input = p.Results.background_flow;
if ~isfield(bg_input, 'type'), bg_input.type = 'off'; end

params.bg_flow.type = bg_input.type;
if ~strcmp(bg_input.type, 'off')
    U_max_phys = bg_input.U_max;
    params.bg_flow.U_max_dimless = U_max_phys / U_0;

    if isfield(bg_input, 'direction')
        d_vec = bg_input.direction(:);
    else
        d_vec = [1; 0; 0];
    end
    params.bg_flow.direction = d_vec / max(norm(d_vec), 1e-30);

    if isfield(bg_input, 'shear_axis')
        s_vec = bg_input.shear_axis(:);
    else
        s_vec = [0; 0; 1];
    end
    params.bg_flow.shear_axis = s_vec / max(norm(s_vec), 1e-30);

    if isfield(bg_input, 'shear_height_m')
        params.bg_flow.shear_height_dimless = bg_input.shear_height_m / L_0;
    else
        params.bg_flow.shear_height_dimless = 1.0;
    end
else
    params.bg_flow.U_max_dimless        = 0;
    params.bg_flow.direction            = [1; 0; 0];
    params.bg_flow.shear_axis           = [0; 0; 1];
    params.bg_flow.shear_height_dimless = 1.0;
end

%% Brownian motion configuration
%  The effective body diffusion tensor is now computed dynamically at each
%  Brownian step via compute_effective_body_brownian, which probes the full
%  coupled body-pili linear system.  No pre-computation is needed here.
%  The old body-only routine compute_body_brownian_params.m has been
%  removed; its functionality is fully superseded by the new method.
%  Thermal drift is handled explicitly in run_simulation.m with the strict
%  body-only Ito finite-difference correction kB*T*div(M).

%% Run token for cross-step persistent caches in ode_coupled_rhs
%  ode_coupled_rhs.m keeps a persistent dpsi_1 cache for motor speed-
%  limiting (free_torsion + constant motor torque).  The cache must be
%  reset between independent simulations sharing the same MATLAB session;
%  comparing this token lets that function detect a fresh `params` and
%  flush the cache automatically.  `tic` is a high-resolution monotonic
%  counter, so collisions across runs are practically impossible.
params.run_token = uint64(tic);

%% Per-pilus dpsi_1 cache (placeholder, see ode_coupled_rhs persistent block)
%  Initialised here so that downstream code can reference the field if
%  needed; the authoritative cross-step state lives in the persistent
%  variable inside ode_coupled_rhs.m.
if M_pili > 0
    params.dpsi1_prev = cell(M_pili, 1);
else
    params.dpsi1_prev = {};
end

end
