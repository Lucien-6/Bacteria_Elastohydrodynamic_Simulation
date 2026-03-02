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
%    'N_body'      Number of body surface points      (default: 1000)
%    'minor_axis'  Semi-minor axis [m]                 (default: 0.4e-6)
%    'major_axis'  Semi-major axis [m]                 (default: 1.25e-6)
%    'mu'          Dynamic viscosity [Pa*s]            (default: 1e-3)
%    'E'           Young's modulus [Pa]                (default: 3.67e8)
%    'r_pili'      Pili radius [m]                     (default: 3e-9)
%    'sigma'       Poisson's ratio                     (default: 0.2)
%    'rho_bac'     Bacteria density [kg/m^3]           (default: 1180)
%    'rho_fluid'   Fluid density [kg/m^3]              (default: 998.232)
%    'N_segments'  Segment counts [M x 1]              (default: 300)
%    'Pili_Matrix' Pili parameter matrix [8 x M]       (default: single pole)
%    'seg_per_um'  Segments per micron                 (default: 15)
%    'delta'       Euler-angle singularity threshold   (default: pi/20)
%    'alpha_0'     Initial body tilt angle [rad]       (default: pi/2)
%                  Angle of body long axis (body X) from +X in XZ plane
%                  pi/2 = vertical (+Z), 0 = horizontal (+X),
%                  -pi/2 = inverted (-Z).  Range: [-pi/2, pi/2]
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
%  Updated: 2026-02-14 — Pure RSS mode (removed RSM fallback)
%           2026-03-01 — eps_body: h_min/2 → h_mean/6 (Cortez 2005 optimal ratio)
%           2026-03-01 — eps_body: h_mean/6 → h_mean/3 (empirically optimal)

%% Parse inputs
p = inputParser;
addParameter(p, 'N_body',      1000,   @isnumeric);
addParameter(p, 'minor_axis',  0.4e-6, @isnumeric);
addParameter(p, 'major_axis',  1.25e-6, @isnumeric);
addParameter(p, 'mu',          1e-3,   @isnumeric);
addParameter(p, 'E',           3.67e8, @isnumeric);
addParameter(p, 'r_pili',      3e-9,   @isnumeric);
addParameter(p, 'sigma',       0.2,    @isnumeric);
addParameter(p, 'rho_bac',     1180.0, @isnumeric);
addParameter(p, 'rho_fluid',   998.232, @isnumeric);
addParameter(p, 'N_segments',  300,    @isnumeric);
addParameter(p, 'Pili_Matrix', [1; 3; 0; 0; NaN; NaN; NaN; 1], @isnumeric);
addParameter(p, 'seg_per_um',  15,     @isnumeric);
addParameter(p, 'delta',       pi/20,  @isnumeric);
addParameter(p, 'alpha_0',    pi/2,   @isnumeric);
parse(p, varargin{:});

%% Store dimensional parameters
params.dim.minor_axis = p.Results.minor_axis;
params.dim.major_axis = p.Results.major_axis;
params.dim.mu         = p.Results.mu;
params.dim.E          = p.Results.E;
params.dim.r_pili     = p.Results.r_pili;
params.dim.rho_bac    = p.Results.rho_bac;
params.dim.rho_fluid  = p.Results.rho_fluid;

%% Characteristic scales for nondimensionalisation
L_0 = p.Results.major_axis;                       % Length scale [m]
I_pili = pi * p.Results.r_pili^4 / 4;             % Second moment of area [m^4]
EI  = p.Results.E * I_pili;                        % Bending stiffness [N*m^2]
T_0 = 8 * pi * p.Results.mu * L_0^4 / EI;         % Time scale [s]
U_0 = L_0 / T_0;                                   % Velocity scale [m/s]
F_0 = EI / L_0^2;                                  % Force scale [N]

params.scales.L_0 = L_0;
params.scales.T_0 = T_0;
params.scales.U_0 = U_0;
params.scales.F_0 = F_0;
params.scales.EI  = EI;

%% Compute sinking force and dimensionless parameters
% Capsule (spherocylinder) volume:
%   V = (4/3)*pi*b^3 + pi*b^2*2*(a-b)
a = p.Results.major_axis;
b = p.Results.minor_axis;
V_capsule = pi * b^2 * (4/3*b + 2*(a - b));

Delta_rho = p.Results.rho_bac - p.Results.rho_fluid;
g         = 9.81;
F_sink    = Delta_rho * V_capsule * g;             % Sinking force [N]

params.Eg      = F_sink * L_0^2 / EI;             % Elastogravitational number
params.epsilon = p.Results.r_pili / L_0;           % Global slenderness
params.sigma   = p.Results.sigma;
params.g_hat   = [0; 0; -1];                       % Gravity direction
params.delta   = p.Results.delta;

% Initial body tilt angle (validated to [-pi/2, pi/2])
alpha_0_input = p.Results.alpha_0;
if alpha_0_input < -pi/2 || alpha_0_input > pi/2
    error('init_coupled_params:InvalidAlpha0', ...
        'alpha_0 = %.4f rad is outside valid range [-pi/2, pi/2].', alpha_0_input);
end
params.alpha_0 = alpha_0_input;

%% Build capsule body geometry
N_body_requested = p.Results.N_body;
[Capsule_str, N_body_actual] = capsule_body_building( ...
    N_body_requested, p.Results.minor_axis, p.Results.major_axis);

% Store body surface points in LOCAL frame (nondimensionalised)
params.N_body    = N_body_actual;
params.body_local = [Capsule_str.rxB'; Capsule_str.ryB'; Capsule_str.rzB'] / L_0;

%% Compute body surface characteristic spacing for RSM regularisation
X_body_local = params.body_local';                  % [N_body x 3]
D2_body = sum((reshape(X_body_local, [N_body_actual, 1, 3]) - ...
               reshape(X_body_local, [1, N_body_actual, 3])).^2, 3);
D2_body(1:N_body_actual+1:end) = Inf;              % Exclude self
min_dists_body = sqrt(min(D2_body, [], 2));

h_min  = min(min_dists_body);
h_mean = mean(min_dists_body);

params.eps_body   = h_mean / 3;                    % Empirically optimal (eps/ds_avg ~ 0.333)
params.body_h_min  = h_min;
params.body_h_mean = h_mean;

%% Compute adaptive gap distance
eps_pili_dimless = params.epsilon;
eps_body_dim     = params.eps_body * L_0;
eps_pili_dim     = eps_pili_dimless * L_0;
dis              = 1.1 * (eps_body_dim + eps_pili_dim);
params.dim.dis   = dis;

fprintf('Computed adaptive gap distance:\n');
fprintf('  eps_body = %.2e m, eps_pili = %.2e m\n', eps_body_dim, eps_pili_dim);
fprintf('  dis = 1.1 * (eps_body + eps_pili) = %.2e m\n\n', dis);

%% Build pili geometry
Pili_Matrix = p.Results.Pili_Matrix;
M_pili      = size(Pili_Matrix, 2);
params.M_pili = M_pili;

if M_pili > 0
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

        % Compute intrinsic curvature at nodes
        d1NC0 = [d1_0; d1_0(end,:)];
        d2NC0 = [d2_0; d2_0(end,:)];
        d3NC0 = [d3_0; d3_0(end,:)];

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
else
    params.N_segments = [];
    params.r_local    = {};
    params.d_local    = {};
    params.D_local    = {};
    params.ell        = [];
    params.rot        = {};
    params.kb0        = {};
    params.pili_x0    = {};
end

%% Print summary
fprintf('\n');
fprintf('=============================================================\n');
fprintf('  Coupled Body-Pili Simulation Parameters\n');
fprintf('=============================================================\n');
fprintf('  Body : N_body = %d, a = %.2f um, b = %.2f um\n', ...
    params.N_body, p.Results.major_axis*1e6, p.Results.minor_axis*1e6);
fprintf('  Pili : M = %d', M_pili);
if M_pili > 0
    fprintf(', N_seg = [');
    fprintf('%d ', params.N_segments);
    fprintf('], ell = [');
    fprintf('%.2f ', params.ell);
    fprintf(']');
end
fprintf('\n');
fprintf('  Scales : L0 = %.2e m, T0 = %.2e s, EI = %.2e N*m^2\n', L_0, T_0, EI);
fprintf('  Eg = %.4f, eps_slender = %.2e, sigma = %.1f\n', ...
    params.Eg, params.epsilon, params.sigma);
fprintf('  alpha_0 = %.4f rad (%.1f deg) — initial body tilt\n', ...
    params.alpha_0, rad2deg(params.alpha_0));
fprintf('  F_sink = %.2e N, F_0 = %.2e N\n', F_sink, F_0);
fprintf('  V_capsule = %.4e m^3\n', V_capsule);
fprintf('  State vector size : %d DOFs\n', 7 + 3*sum(params.N_segments));
fprintf('-------------------------------------------------------------\n');
fprintf('  Mobility Method : RSM (body) + RSS (pili) + Rotlet  [Cortez 2018]\n');
fprintf('  Regularisation Parameters\n');
fprintf('  Body spacing : h_min = %.4f L0, h_mean = %.4f L0\n', ...
    params.body_h_min, params.body_h_mean);
fprintf('  Body  : eps_body = %.4f L0 = %.1f nm (= h_mean/3, eps/ds ~ 0.333)\n', ...
    params.eps_body, params.eps_body * L_0 * 1e9);
if M_pili > 0
    eps_pili_dim_print = params.epsilon * L_0;
    for k = 1:M_pili
        ds_k    = params.ell(k) / params.N_segments(k);
        ds_nm   = ds_k * L_0 * 1e9;
        eps_nm  = eps_pili_dim_print * 1e9;
        fprintf('  Pili %d: N_seg = %d, ds = %.4f L0 (%.2f nm)\n', ...
            k, params.N_segments(k), ds_k, ds_nm);
        fprintf('          eps_pili = %.4f L0 (%.1f nm = r_pili), ds/eps = %.1f\n', ...
            params.epsilon, eps_nm, ds_nm / eps_nm);
    end
    fprintf('  Gap : dis = 1.1*(eps_body + eps_pili) = %.1f nm\n', dis * 1e9);
    fprintf('  Note: eps_ij = eps_j (determined by force point j)\n');
    fprintf('    Body force -> any point: eps = %.1f nm\n', ...
        params.eps_body * L_0 * 1e9);
    fprintf('    Pili force -> any point: eps = %.1f nm (physical radius)\n', ...
        eps_pili_dim_print * 1e9);
end
fprintf('=============================================================\n\n');

end
