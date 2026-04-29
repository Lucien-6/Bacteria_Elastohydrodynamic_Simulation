function Pili = custom_pili_building(custom_pili_cfg, attachment_sites, seg_per_um, L_0, dis, E, r_pili)
%CUSTOM_PILI_BUILDING  Initialize custom pili from user-provided centerline files
%
%  Description:
%    Constructs pili structures from user-provided centerline point files,
%    applying the same output format as pili_building.m for downstream
%    compatibility.  Each pilus is resampled to the required segment density,
%    aligned to its attachment site direction, translated into the body-local
%    frame, and equipped with bending stiffness and intrinsic curvature.
%
%  Syntax:
%    Pili = custom_pili_building(custom_pili_cfg, attachment_sites, ...
%             seg_per_um, L_0, dis, E, r_pili)
%
%  Inputs:
%    custom_pili_cfg   - [M x 1] struct array, one per pilus:
%        .points_file      - Path to CSV/MAT file with centerline [m]
%        .site_index       - Index into attachment_sites struct array
%        .deflection_z     - Additional rotation about Z-axis [rad] (default: 0)
%        .deflection_x     - Additional rotation about X-axis [rad] (default: 0)
%        .bending_stiffness - Custom A [N*m^2] (optional, default: E*I)
%    attachment_sites  - [S x 1] struct array from custom_body_building:
%        .position         - [3 x 1] Root position [m]
%        .direction        - [3 x 1] Unit outward direction
%        .label            - String label
%    seg_per_um        - Segments per micron for discretisation [seg/um]
%    L_0               - Characteristic length scale [m]
%    dis               - Gap between body surface and pili root [m]
%    E                 - Young's modulus [Pa]
%    r_pili            - Pili radius [m]
%
%  Outputs:
%    Pili - [1 x M] Structure array with fields (same as pili_building output):
%        .N          Number of vertices
%        .L          Total length [m]
%        .ell        Natural segment length [m]
%        .x          Vertex positions [3 x N] [m]
%        .v          Vertex velocities [3 x N] [m/s] (zero)
%        .kb         Placeholder (empty); curvature recomputed by init_coupled_params
%        .kb0        Placeholder (empty); intrinsic curvature recomputed by init_coupled_params
%        .A          Bending stiffness [N*m^2]
%        .r          Pili radius [m]
%        .f_bend     Bending forces [3 x N] [N] (zero)
%        .attach_idx Attachment site index (1-based)
%        .root_local Root position in body-local coordinates [3 x 1] [m]
%        .root_dir   Root tangent direction in body-local coords [3 x 1]
%        .type       Morphology type (0 = custom)
%        .is_clamped Clamped boundary condition flag (default: true)
%
%  See also: pili_building, read_geometry_csv_mat, rotation_between_vectors
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-28

%% Input validation
if nargin < 7
    error('custom_pili_building:InsufficientInputs', ...
          'custom_pili_building requires at least 7 input arguments.');
end

%% Calculate default bending stiffness (EI = E * pi * r^4 / 4)
I_pili    = pi * r_pili^4 / 4;   % Second moment of area [m^4]
A_default = E * I_pili;           % Bending stiffness [N*m^2]

%% Get number of pili
NPili = numel(custom_pili_cfg);

%% Pre-allocate output structure array (same fields as pili_building)
Pili(NPili) = struct('N', [], 'L', [], 'ell', [], ...
    'x', [], 'v', [], 'kb', [], 'kb0', [], ...
    'A', [], 'r', [], 'f_bend', [], ...
    'attach_idx', [], 'root_local', [], 'root_dir', [], ...
    'type', [], 'is_clamped', []);

%% Construct each pilus
for j = 1:NPili

    % --- Extract configuration ---
    cfg_j = custom_pili_cfg(j);
    site_idx = cfg_j.site_index;

    % Deflection angles (additional rotation after alignment)
    deflection_z = 0;
    if isfield(cfg_j, 'deflection_z') && ~isempty(cfg_j.deflection_z)
        deflection_z = cfg_j.deflection_z;
    end
    deflection_x = 0;
    if isfield(cfg_j, 'deflection_x') && ~isempty(cfg_j.deflection_x)
        deflection_x = cfg_j.deflection_x;
    end

    % Bending stiffness
    if isfield(cfg_j, 'bending_stiffness') && ~isempty(cfg_j.bending_stiffness)
        A = cfg_j.bending_stiffness;
    else
        A = A_default;
    end

    % --- Read centerline from file ---
    raw_pts = read_geometry_csv_mat(cfg_j.points_file);   % [N_raw x 3] [m]

    % --- Compute actual arc length ---
    seg_vecs = diff(raw_pts, 1, 1);                       % [(N_raw-1) x 3]
    seg_lens = sqrt(sum(seg_vecs.^2, 2));                  % [(N_raw-1) x 1]
    L_actual = sum(seg_lens);                              % Total length [m]

    if L_actual <= 0
        error('custom_pili_building:ZeroLength', ...
              'Pilus %d: computed arc length is zero or negative (%.4e m).', j, L_actual);
    end

    % --- Determine segment count ---
    L_um  = L_actual * 1e6;                               % Length [um]
    N     = max(ceil(L_um * seg_per_um) + 1, 3);          % At least 3 nodes
    ell   = L_actual / (N - 1);                            % Segment length [m]

    % --- Resample to uniform arc length using pchip interpolation ---
    s_raw = [0; cumsum(seg_lens)];                         % [N_raw x 1]
    s_new = linspace(0, L_actual, N)';                     % [N x 1]

    x_resampled = zeros(N, 3);
    for dim = 1:3
        x_resampled(:, dim) = interp1(s_raw, raw_pts(:, dim), s_new, 'pchip');
    end

    % --- Compute root tangent direction ---
    root_tangent = x_resampled(2, :) - x_resampled(1, :);  % [1 x 3]
    root_tangent = root_tangent / max(norm(root_tangent), 1e-30);

    % --- Align root tangent to attachment site direction ---
    site = attachment_sites(site_idx);
    site_dir = site.direction(:);                           % [3 x 1]
    root_tangent_col = root_tangent(:)';                    % [1 x 3]

    % Rotation from current root tangent to site direction
    R_align = rotation_between_vectors(root_tangent_col', site_dir);

    % Apply alignment rotation to all resampled points
    x_aligned = (R_align * x_resampled')';                 % [N x 3]

    % --- Apply additional deflection rotations ---
    if abs(deflection_z) > 1e-14 || abs(deflection_x) > 1e-14
        RZ = [ cos(deflection_z), -sin(deflection_z), 0;
               sin(deflection_z),  cos(deflection_z), 0;
               0,                   0,                  1];
        RX = [1,  0,                  0;
              0,  cos(deflection_x),  sin(deflection_x);
              0, -sin(deflection_x),  cos(deflection_x)];
        R_deflect = RX * RZ;
        x_aligned = (R_deflect * x_aligned')';
    end

    % --- Translate to attachment site position + gap ---
    %  Root is placed at: site.position + dis * site.direction
    root_pos = site.position(:) + dis * site.direction(:);  % [3 x 1]
    offset   = root_pos' - x_aligned(1, :);                 % Shift so node 1 = root_pos
    x_final  = x_aligned + repmat(offset, N, 1);            % [N x 3]

    % --- Pack as [3 x N] column-major format (same as pili_building) ---
    x = x_final';                                           % [3 x N]

    % --- Initialize velocities and forces to zero ---
    v      = zeros(3, N);
    f_bend = zeros(3, N);

    % --- Store in structure ---
    Pili(j).N          = N;
    Pili(j).L          = L_actual;
    Pili(j).ell        = ell;
    Pili(j).x          = x;
    Pili(j).v          = v;
    Pili(j).kb         = zeros(3, 0);                        % Not used downstream; init_coupled_params recomputes
    Pili(j).kb0        = zeros(3, 0);                        % Intrinsic curvature recomputed in init_coupled_params
    Pili(j).A          = A;
    Pili(j).r          = r_pili;
    Pili(j).f_bend     = f_bend;
    Pili(j).attach_idx = site_idx;
    Pili(j).root_local = root_pos;                          % [3 x 1]
    Pili(j).root_dir   = site_dir;                          % [3 x 1]
    Pili(j).type       = 0;                                 % 0 = custom type
    Pili(j).is_clamped = true;

    fprintf('  [Custom Pili %d] L = %.3f um, N = %d, site = %d (''%s''), defl = [%.3f, %.3f] rad\n', ...
        j, L_um, N, site_idx, site.label, deflection_z, deflection_x);

end

end
