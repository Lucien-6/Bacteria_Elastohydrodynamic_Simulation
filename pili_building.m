function Pili = pili_building(Pili_Matrix, seg_per_um, major_axis, minor_axis, dis, E, r_pili)
%PILI_BUILDING  Initialize pili structures for discrete elastic rod simulation
%
%  Description:
%    Constructs pili structures based on the pili geometry generator,
%    adding bending stiffness and computing intrinsic curvature binormal
%    vectors for the undeformed configuration.
%
%  Syntax:
%    Pili = pili_building(Pili_Matrix, seg_per_um, major_axis, minor_axis, dis, E, r_pili)
%
%  Inputs:
%    Pili_Matrix  - [8-9 x M] Pili parameter matrix (M = number of pili)
%        Row 1 : Morphology type (1=linear, 2=parabolic, 3=circular,
%                 4=sinusoidal, 5=conical spiral)
%        Row 2 : Length [um]
%        Row 3 : Deflection angle relative to Z-axis [rad]
%        Row 4 : Deflection angle relative to X-axis [rad]
%        Row 5-7: Shape control parameters (type-dependent)
%        Row 8 : Attachment position on bacterial body (1-14)
%        Row 9 : (optional) Custom bending stiffness A [N*m^2]
%    seg_per_um   - Segments per micron for discretization [seg/um]
%    major_axis   - Semi-major axis of bacterial body [m]
%    minor_axis   - Semi-minor axis of bacterial body [m]
%    dis          - Gap between body surface and pili root [m]
%    E            - Young's modulus of pili material [Pa]
%    r_pili       - Radius of pili cross-section [m]
%
%  Outputs:
%    Pili - [1 x M] Structure array with fields:
%        .N          Number of vertices
%        .L          Total length [m]
%        .ell        Natural segment length [m]
%        .x          Vertex positions [3 x N] [m]
%        .v          Vertex velocities [3 x N] [m/s] (zero)
%        .kb         Curvature binormal vectors [3 x (N-2)] [1/m]
%        .kb0        Intrinsic curvature binormal vectors [3 x (N-2)] [1/m]
%        .A          Bending stiffness [N*m^2]
%        .r          Pili radius [m]
%        .f_bend     Bending forces [3 x N] [N] (zero)
%        .attach_idx Attachment point index
%        .root_local Root position in body-local coordinates [3 x 1] [m]
%        .root_dir   Root tangent direction in body-local coords [3 x 1]
%        .type       Morphology type
%        .is_clamped Clamped boundary condition flag (default: true)
%
%  References:
%    [1] Bergou, M. et al. (2008). Discrete Elastic Rods. ACM TOG 27(3):63.
%
%  See also: pili_type_building, compute_curvature_binormal
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-01-23

%% Input validation
if nargin < 7
    error('pili_building:InsufficientInputs', ...
          'pili_building requires at least 7 input arguments.');
end

%% Calculate default bending stiffness (EI = E * pi * r^4 / 4)
I_pili    = pi * r_pili^4 / 4;   % Second moment of area [m^4]
A_default = E * I_pili;           % Bending stiffness [N*m^2]

%% Get number of pili
NPili = size(Pili_Matrix, 2);

%% Pre-allocate output structure array
Pili(NPili) = struct('N', [], 'L', [], 'ell', [], ...
    'x', [], 'v', [], 'kb', [], 'kb0', [], ...
    'A', [], 'r', [], 'f_bend', [], ...
    'attach_idx', [], 'root_local', [], 'root_dir', [], ...
    'type', [], 'is_clamped', []);

%% Construct each pilus
for j = 1:NPili

    % --- Extract parameters ---
    pili_type  = Pili_Matrix(1, j);
    L_um       = Pili_Matrix(2, j);          % Length [um]
    L          = L_um * 1e-6;                % Length [m]
    attach_idx = Pili_Matrix(8, j);

    % Bending stiffness: use custom value if provided in row 9
    if size(Pili_Matrix, 1) >= 9 && ~isnan(Pili_Matrix(9, j))
        A = Pili_Matrix(9, j);
    else
        A = A_default;
    end

    % --- Compute number of vertices ---
    N   = ceil(L_um * seg_per_um) + 1;
    N   = max(N, 3);        % Minimum 3 vertices for curvature calculation
    ell = L / (N - 1);      % Segment length [m]

    % --- Build pilus geometry ---
    Pili_str = pili_type_building(Pili_Matrix(1:7, j), N);

    % --- Apply location transformation ---
    rxP = zeros(1, N);
    ryP = zeros(1, N);
    rzP = zeros(1, N);

    [rxP, ryP, rzP, SPoint, PVector] = pili_location_single(...
        Pili_str, dis, major_axis, minor_axis, rxP, ryP, rzP, attach_idx, N);

    x = [rxP; ryP; rzP];    % [3 x N]

    % --- Compute curvature binormal vectors ---
    kb  = compute_curvature_binormal(x);   % [3 x (N-2)]
    kb0 = kb;                              % Intrinsic (rest) curvature

    % --- Initialize velocities and forces to zero ---
    v      = zeros(3, N);
    f_bend = zeros(3, N);

    % --- Store in structure ---
    Pili(j).N          = N;
    Pili(j).L          = L;
    Pili(j).ell        = ell;
    Pili(j).x          = x;
    Pili(j).v          = v;
    Pili(j).kb         = kb;
    Pili(j).kb0        = kb0;
    Pili(j).A          = A;
    Pili(j).r          = r_pili;
    Pili(j).f_bend     = f_bend;
    Pili(j).attach_idx = attach_idx;
    Pili(j).root_local = SPoint';    % [3 x 1]
    Pili(j).root_dir   = PVector';   % [3 x 1] unit tangent (body-local)
    Pili(j).type       = pili_type;
    Pili(j).is_clamped = true;

end

end


%% ==================== Local Helper Function ====================

function [rxP, ryP, rzP, SPoint, PV] = pili_location_single( ...
        hTail, dis, major_axis, minor_axis, rxP, ryP, rzP, NSite, Ntail)
%PILI_LOCATION_SINGLE  Apply attachment-site transformation for one pilus
%
%  All coordinates are in body-local frame (body centre at origin).
%  The capsule has cylinder spanning x in [-(a-b), +(a-b)] (a=major, b=minor)
%  with hemispherical caps at the two ends.  The "side" attachments
%  (sites 2-5, 7-10, 11-14) are placed on the capsule cylindrical
%  surface (radial distance b+dis from the long axis) and point in
%  the corresponding outward Cartesian direction.
%
%  Attachment site map (X = long axis):
%     1      : Left pole       SPoint = (-(a+dis), 0, 0)         outward = -x
%     2-5    : Left  cap-cylinder junction (x = -(a-b))          (y+, z-, y-, z+)
%     6      : Right pole      SPoint = (+(a+dis), 0, 0)         outward = +x
%     7-10   : Right cap-cylinder junction (x = +(a-b))          (y+, z-, y-, z+)
%     11-14  : Cylinder midpoint           (x = 0)               (y+, z-, y-, z+)
%
%  Note: sites 2-5 and 7-10 are NOT on the curved hemisphere surface
%  proper — they are at the cylinder/cap junction circle, just inside
%  each cap.  They share the same outward direction as their cylinder
%  midpoint counterparts (11-14); only the X coordinate of the root
%  differs.  This naming is preserved for backward compatibility with
%  earlier Pili_Matrix configs.

switch NSite
    case 1   % Left pole
        rxP(1:Ntail) = -major_axis - dis - hTail.rxtail;
        ryP(1:Ntail) = hTail.rytail;
        rzP(1:Ntail) = hTail.rztail;
        SPoint = [-major_axis - dis, 0, 0];
        PV     = [-1, 0, 0];

    case 2   % Left hemisphere, y+
        rxP(1:Ntail) = hTail.rytail + minor_axis - major_axis;
        ryP(1:Ntail) = dis + hTail.rxtail + minor_axis;
        rzP(1:Ntail) = hTail.rztail;
        SPoint = [minor_axis - major_axis, minor_axis + dis, 0];
        PV     = [0, 1, 0];

    case 3   % Left hemisphere, z-
        rxP(1:Ntail) = hTail.rztail + minor_axis - major_axis;
        ryP(1:Ntail) = hTail.rytail;
        rzP(1:Ntail) = -hTail.rxtail - dis - minor_axis;
        SPoint = [minor_axis - major_axis, 0, -minor_axis - dis];
        PV     = [0, 0, -1];

    case 4   % Left hemisphere, y-
        rxP(1:Ntail) = hTail.rytail + minor_axis - major_axis;
        ryP(1:Ntail) = -dis - hTail.rxtail - minor_axis;
        rzP(1:Ntail) = hTail.rztail;
        SPoint = [minor_axis - major_axis, -minor_axis - dis, 0];
        PV     = [0, -1, 0];

    case 5   % Left hemisphere, z+
        rxP(1:Ntail) = hTail.rztail + minor_axis - major_axis;
        ryP(1:Ntail) = hTail.rytail;
        rzP(1:Ntail) = hTail.rxtail + dis + minor_axis;
        SPoint = [minor_axis - major_axis, 0, minor_axis + dis];
        PV     = [0, 0, 1];

    case 6   % Right pole
        rxP(1:Ntail) = major_axis + dis + hTail.rxtail;
        ryP(1:Ntail) = hTail.rytail;
        rzP(1:Ntail) = hTail.rztail;
        SPoint = [major_axis + dis, 0, 0];
        PV     = [1, 0, 0];

    case 7   % Right hemisphere, y+
        rxP(1:Ntail) = hTail.rytail + major_axis - minor_axis;
        ryP(1:Ntail) = dis + hTail.rxtail + minor_axis;
        rzP(1:Ntail) = hTail.rztail;
        SPoint = [major_axis - minor_axis, minor_axis + dis, 0];
        PV     = [0, 1, 0];

    case 8   % Right hemisphere, z-
        rxP(1:Ntail) = hTail.rztail + major_axis - minor_axis;
        ryP(1:Ntail) = hTail.rytail;
        rzP(1:Ntail) = -hTail.rxtail - dis - minor_axis;
        SPoint = [major_axis - minor_axis, 0, -minor_axis - dis];
        PV     = [0, 0, -1];

    case 9   % Right hemisphere, y-
        rxP(1:Ntail) = hTail.rytail + major_axis - minor_axis;
        ryP(1:Ntail) = -dis - hTail.rxtail - minor_axis;
        rzP(1:Ntail) = hTail.rztail;
        SPoint = [major_axis - minor_axis, -minor_axis - dis, 0];
        PV     = [0, -1, 0];

    case 10  % Right hemisphere, z+
        rxP(1:Ntail) = hTail.rztail + major_axis - minor_axis;
        ryP(1:Ntail) = hTail.rytail;
        rzP(1:Ntail) = hTail.rxtail + dis + minor_axis;
        SPoint = [major_axis - minor_axis, 0, minor_axis + dis];
        PV     = [0, 0, 1];

    case 11  % Middle cylinder, y+
        rxP(1:Ntail) = hTail.rytail;
        ryP(1:Ntail) = dis + hTail.rxtail + minor_axis;
        rzP(1:Ntail) = hTail.rztail;
        SPoint = [0, minor_axis + dis, 0];
        PV     = [0, 1, 0];

    case 12  % Middle cylinder, z-
        rxP(1:Ntail) = hTail.rztail;
        ryP(1:Ntail) = hTail.rytail;
        rzP(1:Ntail) = -hTail.rxtail - dis - minor_axis;
        SPoint = [0, 0, -minor_axis - dis];
        PV     = [0, 0, -1];

    case 13  % Middle cylinder, y-
        rxP(1:Ntail) = hTail.rytail;
        ryP(1:Ntail) = -dis - hTail.rxtail - minor_axis;
        rzP(1:Ntail) = hTail.rztail;
        SPoint = [0, -minor_axis - dis, 0];
        PV     = [0, -1, 0];

    case 14  % Middle cylinder, z+
        rxP(1:Ntail) = hTail.rztail;
        ryP(1:Ntail) = hTail.rytail;
        rzP(1:Ntail) = hTail.rxtail + dis + minor_axis;
        SPoint = [0, 0, minor_axis + dis];
        PV     = [0, 0, 1];

    otherwise
        error('pili_location_single:InvalidSite', ...
              'Attachment position %d is not available.', NSite);
end

end
