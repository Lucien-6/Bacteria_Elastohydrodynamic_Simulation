function Pili = pili_building(Pili_Matrix, seg_per_um, major_axis, minor_axis, dis, E, r_pili)
%PILI_BUILDING  Initialize pili structures for discrete elastic rod simulation
%
%  Description:
%    Constructs pili structures based on the pili geometry generator,
%    adding bending stiffness and empty curvature-binormal placeholders.
%    init_coupled_params recomputes the operational curvature data from the
%    initial pilus geometry.
%
%  Syntax:
%    Pili = pili_building(Pili_Matrix, seg_per_um, major_axis, minor_axis, dis, E, r_pili)
%
%  Inputs:
%    Pili_Matrix  - [8 x M] Pili parameter matrix (M = number of pili)
%        Row 1 : Morphology type (1=linear, 2=parabolic, 3=circular,
%                 4=sinusoidal, 5=conical spiral)
%        Row 2 : Length [um]
%        Row 3 : Deflection angle relative to Z-axis [rad]
%        Row 4 : Deflection angle relative to X-axis [rad]
%        Row 5-7: Shape control parameters (type-dependent)
%        Row 8 : Attachment position on bacterial body (1-14)
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
%        .kb         Placeholder [3 x 0]; curvature recomputed by init_coupled_params
%        .kb0        Placeholder [3 x 0]; curvature recomputed by init_coupled_params
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
%  See also: pili_type_building, init_coupled_params
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-01-23
%  Updated: 2026-09-04 (v8.3.0: attachment-site maps are proper rotations;
%                       reflections at sites 1,2,5,7,10,11,14 removed)

%% Input validation
if nargin < 7
    error('pili_building:InsufficientInputs', ...
          'pili_building requires at least 7 input arguments.');
end
if size(Pili_Matrix, 1) ~= 8
    error('pili_building:InvalidPiliMatrixRows', ...
        ['Pili_Matrix must have exactly 8 rows; got %d. ', ...
         'Pilus stiffness is derived from the global E and r_pili values.'], ...
        size(Pili_Matrix, 1));
end

%% Calculate default bending stiffness (EI = E * pi * r^4 / 4)
I_pili    = pi * r_pili^4 / 4;   % Second moment of area [m^4]
A_default = E * I_pili;           % Bending stiffness [N*m^2]

%% Get number of pili
NPili = size(Pili_Matrix, 2);

%% Pre-allocate output structure array
template = struct('N', [], 'L', [], 'ell', [], ...
    'x', [], 'v', [], 'kb', [], 'kb0', [], ...
    'A', [], 'r', [], 'f_bend', [], ...
    'attach_idx', [], 'root_local', [], 'root_dir', [], ...
    'type', [], 'is_clamped', []);
Pili = repmat(template, 1, NPili);

%% Construct each pilus
for j = 1:NPili

    % --- Extract parameters ---
    pili_type  = Pili_Matrix(1, j);
    L_um       = Pili_Matrix(2, j);          % Length [um]
    L          = L_um * 1e-6;                % Length [m]
    attach_idx = Pili_Matrix(8, j);

    % Bending stiffness is uniform global material metadata.
    A = A_default;

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

    % --- Curvature binormal placeholders ---
    % All curvature information is recomputed from scratch by
    % init_coupled_params (params.kb0{k}, [N+1 x 3], Euler-angle based);
    % these struct fields are kept only as placeholders for shape
    % compatibility, mirroring custom_pili_building.m.
    kb  = zeros(3, 0);
    kb0 = zeros(3, 0);

    % --- Initialize velocities and forces to zero ---
    v      = zeros(3, N);
    f_bend = zeros(3, N);

    % --- Store in structure ---
    Pili(j).N          = N;
    Pili(j).L          = L;
    Pili(j).ell        = ell;
    Pili(j).x          = x;
    Pili(j).v          = v;
    Pili(j).kb         = kb;            % placeholder; recomputed by init_coupled_params
    Pili(j).kb0        = kb0;           % placeholder; recomputed by init_coupled_params
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
%
%  Transformation (v8.3.0):
%    The template pilus from pili_type_building is built along +x
%    (tangent = +x, in-plane bending direction = +y, out-of-plane = +z).
%    Each site applies a PROPER rotation R (det R = +1) followed by a
%    translation to the root point:
%        p_body = R * p_template + SPoint,   R = [c1, c2, c1 x c2]
%    where c1 = PV (image of the template tangent +x) and c2 is the
%    image of the template bending direction +y.  The third column is
%    forced to c1 x c2 so that R is always a rotation.  The pre-v8.3.0
%    table hard-coded coordinate swaps; for sites 1, 2, 5, 7, 10, 11 and
%    14 those maps had det = -1 (reflections) and silently inverted the
%    handedness of chiral templates (conical spiral, type 5) and of the
%    deflection_x twist.  The rotations below agree with the legacy maps
%    on the template's bending plane z = 0, so every planar template
%    (types 1-4 with deflection_x = 0) is placed identically to before;
%    only the out-of-plane component changes sign at the seven affected
%    sites, which is exactly the chirality correction.

ex = [1, 0, 0];
ey = [0, 1, 0];
ez = [0, 0, 1];

switch NSite
    % ---- Poles ------------------------------------------------------
    case 1   % Left pole
        c1 = -ex;   c2 =  ey;
        SPoint = -(major_axis + dis) * ex;
    case 6   % Right pole
        c1 =  ex;   c2 =  ey;
        SPoint =  (major_axis + dis) * ex;

    % ---- Left cap/cylinder junction (x = -(a-b)) --------------------
    case 2   % Left hemisphere, y+
        c1 =  ey;   c2 =  ex;
        SPoint = -(major_axis - minor_axis) * ex + (minor_axis + dis) * c1;
    case 3   % Left hemisphere, z-
        c1 = -ez;   c2 =  ey;
        SPoint = -(major_axis - minor_axis) * ex + (minor_axis + dis) * c1;
    case 4   % Left hemisphere, y-
        c1 = -ey;   c2 =  ex;
        SPoint = -(major_axis - minor_axis) * ex + (minor_axis + dis) * c1;
    case 5   % Left hemisphere, z+
        c1 =  ez;   c2 =  ey;
        SPoint = -(major_axis - minor_axis) * ex + (minor_axis + dis) * c1;

    % ---- Right cap/cylinder junction (x = +(a-b)) -------------------
    case 7   % Right hemisphere, y+
        c1 =  ey;   c2 =  ex;
        SPoint =  (major_axis - minor_axis) * ex + (minor_axis + dis) * c1;
    case 8   % Right hemisphere, z-
        c1 = -ez;   c2 =  ey;
        SPoint =  (major_axis - minor_axis) * ex + (minor_axis + dis) * c1;
    case 9   % Right hemisphere, y-
        c1 = -ey;   c2 =  ex;
        SPoint =  (major_axis - minor_axis) * ex + (minor_axis + dis) * c1;
    case 10  % Right hemisphere, z+
        c1 =  ez;   c2 =  ey;
        SPoint =  (major_axis - minor_axis) * ex + (minor_axis + dis) * c1;

    % ---- Cylinder midpoint (x = 0) ----------------------------------
    case 11  % Middle cylinder, y+
        c1 =  ey;   c2 =  ex;
        SPoint = (minor_axis + dis) * c1;
    case 12  % Middle cylinder, z-
        c1 = -ez;   c2 =  ey;
        SPoint = (minor_axis + dis) * c1;
    case 13  % Middle cylinder, y-
        c1 = -ey;   c2 =  ex;
        SPoint = (minor_axis + dis) * c1;
    case 14  % Middle cylinder, z+
        c1 =  ez;   c2 =  ey;
        SPoint = (minor_axis + dis) * c1;

    otherwise
        error('pili_location_single:InvalidSite', ...
              'Attachment position %d is not available.', NSite);
end

% Assemble the proper rotation and apply it to the template
R  = [c1(:), c2(:), cross(c1(:), c2(:))];      % columns = images of ex, ey, ez
if abs(det(R) - 1) > 1e-12
    error('pili_location_single:NotARotation', ...
          'Internal error: site %d transform has det = %.3g (expected +1).', ...
          NSite, det(R));
end

P_template = [hTail.rxtail(:).'; hTail.rytail(:).'; hTail.rztail(:).'];  % [3 x Ntail]
P_body     = R * P_template + SPoint(:);                                 % [3 x Ntail]

rxP(1:Ntail) = P_body(1, :);
ryP(1:Ntail) = P_body(2, :);
rzP(1:Ntail) = P_body(3, :);
PV = c1;

end
