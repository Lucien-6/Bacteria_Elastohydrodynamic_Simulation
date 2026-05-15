function [Capsule_str, NBody] = capsule_body_building(NBody, minor_axis, major_axis)
%CAPSULE_BODY_BUILDING  Construct capsule-shaped bacterial body surface mesh
%
%  Description:
%    Constructs a capsule-shaped bacterial body model consisting of two
%    hemispherical caps connected by a cylindrical middle section. Points
%    use Anton Semechko S2 Sampling Toolbox on the sphere caps (see
%    capsule_sphere_cap_points) and a regular grid for the cylinder.
%    Cylinder circumferential resolution TN is snapped to a multiple of 4
%    via TN = max(round(TN/4)*4, 4) after the area-based TN estimate.
%    Sphere points with x==0 (equator, strict) are duplicated so each
%    hemispherical cap receives a copy after translation.
%
%    In the Regularised Stokeslet Method (RSM), force points and evaluation
%    points are collocated, so only one set of surface coordinates is needed.
%
%  Syntax:
%    [Capsule_str, NBody] = capsule_body_building(NBody, minor_axis, major_axis)
%
%  Inputs:
%    NBody      - Desired number of surface points (approximate)
%    minor_axis - Semi-minor axis (radius of hemispherical caps) [m]
%    major_axis - Semi-major axis (half total body length) [m]
%
%  Outputs:
%    Capsule_str - Structure with surface point coordinates:
%                    .rxB, .ryB, .rzB  Surface point coordinates [NBody x 1] [m]
%    NBody       - Actual number of surface points (may differ from input)
%
%  Geometry:
%    Total length = 2 * major_axis, centered at origin along x-axis.
%    Capsule = two hemispheres (radius = minor_axis) + cylinder
%              (radius = minor_axis, length = 2*(major_axis - minor_axis)).
%
%  References:
%    [1] Cortez et al. (2005), Phys. Fluids 17, 031504.
%         (RSM: force and evaluation points collocated on surface)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2024-10-24
%  Updated: 2026-04-17 — S2 toolbox sphere caps; optional snap_sphere_cap_count
%  Updated: 2026-04-17 — duplicate points with x==0 (equator) onto both caps
%  Updated: 2026-04-20 — cylinder TN snapped to multiple of 4 (round(TN/4)*4)

%% Compute surface point allocation between sphere caps and cylinder
SA  = 4 * pi * minor_axis^2;                       % Total sphere surface area
CA  = 4 * pi * minor_axis * (major_axis - minor_axis);  % Cylinder lateral area
SPN_alloc = round(NBody * SA / (SA + CA));          % Area-based cap allocation
SPN = snap_sphere_cap_count(SPN_alloc);             % Snap to Ico/Quad native N if |diff|<=5
CPN = NBody - SPN;                                  % Remainder on cylinder

%% Generate sphere cap points (S2 toolbox; see capsule_sphere_cap_points)
V_unit = capsule_sphere_cap_points(SPN);
SLoc = minor_axis .* V_unit;

% Equator x==0 (strict): duplicate so left/right caps each receive one copy
xe = SLoc(:, 1);
idx_lt = xe < 0;
idx_gt = xe > 0;
idx_eq = xe == 0;
LocL = [SLoc(idx_lt, :); SLoc(idx_eq, :)];
LocR = [SLoc(idx_gt, :); SLoc(idx_eq, :)];
dx_cap = major_axis - minor_axis;
LocL(:, 1) = LocL(:, 1) - dx_cap;
LocR(:, 1) = LocR(:, 1) + dx_cap;

%% Generate cylinder surface points (regular grid)
if CPN > 0
    CDis  = sqrt(CA / CPN);
    Turns = ceil(2 * (major_axis - minor_axis) / CDis);
    CDis  = 2 * (major_axis - minor_axis) / Turns;
    TN    = ceil(CPN / Turns);
    TN    = max(round(TN / 4) * 4, 4);   % circumferential points per ring: multiple of 4
    Alpha = 0 : 2*pi/TN : (2*pi - 2*pi/TN);

    TempX = (minor_axis - major_axis + CDis/2) : CDis : (major_axis - minor_axis - CDis/2);
    TempY = minor_axis * sin(Alpha);
    TempZ = minor_axis * cos(Alpha);

    CLoc = zeros(TN * Turns, 3);
    CLoc(:, 1) = repelem(TempX, TN);
    CLoc(:, 2) = repmat(TempY', Turns, 1);
    CLoc(:, 3) = repmat(TempZ', Turns, 1);
else
    CDis  = NaN;
    CLoc  = NaN(3, 3);
    Turns = 0;
    TN    = 0;
end

%% Concatenate all surface points (left cap, cylinder, right cap)
rxB = [LocL(:, 1);  CLoc(:, 1);  LocR(:, 1)];
ryB = [LocL(:, 2);  CLoc(:, 2);  LocR(:, 2)];
rzB = [LocL(:, 3);  CLoc(:, 3);  LocR(:, 3)];

%% Remove NaN entries (from degenerate cylinder case)
rxB(isnan(rxB)) = [];
ryB(isnan(ryB)) = [];
rzB(isnan(rzB)) = [];

%% Final body point count (actual rows after snap, equator duplication, cylinder grid)
NBody = numel(rxB);

%% Pack output structure
Capsule_str = struct('rxB', rxB, 'ryB', ryB, 'rzB', rzB);

end
