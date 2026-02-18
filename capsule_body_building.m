function [Capsule_str, NBody] = capsule_body_building(NBody, minor_axis, major_axis)
%CAPSULE_BODY_BUILDING  Construct capsule-shaped bacterial body surface mesh
%
%  Description:
%    Constructs a capsule-shaped bacterial body model consisting of two
%    hemispherical caps connected by a cylindrical middle section. Points
%    are distributed uniformly on the surface using spiral point placement
%    for spheres and a regular grid for the cylinder.
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
%  Updated: 2026-02-11 — Removed legacy force-point (shift) code for RSM

%% Compute surface point allocation between sphere caps and cylinder
SA  = 4 * pi * minor_axis^2;                       % Total sphere surface area
CA  = 4 * pi * minor_axis * (major_axis - minor_axis);  % Cylinder lateral area
SPN = round(NBody * SA / (SA + CA));                % Points on sphere caps
CPN = NBody - SPN;                                  % Points on cylinder

%% Generate sphere cap points (spiral distribution)
SLoc = zeros(SPN, 3);
for ii = 1:SPN
    cos_phi = -1.0 + (2.0 * ii - 1.0) / SPN;
    cos_phi = max(-1, min(1, cos_phi));
    phi   = acos(cos_phi);
    theta = sqrt(SPN * pi) * phi;
    SLoc(ii, 3) = minor_axis * cos(theta) * sin(phi);
    SLoc(ii, 2) = minor_axis * sin(theta) * sin(phi);
    SLoc(ii, 1) = minor_axis * cos(phi);
end

%% Generate cylinder surface points (regular grid)
if CPN > 0
    CDis  = sqrt(CA / CPN);
    Turns = ceil(2 * (major_axis - minor_axis) / CDis);
    CDis  = 2 * (major_axis - minor_axis) / Turns;
    TN    = ceil(CPN / Turns);
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

%% Assemble hemisphere caps onto cylinder ends
SLoc(SLoc(:, 1) <= 0, 1) = SLoc(SLoc(:, 1) <= 0, 1) - major_axis + minor_axis;
SLoc(SLoc(:, 1) >  0, 1) = SLoc(SLoc(:, 1) >  0, 1) + major_axis - minor_axis;

%% Concatenate all surface points (left cap, cylinder, right cap)
rxB = [SLoc(SLoc(:, 1) <= 0, 1);  CLoc(:, 1);  SLoc(SLoc(:, 1) > 0, 1)];
ryB = [SLoc(SLoc(:, 1) <= 0, 2);  CLoc(:, 2);  SLoc(SLoc(:, 1) > 0, 2)];
rzB = [SLoc(SLoc(:, 1) <= 0, 3);  CLoc(:, 3);  SLoc(SLoc(:, 1) > 0, 3)];

%% Remove NaN entries (from degenerate cylinder case)
rxB(isnan(rxB)) = [];
ryB(isnan(ryB)) = [];
rzB(isnan(rzB)) = [];

%% Final body point count
NBody = SPN + Turns * TN;

%% Pack output structure
Capsule_str = struct('rxB', rxB, 'ryB', ryB, 'rzB', rzB);

end
