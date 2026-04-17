function [Capsule_str, NBody] = capsule_body_building(NBody, minor_axis, major_axis)
%CAPSULE_BODY_BUILDING  Construct capsule/sphere-shaped bacterial body surface mesh
%
%  Description:
%    Constructs a capsule-shaped bacterial body model consisting of two
%    hemispherical caps connected by a cylindrical middle section.  When
%    major_axis == minor_axis the geometry degenerates to a sphere and a
%    dedicated branch is used (see below).
%
%    Sphere caps use Lebedev quadrature points (Laikov & Lebedev 1999),
%    which possess exact octahedral O_h symmetry (48 symmetry operations
%    including mirror symmetry about all three coordinate planes).  This
%    eliminates the systematic lateral drift caused by the asymmetric
%    Fibonacci spiral distribution previously used.
%
%    The cylinder uses a regular grid with circumferential resolution
%    forced to a multiple of 4, ensuring Y=0 and Z=0 mirror symmetry.
%
%    Capsule branch (major_axis > minor_axis):
%      Left and right caps are taken directly from the natural lx<0 and
%      lx>0 sub-sets of the Lebedev sphere (after discarding equatorial
%      lx=0 points, which would fall on the cap-cylinder junction).
%      Since Lebedev's generator `gen_oh` produces the two sub-sets by
%      explicit IEEE-754 sign flips, they are bit-exact mirror images of
%      each other, so left/right cap symmetry is inherited without any
%      additional mirroring step.
%
%    Sphere branch (major_axis == minor_axis):
%      All Lebedev points lie directly on the sphere; none is discarded.
%      The Lebedev rule is selected by ceiling against the *total* point
%      count (rather than the cap count used in the capsule branch), so
%      the chosen order accurately reflects the requested NBody.
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
%                 Must satisfy major_axis >= minor_axis.
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
%    When major_axis == minor_axis the cylinder degenerates to zero length
%    and the geometry is a sphere of radius minor_axis.
%
%  References:
%    [1] Cortez et al. (2005), Phys. Fluids 17, 031504.
%         (RSM: force and evaluation points collocated on surface)
%    [2] Lebedev V.I. & Laikov D.N. (1999), Doklady Math. 59(3), 477-481.
%         (Lebedev quadrature with O_h symmetry)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2024-10-24
%  Modified: 2026-04-16 — replaced Fibonacci spiral with Lebedev quadrature;
%                         use Lebedev's natural lx<0/lx>0 subsets directly
%                         (removed redundant explicit X-mirror);
%                         added dedicated sphere branch that keeps all
%                         Lebedev points when major_axis == minor_axis.

%% Input validation
if ~isscalar(minor_axis) || minor_axis <= 0
    error('capsule_body_building:InvalidMinorAxis', ...
        'minor_axis must be a positive scalar; got %.4e.', minor_axis);
end
if ~isscalar(major_axis) || major_axis < minor_axis
    error('capsule_body_building:InvalidMajorAxis', ...
        'major_axis (%.4e) must be >= minor_axis (%.4e).', ...
        major_axis, minor_axis);
end

%% Ensure Lebedev library is on path
persistent lebedev_path_added
if isempty(lebedev_path_added)
    leb_dir = fullfile(fileparts(mfilename('fullpath')), 'lebedev');
    if isfolder(leb_dir)
        addpath(leb_dir);
    end
    lebedev_path_added = true;
end

cap_offset = major_axis - minor_axis;
is_sphere  = (cap_offset == 0);

if is_sphere
    %% ============== Sphere branch (major == minor) ==============
    %  No cap-cylinder junction exists, so equatorial Lebedev points are
    %  valid sphere points and must NOT be discarded.  The Lebedev order
    %  is selected by ceiling against the total point count (mode='total').
    [SPN_lebedev, SPN] = select_lebedev_order(NBody, 'total');
    [lx, ly, lz, ~]    = ld_by_order(SPN_lebedev);

    rxB = lx(:) * minor_axis;
    ryB = ly(:) * minor_axis;
    rzB = lz(:) * minor_axis;

    NBody = SPN;
else
    %% ============== Capsule branch (major > minor) ==============
    %% Surface point allocation between sphere caps and cylinder
    SA = 4 * pi * minor_axis^2;
    CA = 4 * pi * minor_axis * cap_offset;
    SPN_target = round(NBody * SA / (SA + CA));

    %% Select Lebedev rule (ceiling strategy against cap point count)
    %  For each Lebedev order, the actual cap point count is 2*N_hemisphere
    %  where N_hemisphere = number of points with lx > 0 (equatorial lx=0
    %  points are discarded — that region is covered by the cylinder).
    %  Ceiling guarantees the cap is never undersampled relative to the
    %  area-proportional target, and all cylinder allocations stay >= 0.
    [SPN_lebedev, SPN] = select_lebedev_order(SPN_target, 'cap');
    CPN = max(NBody - SPN, 0);

    %% Generate sphere cap points (Lebedev quadrature, O_h symmetry)
    %  Lebedev points are bit-exactly symmetric about the X=0, Y=0, Z=0
    %  planes by construction (gen_oh uses IEEE-754 sign flips).  We take
    %  the natural lx<0 and lx>0 sub-sets directly — the two sets are
    %  exact mirror images of each other, so left/right cap symmetry is
    %  inherited without any additional mirroring step.  Equatorial
    %  (lx=0) points are discarded: they would fall on the cap-cylinder
    %  junction and create self-duplicates under any mirror assignment.
    [lx, ly, lz, ~] = ld_by_order(SPN_lebedev);

    mask_right = (lx > 0);
    mask_left  = (lx < 0);

    rx_right = lx(mask_right) * minor_axis + cap_offset;
    ry_right = ly(mask_right) * minor_axis;
    rz_right = lz(mask_right) * minor_axis;

    rx_left  = lx(mask_left)  * minor_axis - cap_offset;
    ry_left  = ly(mask_left)  * minor_axis;
    rz_left  = lz(mask_left)  * minor_axis;

    %% Generate cylinder surface points (regular grid, TN multiple of 4)
    if CPN > 0 && CA > 0
        CDis  = sqrt(CA / CPN);
        Turns = ceil(2 * cap_offset / CDis);
        Turns = max(Turns, 1);
        CDis  = 2 * cap_offset / Turns;
        TN    = ceil(CPN / Turns);
        TN    = max(ceil(TN / 4) * 4, 4);
        Alpha = (0 : TN-1) * (2*pi / TN);

        TempX = (-cap_offset + CDis/2) : CDis : (cap_offset - CDis/2);
        Turns = numel(TempX);
        TempY = minor_axis * sin(Alpha);
        TempZ = minor_axis * cos(Alpha);

        CLoc = zeros(TN * Turns, 3);
        CLoc(:, 1) = repelem(TempX, TN)';
        CLoc(:, 2) = repmat(TempY', Turns, 1);
        CLoc(:, 3) = repmat(TempZ', Turns, 1);
    else
        CLoc  = zeros(0, 3);
        Turns = 0;
        TN    = 0;
    end

    %% Concatenate all surface points (left cap, cylinder, right cap)
    rxB = [rx_left;  CLoc(:, 1);  rx_right];
    ryB = [ry_left;  CLoc(:, 2);  ry_right];
    rzB = [rz_left;  CLoc(:, 3);  rz_right];

    NBody = SPN + Turns * TN;
end

%% Pack output structure
Capsule_str = struct('rxB', rxB, 'ryB', ryB, 'rzB', rzB);

end


%% ==================== Local Helper Function ====================

function [leb_order, spn_actual] = select_lebedev_order(spn_target, mode)
%SELECT_LEBEDEV_ORDER  Choose the smallest Lebedev order whose point
%  count is at least equal to the target (ceiling strategy).  If the
%  target exceeds the largest available Lebedev rule, the maximum
%  available order is returned.
%
%  Inputs:
%    spn_target - Target point count (integer >= 0)
%    mode       - 'cap'   : match against cap count (2*sum(lx>0));
%                           used for the capsule branch where equatorial
%                           Lebedev points are discarded.
%                 'total' : match against total Lebedev point count;
%                           used for the sphere branch where all Lebedev
%                           points are kept.
%
%  Since Lebedev's O_h symmetry guarantees #{lx>0} == #{lx<0}, the cap
%  count is exactly 2 * sum(lx > 0).

persistent LEBEDEV_LUT

if isempty(LEBEDEV_LUT)
    orders = [6 14 26 38 50 74 86 110 146 170 ...
              194 230 266 302 350 434 590 770];
    n_cap  = zeros(size(orders));
    for k = 1:numel(orders)
        [lx, ~, ~, ~] = ld_by_order(orders(k));
        n_cap(k) = 2 * sum(lx > 0);
    end
    % Column 1 = total Lebedev point count (also = order label by construction)
    % Column 2 = cap point count after discarding equator
    LEBEDEV_LUT = [orders(:), n_cap(:)];
end

switch mode
    case 'total'
        col = 1;
    case 'cap'
        col = 2;
    otherwise
        error('select_lebedev_order:InvalidMode', ...
            'mode must be ''cap'' or ''total''; got ''%s''.', mode);
end

idx = find(LEBEDEV_LUT(:, col) >= spn_target, 1, 'first');
if isempty(idx)
    idx = size(LEBEDEV_LUT, 1);
end
leb_order  = LEBEDEV_LUT(idx, 1);
spn_actual = LEBEDEV_LUT(idx, col);

end
