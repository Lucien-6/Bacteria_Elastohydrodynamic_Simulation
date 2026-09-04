function Pili_str = pili_type_building(PM, Ntail)
%PILI_TYPE_BUILDING  Generate single-pilus 3D morphology coordinates
%
%  Description:
%    Generates the 3D coordinates for a single pilus with the specified
%    morphology type.  The pilus is created in a local coordinate system
%    starting from the origin, then rotated according to deflection angles.
%
%  Syntax:
%    Pili_str = pili_type_building(PM, Ntail)
%
%  Inputs:
%    PM    - [7 x 1] Parameter vector:
%              PM(1) - Morphology type (1=linear, 2=parabolic, 3=circular,
%                       4=sinusoidal, 5=conical spiral)
%              PM(2) - Length [um]
%              PM(3) - Deflection angle relative to Z-axis [rad]
%              PM(4) - Deflection angle relative to X-axis [rad]
%              PM(5) - Shape parameter 1 (type-dependent)
%              PM(6) - Shape parameter 2 (type-dependent)
%              PM(7) - Shape parameter 3 (type-dependent)
%    Ntail - Number of discretization points
%
%  Outputs:
%    Pili_str - Structure with pilus coordinates:
%                 .rxtail, .rytail, .rztail  Surface points [1 x Ntail] [m]
%
%  Shape parameters by type:
%    Type 1 (Linear)   : No additional parameters needed.
%    Type 2 (Parabolic): PM(5) = curvature scaling factor K.
%    Type 3 (Circular) : PM(5) = arc angle [rad].
%    Type 4 (Sinusoid) : PM(5) = amplitude factor a, PM(6) = frequency b.
%    Type 5 (Spiral)   : PM(5) = a, PM(6) = b, PM(7) = c (spiral params).
%
%  Sampling strategy:
%    Types 1, 3 naturally produce arc-length-uniform nodes.
%    Type 2 uses an explicit arc-length → x polynomial inversion.
%    Types 4, 5 use arc-length uniform resampling via PCHIP interpolation
%      to ensure all nodes are equally spaced in arc length, which is
%      required by the inextensible Kirchhoff rod framework (constant ds).
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2024-10-10
%  Updated: 2026-09-04 (v8.3.0: RX made right-handed)
%  Updated: 2026-09-04 (v8.5.0: type-2 parabola anchored at the origin)

%% Extract common parameters
Pili_Type = PM(1);   % Morphology type
L_helix   = PM(2);   % Length of pilus [um]
theta     = PM(3);   % Deflection angle from Z-axis [rad]
gama      = PM(4);   % Deflection angle from X-axis [rad]

% Build deflection rotation matrix.
% Both RZ and RX are standard RIGHT-HANDED rotations (positive angle =
% counter-clockwise when looking down the axis toward the origin), so
% RM = RX * RZ is a proper rotation with det = +1 and the same sign
% convention as pili_directors / quat_to_rotmat.  (Before v8.3.0 RX used
% a left-handed sign for gama; configurations with a non-zero
% deflection_x therefore see the sense of that twist reversed.)
%
% Parametric-mode semantics (for cross-reference with custom_pili_building.m):
%   - The pilus is constructed along the +X axis, so X IS the initial
%     tangent direction here.
%   - deflection_z (theta) rotates about Z (perpendicular to the tangent)
%     => TILT of the tangent direction.
%   - deflection_x (gama) rotates about X (the tangent itself)
%     => TWIST around the pilus axis (a linear pilus is invariant).
% custom_pili_building.m applies rotations in a local basis with e3 equal
% to the aligned tangent, so its deflection_z/_x axes are SWAPPED relative
% to the parametric names.  When porting a non-zero deflection between
% modes, swap deflection_x <-> deflection_z.
RZ = [ cos(theta), -sin(theta), 0;
       sin(theta),  cos(theta), 0;
       0,           0,          1];
RX = [1, 0,          0;
      0, cos(gama), -sin(gama);
      0, sin(gama),  cos(gama)];
RM = RX * RZ;

%% Generate morphology-specific coordinates
switch Pili_Type

    case 1  % Linear (straight line)
        L_helix = L_helix * 1e-6;
        rxtail  = 0.0 : L_helix/(Ntail-1) : L_helix;
        rytail  = zeros(1, Ntail);
        rztail  = zeros(1, Ntail);

    case 2  % Parabolic
        % K is a curvature scaling factor that controls how strongly the
        % pilus curves.  The arc-length parameter is scaled by K so that
        % different curvatures can be generated from a single template curve.
        %
        % The 6th-order polynomial below maps the scaled arc-length parameter
        % (range [0, K*L_um]) to the projected x-coordinate of the parabola.
        % Coefficients were obtained by least-squares fitting of the analytic
        % arc-length integral for y = x^2 (parabola) using MATLAB's polyfit
        % over the range [0, 4] (unitless).  Valid input range: K*L_um <= 4.
        % The polynomial is then inverted: given arc length s, find x such
        % that arc_length(0->x) = s.
        K       = PM(5);
        L_helix = 0 : K*L_helix/(Ntail-1) : K*L_helix;
        XEnd    = -0.0002146*L_helix.^6 + 0.004265*L_helix.^5 ...
                  - 0.0352*L_helix.^4 + 0.1586*L_helix.^3 ...
                  - 0.4459*L_helix.^2 + 1.087*L_helix - 0.004127;
        % The least-squares fit carries a residual constant term (-0.004127)
        % so that x(s=0) != 0 and the template's first node would sit
        % ~4.1 nm/K off the origin.  Every downstream consumer
        % (pili_building, custom_pili_building, enforce_clamped_geometry)
        % assumes the template starts exactly at the origin, so anchor the
        % inverted curve there (v8.5.0).  The shift is uniform in x and does
        % not change the node spacing.
        XEnd    = XEnd - XEnd(1);
        rxtail  = XEnd / K;
        rytail  = rxtail.^2.0;
        rxtail  = rxtail .* 1e-6;
        rytail  = rytail .* 1e-6;
        rztail  = zeros(1, Ntail);

    case 3  % Circular arc
        L_helix      = L_helix * 1e-6;
        Angle_Circle = PM(5);
        Dtheta = 0 : Angle_Circle/(Ntail-1) : Angle_Circle;
        rxtail = L_helix / Angle_Circle .* sin(Dtheta);
        rytail = L_helix / Angle_Circle .* (1.0 - cos(Dtheta));
        rztail = zeros(1, Ntail);

    case 4  % Expansive sinusoidal — arc-length uniform resampling
        a = PM(5);
        b = PM(6);
        L_helix_m = L_helix * 1e-6;         % Target arc length [m]

        % Fine parameter sweep with sub-nanometer arc-length resolution
        x_fine = 0 : 1e-4 : 6;
        y_fine = a .* x_fine .* cos(b .* x_fine);
        n_fine = length(x_fine);

        % Cumulative arc length [m]
        L_fine_m = zeros(1, n_fine);
        for i = 2:n_fine
            L_fine_m(i) = L_fine_m(i-1) + 1e-6 * sqrt((x_fine(i)-x_fine(i-1))^2 + ...
                (y_fine(i)-y_fine(i-1))^2);
        end

        % Truncate to the region where arc length >= target (ensures extrapolation-free)
        idx_past = find(L_fine_m >= L_helix_m, 1, 'first');
        if isempty(idx_past)
            error('pili_type_building:SpiralTooShort', ...
                'Type 4: arc length at max parameter (6) is %.4f um, shorter than target %.4f um.', ...
                L_fine_m(end)/1e-6, L_helix);
        end
        idx_use = idx_past;                  % Include the point just past target

        % Uniform arc-length targets [m] — use exact target L_helix_m
        s_uniform = (0 : Ntail-1) * L_helix_m / (Ntail-1);

        % Directly interpolate coordinates from arc length (PCHIP)
        rxtail = interp1(L_fine_m(1:idx_use), x_fine(1:idx_use) * 1e-6, ...
            s_uniform, 'pchip');
        rytail = interp1(L_fine_m(1:idx_use), y_fine(1:idx_use) * 1e-6, ...
            s_uniform, 'pchip');
        rztail = zeros(1, Ntail);

    case 5  % Conical spiral — arc-length uniform resampling
        a = PM(5);
        b = PM(6);
        c = PM(7);
        L_helix_m = L_helix * 1e-6;         % Target arc length [m]

        % Fine parameter sweep with sub-nanometer arc-length resolution
        t_fine = 0 : 1e-4 : 6;
        x_fine = b .* t_fine;                % axial [μm]
        y_fine = a .* t_fine .* cos(c .* t_fine);  % transverse [μm]
        z_fine = a .* t_fine .* sin(c .* t_fine);  % transverse [μm]

        % Cumulative arc length [m] — compute directly in meters
        n_fine = length(t_fine);
        L_fine_m = zeros(1, n_fine);
        for i = 2:n_fine
            L_fine_m(i) = L_fine_m(i-1) + 1e-6 * sqrt((x_fine(i)-x_fine(i-1))^2 + ...
                (y_fine(i)-y_fine(i-1))^2 + (z_fine(i)-z_fine(i-1))^2);
        end

        % Truncate to the region where arc length >= target (ensures extrapolation-free)
        idx_past = find(L_fine_m >= L_helix_m, 1, 'first');
        if isempty(idx_past)
            error('pili_type_building:SpiralTooShort', ...
                'Type 5: arc length at max parameter (6) is %.4f um, shorter than target %.4f um.', ...
                L_fine_m(end)/1e-6, L_helix);
        end
        idx_use = idx_past;                  % Include the point just past target

        % Uniform arc-length targets [m] — use exact target L_helix_m
        s_uniform = (0 : Ntail-1) * L_helix_m / (Ntail-1);

        % Directly interpolate coordinates from arc length (PCHIP)
        rxtail = interp1(L_fine_m(1:idx_use), x_fine(1:idx_use) * 1e-6, ...
            s_uniform, 'pchip');
        rytail = interp1(L_fine_m(1:idx_use), y_fine(1:idx_use) * 1e-6, ...
            s_uniform, 'pchip');
        rztail = interp1(L_fine_m(1:idx_use), z_fine(1:idx_use) * 1e-6, ...
            s_uniform, 'pchip');

    otherwise
        error('pili_type_building:InvalidType', ...
              'Pili morphology type %d is not available.', Pili_Type);
end

%% Apply deflection rotation
T      = RM * [rxtail; rytail; rztail];
rxtail = T(1, :);
rytail = T(2, :);
rztail = T(3, :);

%% Pack output structure
Pili_str = struct('rxtail', rxtail, 'rytail', rytail, 'rztail', rztail);

end
