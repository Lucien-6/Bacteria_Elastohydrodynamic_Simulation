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
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2024-10-10
%  Updated: 2026-02-11 — Removed legacy force-point fields for RSM

%% Extract common parameters
Pili_Type = PM(1);   % Morphology type
L_helix   = PM(2);   % Length of pilus [um]
theta     = PM(3);   % Deflection angle from Z-axis [rad]
gama      = PM(4);   % Deflection angle from X-axis [rad]

% Build deflection rotation matrix
RZ = [ cos(theta), -sin(theta), 0;
       sin(theta),  cos(theta), 0;
       0,           0,          1];
RX = [1, 0,          0;
      0, cos(gama),  sin(gama);
      0, -sin(gama), cos(gama)];
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

    case 4  % Expansive sinusoidal
        a = PM(5);
        b = PM(6);
        x = 0 : 0.001 : 4;
        y = a .* x .* cos(b * x);
        n = length(x);
        L = zeros(1, n);
        for i = 2:n
            L(i) = L(i-1) + sqrt((x(i)-x(i-1))^2 + (y(i)-y(i-1))^2);
        end
        [~, index] = min(abs(L_helix - L));
        XEnd   = x(index);
        rxtail = 0.0 : XEnd/(Ntail-1) : XEnd;
        rytail = a .* rxtail .* cos(b .* rxtail);
        rztail = zeros(1, Ntail);
        rxtail = rxtail .* 1e-6;
        rytail = rytail .* 1e-6;
        rztail = rztail .* 1e-6;

    case 5  % Conical spiral
        a = PM(5);
        b = PM(6);
        c = PM(7);
        t = 0 : 0.001 : 4;
        x = a .* t .* cos(c .* t);
        y = a .* t .* sin(c .* t);
        z = b .* t;
        n = length(t);
        L = zeros(1, n);
        for i = 2:n
            L(i) = L(i-1) + sqrt((x(i)-x(i-1))^2 + (y(i)-y(i-1))^2 + (z(i)-z(i-1))^2);
        end
        [~, index] = min(abs(L_helix - L));
        tEnd   = t(index);
        t      = 0.0 : tEnd/(Ntail-1) : tEnd;
        rytail = a .* t .* cos(c .* t);
        rztail = a .* t .* sin(c .* t);
        rxtail = b .* t;
        rxtail = rxtail .* 1e-6;
        rytail = rytail .* 1e-6;
        rztail = rztail .* 1e-6;

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
