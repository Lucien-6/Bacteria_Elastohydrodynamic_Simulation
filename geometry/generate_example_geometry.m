%% generate_example_geometry.m
%  Regenerate the example custom geometry CSV files shipped with the
%  project (example_body_sphere.csv, example_body_capsule.csv,
%  example_pilus_linear.csv, example_pilus_helix.csv).
%
%  Description:
%    Produces a small, representative set of body point clouds and pili
%    centerlines so that config_custom_example.json can run out of the
%    box.  The generated geometries are intentionally compact (~300 body
%    points, ~40 pilus centerline points) to keep demo runtime short.
%
%    The pili centerlines are generated with a UNIFORM arc-length
%    discretisation, as required by custom_pili_building.m.
%
%  Usage:
%    >> cd geometry
%    >> generate_example_geometry
%
%  Outputs (written to the current folder):
%    example_body_sphere.csv        - 200-point sphere, R = 1 um
%    example_body_capsule.csv       - 300-point capsule (R=0.4um, L_cyl=1.7um)
%    example_pilus_linear.csv       - 41-point straight pilus, 4 um long
%    example_pilus_helix.csv        - 41-point helical pilus, 4 um long
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-05-04

close all; clear; clc;

out_dir = fileparts(mfilename('fullpath'));
if isempty(out_dir), out_dir = pwd; end
fprintf('Writing example geometry files to: %s\n', out_dir);

%% --- 1) Body: sphere (Fibonacci lattice, R = 1 um) ---
N_sphere = 200;
R_sphere = 1.0e-6;
idx      = (0:N_sphere-1)';
phi      = acos(1 - 2 * (idx + 0.5) / N_sphere);
theta    = pi * (1 + sqrt(5)) * (idx + 0.5);
sphere_pts = R_sphere * [sin(phi) .* cos(theta), ...
                         sin(phi) .* sin(theta), ...
                         cos(phi)];
write_csv(fullfile(out_dir, 'example_body_sphere.csv'), sphere_pts, ...
    'Example body: Fibonacci-lattice sphere, R = 1 um');

%% --- 2) Body: capsule (cylinder + hemispheres, long axis along body-X) ---
R_cap    = 0.4e-6;      % minor axis
L_cyl    = 1.7e-6;      % cylinder length (half-length 0.85 um)
N_cyl    = 180;         % lateral lattice points
N_cap    = 60;          % per hemisphere
nu       = sqrt(N_cyl); % aspect ratio for uniform-ish coverage
n_theta  = round(nu * sqrt(L_cyl / R_cap));
n_z      = max(round(N_cyl / max(n_theta, 1)), 4);
[theta_c, z_c] = ndgrid(linspace(0, 2*pi, n_theta+1), linspace(-L_cyl/2, L_cyl/2, n_z));
theta_c(end, :) = [];     % avoid duplicate column at 0 / 2*pi
cyl_x = z_c(:);
cyl_y = R_cap * cos(theta_c(:));
cyl_z = R_cap * sin(theta_c(:));
cyl_pts = [cyl_x, cyl_y, cyl_z];

idx2 = (0:N_cap-1)';
phi2 = acos(1 - (idx2 + 0.5) / N_cap);
th2  = pi * (1 + sqrt(5)) * (idx2 + 0.5);
hemi_unit = [cos(phi2), sin(phi2) .* cos(th2), sin(phi2) .* sin(th2)];   % points on +X hemisphere
hemi_plus  = R_cap * hemi_unit + [ L_cyl/2, 0, 0];
hemi_minus = R_cap * [-hemi_unit(:,1), hemi_unit(:,2), hemi_unit(:,3)] + [-L_cyl/2, 0, 0];

capsule_pts = [cyl_pts; hemi_plus; hemi_minus];
write_csv(fullfile(out_dir, 'example_body_capsule.csv'), capsule_pts, ...
    'Example body: capsule (R = 0.4 um, cylinder L = 1.7 um, long axis = X)');

%% --- 3) Pilus: straight, 4 um long, 41 UNIFORM samples ---
L_pilus = 4.0e-6;
N_seg   = 40;
s_lin   = linspace(0, L_pilus, N_seg + 1)';
lin_pts = [s_lin, zeros(N_seg+1, 1), zeros(N_seg+1, 1)];
write_csv(fullfile(out_dir, 'example_pilus_linear.csv'), lin_pts, ...
    'Example pilus: straight line along +X, L = 4 um, 41 uniform samples');

%% --- 4) Pilus: helix, 4 um long, 41 UNIFORM arc-length samples ---
R_h   = 0.15e-6;       % helix radius
P_h   = 1.0e-6;        % helix pitch (per turn, um)
% Arc length per turn: s = sqrt((2 pi R)^2 + P^2)
%   2 pi R = 0.9425e-6 m
%   P      = 1.0e-6 m
%   per-turn arc length ~ 1.3743e-6 m
% Parameterise by arc length s in [0, L_pilus]:
s_h    = linspace(0, L_pilus, N_seg + 1)';
per_turn_s = hypot(2*pi*R_h, P_h);
omega  = (2*pi) / per_turn_s;      % [rad / m]
pitch_slope = P_h / per_turn_s;    % [-]
helix_pts = [ pitch_slope * s_h, ...
              R_h * cos(omega * s_h), ...
              R_h * sin(omega * s_h) ];
write_csv(fullfile(out_dir, 'example_pilus_helix.csv'), helix_pts, ...
    'Example pilus: helix (R = 0.15 um, pitch = 1 um) along +X, L = 4 um, 41 uniform samples');

fprintf('Done. Example geometry files regenerated.\n');

% ---- local helper ----
function write_csv(path, pts, descr)
%WRITE_CSV  Write an N x 3 point cloud with a 2-line comment header.
    [fid, msg] = fopen(path, 'w', 'n', 'UTF-8');
    if fid < 0
        error('generate_example_geometry:WriteFail', ...
            'Failed to open %s: %s', path, msg);
    end
    fprintf(fid, '# %s\n', descr);
    fprintf(fid, '# Columns: X [m], Y [m], Z [m]\n');
    fprintf(fid, '%.9e,%.9e,%.9e\n', pts.');
    fclose(fid);
    fprintf('  Wrote %d points -> %s\n', size(pts, 1), path);
end
