%% generate_ecoli_geometry.m
%  Generate a realistic E. coli geometry for the single-polar-flagellum
%  swimming case (configs/config_ecoli_single_flagellum.json).
%
%  Description:
%    Writes two CSV files beside this script:
%      ecoli_body_capsule.csv      - cell body point cloud (sphero-cylinder,
%                                    ~300 points, laid out by the project's
%                                    capsule_body_building scheme)
%      ecoli_flagellum_helix.csv   - flagellar filament centreline
%
%    Cell body (Berg, "E. coli in Motion", 2004; typical swimming cell):
%      total length   L_cell = 2.5 um   (pole to pole)
%      diameter       D_cell = 0.88 um
%    The long axis is body-X.  The tail pole is at x = -L_cell/2.
%
%    Flagellar filament (left-handed "normal" polymorph; Turner, Ryu & Berg,
%    J. Bacteriol. 182, 2793, 2000):
%      helical wavelength (pitch)   2.22 um
%      helix diameter               0.41 um
%      contour length               n_turns turns (default 3.0 -> ~7.5 um;
%                                   filament lengths of 5-10 um are typical)
%    The helix amplitude is ramped smoothly from 0 to its full value over the
%    first RAMP_TURNS turns (sin^2 envelope).  The centreline therefore leaves
%    the cell pole exactly along the cell axis with zero curvature, which is
%    how the flexible hook + filament root is represented when the model
%    imposes a bending-locked root (see the notes in the JSON config).
%
%    The filament is sampled with a UNIFORM arc-length discretisation, as
%    required by custom_pili_building.m.  The first point is the origin (root)
%    and the root tangent is +X; the solver rotates it onto the attachment
%    site direction (-X at the tail pole), which preserves handedness.
%
%  Usage:
%    >> cd geometry
%    >> generate_ecoli_geometry
%
%  Outputs:
%    ecoli_body_capsule.csv
%    ecoli_flagellum_helix.csv
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-19

close all; clear; clc;

out_dir = fileparts(mfilename('fullpath'));
if isempty(out_dir), out_dir = pwd; end
fprintf('Writing E. coli geometry files to: %s\n', out_dir);

%% --- User-adjustable parameters ---------------------------------------
L_cell      = 2.5e-6;     % [m] cell length, pole to pole
D_cell      = 0.88e-6;    % [m] cell diameter
lambda_h    = 2.22e-6;    % [m] helical pitch (wavelength along axis)
D_helix     = 0.41e-6;    % [m] helix diameter (centreline)
n_turns     = 3.0;        % [-] number of helical turns
ramp_turns  = 0.5;        % [-] turns over which the helix amplitude ramps up
seg_len     = 0.10e-6;    % [m] source chord length (the solver re-samples to seg_per_um)
handedness  = -1;         % [-] -1 = left-handed (normal E. coli filament)
N_body      = 300;        % [-] requested body surface points (actual count may differ slightly)

%% --- 1) Body: sphero-cylinder, long axis along body-X ------------------
%  Points are laid out with the project's native capsule scheme
%  (capsule_body_building: S2-toolbox sphere sample for the two caps +
%  regular ring grid on the cylinder), identical to the parametric mode.
if L_cell <= D_cell
    error('generate_ecoli_geometry:BadBody', 'L_cell must exceed D_cell.');
end
proj_root = fileparts(out_dir);
addpath(genpath(fullfile(proj_root, 'src')));
[Cap, N_body_actual] = capsule_body_building(N_body, D_cell/2, L_cell/2);
body_pts = [Cap.rxB(:), Cap.ryB(:), Cap.rzB(:)];
fprintf('  Body: requested %d points, generated %d points\n', N_body, N_body_actual);
write_csv(fullfile(out_dir, 'ecoli_body_capsule.csv'), body_pts, ...
    sprintf(['E. coli body: capsule via capsule_body_building (%d points), ', ...
    'L = %.2f um, D = %.2f um, long axis = X'], ...
    N_body_actual, L_cell*1e6, D_cell*1e6));

%% --- 2) Flagellar filament: left-handed helix with root ramp -----------
R_h      = D_helix / 2;
t_end    = 2*pi * n_turns;
t_ramp   = 2*pi * ramp_turns;
t_fine   = linspace(0, t_end, 20001)';
env      = ones(size(t_fine));
in_ramp  = t_fine < t_ramp;
env(in_ramp) = sin(0.5*pi * t_fine(in_ramp) / t_ramp).^2;   % C1 smooth, env'(0)=0

x_f = lambda_h * t_fine / (2*pi);
y_f = R_h * env .* cos(t_fine);
z_f = handedness * R_h * env .* sin(t_fine);
pts_fine = [x_f, y_f, z_f];

s_fine = [0; cumsum(vecnorm(diff(pts_fine, 1, 1), 2, 2))];
L_total = s_fine(end);
N_seg   = max(round(L_total / seg_len), 4);
s_new   = linspace(0, L_total, N_seg + 1)';
helix_pts = [interp1(s_fine, x_f, s_new, 'pchip'), ...
             interp1(s_fine, y_f, s_new, 'pchip'), ...
             interp1(s_fine, z_f, s_new, 'pchip')];

fprintf('  Filament: contour length %.3f um, axial length %.3f um, %d segments\n', ...
    L_total*1e6, x_f(end)*1e6, N_seg);
write_csv(fullfile(out_dir, 'ecoli_flagellum_helix.csv'), helix_pts, ...
    sprintf(['E. coli flagellar filament: left-handed helix, pitch %.2f um, ', ...
    'diameter %.2f um, %.1f turns, %.2f-turn root ramp, contour %.2f um'], ...
    lambda_h*1e6, D_helix*1e6, n_turns, ramp_turns, L_total*1e6));

fprintf('Done.  Attachment site (tail pole): x = %.4e m, direction [-1 0 0].\n', ...
    -L_cell/2);

% ---- local helper ----
function write_csv(path, pts, descr)
%WRITE_CSV  Write an N x 3 point cloud with a 2-line comment header.
    [fid, msg] = fopen(path, 'w', 'n', 'UTF-8');
    if fid < 0
        error('generate_ecoli_geometry:WriteFail', ...
            'Failed to open %s: %s', path, msg);
    end
    fprintf(fid, '# %s\n', descr);
    fprintf(fid, '# Columns: X [m], Y [m], Z [m]\n');
    fprintf(fid, '%.9e,%.9e,%.9e\n', pts.');
    fclose(fid);
    fprintf('  Wrote %d points -> %s\n', size(pts, 1), path);
end
