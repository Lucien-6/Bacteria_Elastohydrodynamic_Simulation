%% generate_chlamydomonas_geometry.m
%  Generate the geometry of a biflagellate Chlamydomonas reinhardtii
%  (wild type CC-124 / CC-125) for the active-moment swimming cases
%  (configs/config_chlamy_*.json).
%
%  Description:
%    Writes the following CSV files beside this script (SHAPE_MODE selects
%    the flagellar rest shape and the file-name suffix):
%      chlamy_body_spheroid.csv        - cell body point cloud (prolate spheroid)
%      chlamy_flagellum_ccw<sfx>.csv   - rest-shape centreline of the flagellum
%                                        on the +Y side (curls counter-clockwise
%                                        in the XY plane)
%      chlamy_flagellum_cw<sfx>.csv    - mirror image of the above (flagellum on
%                                        the -Y side, curls clockwise)
%    sfx = ''    for SHAPE_MODE = 'uniform'
%    sfx = '_gz' for SHAPE_MODE = 'gutierrez'
%
%    Frame: body-local, origin at the body centre.  +X is the swimming
%    direction (anterior = flagellar pole), the beat plane is the XY plane,
%    and Z is the beat-plane normal.
%
%    Cell body (Pradipta et al., arXiv:2509.06827, CC-125):
%      semi-major axis c = 5 um (along X), semi-minor axis a = 4 um.
%
%    Flagellum rest shape = STATIC (time-averaged) component of the
%    breaststroke beat.  The dynamic component is supplied by the
%    traveling-wave active moment of chlamydomonas/chlamy_moment_fun.m.
%
%    SHAPE_MODE = 'uniform'   (original baseline)
%      Planar circular arc of constant curvature C0 = 0.24 rad/um,
%      L = 11.7 um (Geyer et al., Curr. Biol. 26, 1098, 2016; Sartori et al.,
%      eLife 5, e13258, 2016; isolated axonemes).  Root half-angle 40 deg and
%      basal-body separation 0.5 um are ASSUMED.
%
%    SHAPE_MODE = 'gutierrez' (measured in-cell mean shape)
%      Curvature profile kappa(s) of the time-averaged tangent angle of
%      swimming cells in water, measured by Qin et al., Sci. Rep. 5, 9190
%      (2015) and re-analysed by Gutierrez et al., Phil. Trans. R. Soc. A 383,
%      20240505 (2025); data: Dryad doi:10.5061/dryad.sbcc2frj6.  The table
%      KL below is the mean over the three 1.00 cP cells (A, B, C) of
%      kappa_k * L_cell [rad per flagellar length], with kappa_k = d(phibar)/ds
%      sampled at s/L = k/31, k = 1..30 (phibar = time mean of the measured
%      tangent angle in the body frame).  Its integral gives a total turning
%      of about 124 deg: strongly bent over the first ~2 um (kappa ~ 1 rad/um),
%      a plateau ~0.25 rad/um up to ~0.4 L, then a monotone decrease to a
%      nearly straight distal part.  The profile is scaled to L_flag so that
%      the turning angle is preserved.
%      The basal-body position follows the data definition of the attachment
%      angle: the site is at the polar angle ATCH_DEG from the cell axis,
%      (x, y) = (c cos(atch), a sin(atch)).  The root tangent half-angle BETA
%      (angle between the root tangent and the cell axis, pointing away from
%      the axis) is taken equal to ATCH_DEG; the measured value of the three
%      1.00 cP cells is 9, 32 and 2 deg (mean 14 deg), with a large beat-to-
%      beat scatter.
%
%    Each file starts at the origin and its first chord already points along
%    the intended root tangent (+-BETA from the cell axis).  The solver
%    therefore applies a (near-)identity alignment rotation, so the beat
%    plane stays exactly the XY plane and the root director frame is
%    d1 = in-plane left normal, d2 = +Z (verified by
%    chlamydomonas/chlamy_static_check.m).
%
%  Usage:
%    >> cd geometry
%    >> generate_chlamydomonas_geometry
%
%  Outputs:
%    chlamy_body_spheroid.csv, chlamy_flagellum_ccw<sfx>.csv,
%    chlamy_flagellum_cw<sfx>.csv
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-06
%  Updated: 2026-11-09  added SHAPE_MODE = 'gutierrez' (measured mean shape)

close all; clear; clc;

out_dir = fileparts(mfilename('fullpath'));
if isempty(out_dir), out_dir = pwd; end
proj_root = fileparts(out_dir);
addpath(genpath(fullfile(proj_root, 'src')));
fprintf('Writing Chlamydomonas geometry files to: %s\n', out_dir);

%% --- User-adjustable parameters ---------------------------------------
shape_mode = 'gutierrez';   % 'uniform' | 'gutierrez'

c_body    = 5.0e-6;     % [m] semi-major axis (along +X)
a_body    = 4.0e-6;     % [m] semi-minor axis
N_body    = 400;        % [-] body surface points; must be one of the cached
                        %     sphere sets (e.g. 200, 300, 400, 600, 800)
L_flag    = 11.7e-6;    % [m] flagellum length (Geyer 2016 WT axoneme)
ds        = 0.1e-6;     % [m] source chord length (solver re-samples)

switch shape_mode
    case 'uniform'
        C0       = 0.24e6;   % [1/m] static curvature (0.24 rad/um)
        beta_deg = 40;       % [deg] root half-angle from the cell axis (ASSUMED)
        sep_bb   = 0.5e-6;   % [m] basal-body separation (ASSUMED)
        file_sfx = '';
    case 'gutierrez'
        atch_deg = 13;       % [deg] attachment-site polar angle from the cell axis
        beta_deg = 13;       % [deg] root half-angle from the cell axis
        file_sfx = '_gz';
        % mean of kappa*L [rad per L] over the three 1.00 cP cells, s/L = k/31
        KL = [ 6.224, 11.341,  8.309,  4.222,  2.901,  2.678,  2.644,  2.730, ...
               2.720,  2.533,  2.342,  2.166,  1.958,  1.719,  1.526,  1.366, ...
               1.160,  1.058,  1.131,  1.127,  0.961,  0.767,  0.665,  0.612, ...
               0.529,  0.382,  0.186, -0.158, -0.472, -0.485];
        sKL = (1:numel(KL)) / (numel(KL) + 1);
    otherwise
        error('generate_chlamydomonas_geometry:ShapeMode', ...
            'Unknown shape_mode "%s".', shape_mode);
end

%% --- 1) Body: prolate spheroid from the cached Riesz sphere sample -----
V = sphere_particle_sample_cached(N_body);         % [N x 3] unit sphere
body_pts = [c_body * V(:, 3), a_body * V(:, 1), a_body * V(:, 2)];

D2 = sum((reshape(body_pts, [N_body, 1, 3]) - ...
          reshape(body_pts, [1, N_body, 3])).^2, 3);
D2(1:N_body+1:end) = Inf;
nn = sqrt(min(D2, [], 2));
fprintf('  Body: %d points, h_min = %.3f um, h_mean = %.3f um\n', ...
    N_body, min(nn)*1e6, mean(nn)*1e6);
write_csv(fullfile(out_dir, 'chlamy_body_spheroid.csv'), body_pts, ...
    sprintf(['Chlamydomonas body: prolate spheroid, c = %.2f um (X), ', ...
             'a = %.2f um, %d points (scaled cached Riesz sphere sample)'], ...
             c_body*1e6, a_body*1e6, N_body));

%% --- 2) Flagellum rest shape --------------------------------------------
N_seg = max(round(L_flag / ds), 4);
s = linspace(0, L_flag, N_seg + 1)';
beta = deg2rad(beta_deg);

switch shape_mode
    case 'uniform'
        % Unrotated arc, start tangent +X, curling counter-clockwise (+C0)
        x0 = sin(C0 * s) / C0;
        y0 = (1 - cos(C0 * s)) / C0;
        descr = sprintf(['planar arc, C0 = %.3f rad/um, L = %.2f um, ', ...
                         'root half-angle %.1f deg, concave toward the body'], ...
                         C0*1e-6, L_flag*1e6, beta_deg);
        turn_deg = rad2deg(C0 * L_flag);
    case 'gutierrez'
        % kappa(s) [1/m] from the tabulated kappa*L profile (end values held)
        kap = interp1([0, sKL, 1], [KL(1), KL, KL(end)], s / L_flag, 'pchip') / L_flag;
        th  = cumtrapz(s, kap);                          % tangent angle from the root
        x0  = cumtrapz(s, cos(th));
        y0  = cumtrapz(s, sin(th));
        descr = sprintf(['measured mean shape (Qin 2015 / Gutierrez 2025, 1 cP), ', ...
                         'L = %.2f um, total turning %.1f deg, ', ...
                         'root half-angle %.1f deg, concave toward the body'], ...
                         L_flag*1e6, rad2deg(th(end)), beta_deg);
        turn_deg = rad2deg(th(end));
end

% +Y flagellum: rotate by +beta  (tangent (cos b,  sin b, 0), CCW curl)
R_p = [cos(beta), -sin(beta); sin(beta), cos(beta)];
P_ccw = [x0, y0] * R_p.';
% -Y flagellum: mirror of the above about the XZ plane
P_cw = [P_ccw(:, 1), -P_ccw(:, 2)];

write_csv(fullfile(out_dir, ['chlamy_flagellum_ccw', file_sfx, '.csv']), ...
    [P_ccw, zeros(N_seg + 1, 1)], ['Chlamydomonas flagellum (+Y side): ', descr]);
write_csv(fullfile(out_dir, ['chlamy_flagellum_cw', file_sfx, '.csv']), ...
    [P_cw, zeros(N_seg + 1, 1)], ...
    ['Chlamydomonas flagellum (-Y side), mirror image: ', descr]);

fprintf('  Flagellum (%s): L = %.2f um, total turning = %.1f deg, %d segments\n', ...
    shape_mode, L_flag*1e6, turn_deg, N_seg);

%% --- 3) Attachment sites for the JSON config ---------------------------
switch shape_mode
    case 'uniform'
        y_site = sep_bb / 2;
        x_site = c_body * sqrt(1 - (y_site / a_body)^2);
    case 'gutierrez'
        x_site = c_body * cosd(atch_deg);
        y_site = a_body * sind(atch_deg);
end
fprintf('\nAttachment sites for custom_geometry.attachment_sites:\n');
fprintf('  +Y: position = [%.6e, %+.6e, 0], direction = [%.6f, %+.6f, 0]\n', ...
    x_site, y_site, cos(beta), sin(beta));
fprintf('  -Y: position = [%.6e, %+.6e, 0], direction = [%.6f, %+.6f, 0]\n', ...
    x_site, -y_site, cos(beta), -sin(beta));
fprintf('  basal-body separation = %.2f um\n', 2*y_site*1e6);

%% --- Bending stiffness helper ------------------------------------------
EI_target = 840e-24;     % [N m^2] 840 pN um^2 (Xu et al. 2016, CC-124)
r_pili    = 0.125e-6;    % [m] hydrodynamic radius (axoneme + membrane)
E_eff     = 4 * EI_target / (pi * r_pili^4);
fprintf('\nMaterial: r_pili = %.3f um, E = %.4e Pa gives EI = %.1f pN um^2\n', ...
    r_pili*1e6, E_eff, E_eff * pi * r_pili^4 / 4 * 1e24);

% ---- local helper ----
function write_csv(path, pts, descr)
%WRITE_CSV  Write an N x 3 point cloud with a 2-line comment header.
    [fid, msg] = fopen(path, 'w', 'n', 'UTF-8');
    if fid < 0
        error('generate_chlamydomonas_geometry:WriteFail', ...
            'Failed to open %s: %s', path, msg);
    end
    fprintf(fid, '# %s\n', descr);
    fprintf(fid, '# Columns: X [m], Y [m], Z [m]\n');
    fprintf(fid, '%.9e,%.9e,%.9e\n', pts.');
    fclose(fid);
    fprintf('  Wrote %d points -> %s\n', size(pts, 1), path);
end
