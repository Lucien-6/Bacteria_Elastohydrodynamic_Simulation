function generate_chlamydomonas_geometry()
%GENERATE_CHLAMYDOMONAS_GEOMETRY  Generate measured Chlamydomonas CSV geometry
%
%  Description:
%    Writes the biflagellate Chlamydomonas reinhardtii (wild type CC-124 /
%    CC-125) geometry CSV files beside this script for
%    configs/config_chlamy_biflagellate_gz.json:
%      chlamy_body_spheroid.csv      - prolate spheroid body point cloud
%      chlamy_flagellum_ccw_gz.csv   - +Y flagellum measured mean shape
%      chlamy_flagellum_cw_gz.csv    - -Y flagellum mirror image
%
%    Frame: body-local, origin at the body centre.  +X is the swimming
%    direction (anterior = flagellar pole), the reference beat plane is XY,
%    and +Z is its normal.
%
%    Cell body (Pradipta et al., arXiv:2509.06827, CC-125):
%      semi-major axis c = 5 um (along X), semi-minor axis a = 4 um,
%      N_body = 400 cached Riesz sphere sample points.
%
%    Flagellum rest shape (measured in-cell time-averaged shape in water):
%      Curvature profile kappa(s) averaged over the three 1.00 cP cells
%      (A, B, C) of Qin et al., Sci. Rep. 5, 9190 (2015) / Gutierrez et al.,
%      Phil. Trans. R. Soc. A 383, 20240505 (2025), Dryad
%      doi:10.5061/dryad.sbcc2frj6.  L_flag = 11.7 um, attachment polar
%      angle = root half-angle = 13 deg, total static turning ~ 124 deg.
%      The dynamic 2-D/3-D beat is supplied by chlamy_table_moment_fun and
%      chlamy_make_perp_envelope.
%
%  Usage:
%    generate_chlamydomonas_geometry()
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-06
%  Updated: 2026-11-24  function form, drop obsolete uniform circular-arc mode

out_dir = fileparts(mfilename('fullpath'));
if isempty(out_dir), out_dir = pwd; end
proj_root = fileparts(out_dir);
addpath(genpath(fullfile(proj_root, 'src')));
fprintf('Writing Chlamydomonas geometry files to: %s\n', out_dir);

%% --- Physical geometry parameters --------------------------------------
c_body   = 5.0e-6;     % [m] semi-major axis (along +X)
a_body   = 4.0e-6;     % [m] semi-minor axis
N_body   = 400;        % [-] cached Riesz sphere set (200,300,400,600,800)
L_flag   = 11.7e-6;    % [m] flagellum length (WT axoneme)
ds       = 0.1e-6;     % [m] source chord length (solver re-samples)
atch_deg = 13;         % [deg] attachment-site polar angle from cell axis
beta_deg = 13;         % [deg] root half-angle from cell axis

% Mean kappa*L [rad per L] over the three 1.00 cP cells at s/L = k/31
KL = [ 6.224, 11.341,  8.309,  4.222,  2.901,  2.678,  2.644,  2.730, ...
       2.720,  2.533,  2.342,  2.166,  1.958,  1.719,  1.526,  1.366, ...
       1.160,  1.058,  1.131,  1.127,  0.961,  0.767,  0.665,  0.612, ...
       0.529,  0.382,  0.186, -0.158, -0.472, -0.485];
sKL = (1:numel(KL)) / (numel(KL) + 1);

%% --- 1) Body: prolate spheroid from cached Riesz sphere sample ---------
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

%% --- 2) Flagellum rest shape (measured mean curvature) -----------------
N_seg = max(round(L_flag / ds), 4);
s     = linspace(0, L_flag, N_seg + 1).';
beta  = deg2rad(beta_deg);

kap = interp1([0, sKL, 1], [KL(1), KL, KL(end)], s / L_flag, 'pchip') / L_flag;
th  = cumtrapz(s, kap);
x0  = cumtrapz(s, cos(th));
y0  = cumtrapz(s, sin(th));
descr = sprintf(['measured mean shape (Qin 2015 / Gutierrez 2025, 1 cP), ', ...
                 'L = %.2f um, total turning %.1f deg, ', ...
                 'root half-angle %.1f deg, concave toward the body'], ...
                 L_flag*1e6, rad2deg(th(end)), beta_deg);

% +Y flagellum rotated by +beta; -Y flagellum is its XZ mirror image
R_p   = [cos(beta), -sin(beta); sin(beta), cos(beta)];
P_ccw = [x0, y0] * R_p.';
P_cw  = [P_ccw(:, 1), -P_ccw(:, 2)];

write_csv(fullfile(out_dir, 'chlamy_flagellum_ccw_gz.csv'), ...
    [P_ccw, zeros(N_seg + 1, 1)], ['Chlamydomonas flagellum (+Y side): ', descr]);
write_csv(fullfile(out_dir, 'chlamy_flagellum_cw_gz.csv'), ...
    [P_cw, zeros(N_seg + 1, 1)], ...
    ['Chlamydomonas flagellum (-Y side), mirror image: ', descr]);

fprintf('  Flagellum (gutierrez): L = %.2f um, total turning = %.1f deg, %d segments\n', ...
    L_flag*1e6, rad2deg(th(end)), N_seg);

%% --- 3) Attachment sites & stiffness summary ---------------------------
x_site = c_body * cosd(atch_deg);
y_site = a_body * sind(atch_deg);
fprintf('\nAttachment sites for configs/config_chlamy_biflagellate_gz.json:\n');
fprintf('  +Y: position = [%.6e, %+.6e, 0], direction = [%.6f, %+.6f, 0]\n', ...
    x_site, y_site, cos(beta), sin(beta));
fprintf('  -Y: position = [%.6e, %+.6e, 0], direction = [%.6f, %+.6f, 0]\n', ...
    x_site, -y_site, cos(beta), -sin(beta));
fprintf('  basal-body separation = %.2f um\n', 2*y_site*1e6);

EI_target = 840e-24;     % [N m^2] 840 pN um^2 (Xu et al. 2016, CC-124)
r_pili    = 0.125e-6;    % [m] hydrodynamic radius (axoneme + membrane)
E_eff     = 4 * EI_target / (pi * r_pili^4);
fprintf('\nMaterial: r_pili = %.3f um, E = %.4e Pa gives EI = %.1f pN um^2\n', ...
    r_pili*1e6, E_eff, E_eff * pi * r_pili^4 / 4 * 1e24);
end


function write_csv(path, pts, descr)
%WRITE_CSV  Write an N x 3 point cloud with a 2-line comment header.
    [fid, msg] = fopen(path, 'w', 'n', 'UTF-8');
    if fid < 0
        error('generate_chlamydomonas_geometry:WriteFail', ...
            'Failed to open %s: %s', path, msg);
    end
    cleanup_fid = onCleanup(@() fclose(fid));
    fprintf(fid, '# %s\n', descr);
    fprintf(fid, '# Columns: X [m], Y [m], Z [m]\n');
    fprintf(fid, '%.9e,%.9e,%.9e\n', pts.');
    fprintf('  Wrote %d points -> %s\n', size(pts, 1), path);
end
