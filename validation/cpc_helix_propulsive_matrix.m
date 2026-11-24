function summary = cpc_helix_propulsive_matrix(varargin)
%CPC_HELIX_PROPULSIVE_MATRIX  Rigid-helix thrust, torque and drag
%
%  Description:
%    External hydrodynamic check for the CPC program paper.  A rigid
%    circular helix (constant radius; not the conical type-5 template) is
%    moved without translation or without rotation.  The project RSS kernel
%    supplies the three propulsive coefficients
%
%        F / (mu * omega * R^2),   T / (mu * omega * R^3),   D / (mu * U * R).
%
%    The pass line is a slender-body operator on the same centreline:
%    Johnson's prolate logarithm plus the finite-part Stokeslet.  The a^2
%    dipole of Rodenborn et al., Proc. Natl. Acad. Sci. U.S.A. 110, E338
%    (2013), is not added on neighbour segments, because that raw dipole
%    does not converge with segment length.  Gray-Hancock resistive-force
%    theory is stored as a negative control.  The naked helix is the
%    comparison.  One body-attached repeat at the middle pitch only reports
%    the shift and is not part of the pass line.
%
%    Geometry, in units of the helical radius R: baseline L/R = 20 and
%    lambda/R = 6.5 match the benchmark helical geometry of Rodenborn et al.
%    (2013), Fig. 7, with representative slender cross-section a/R = 1/16
%    evaluated consistently across RSS, Johnson SBT, and Gray-Hancock RFT.
%    Additional pitches lambda/R = 3 and 12 and extended length L/R = 40
%    bracket the pitch and aspect-ratio parameter space.
%
%    Segment density starts at 16 and 32 nodes per turn.  A coefficient
%    that still moves by more than 2 percent is repeated at 64 per turn.
%    beta = sqrt(2) is the reported radius.  beta = 1 and beta = sqrt(e)
%    are stored only for the middle pitch.
%
%  Syntax:
%    summary = cpc_helix_propulsive_matrix()
%    summary = cpc_helix_propulsive_matrix('OutDir', outdir)
%
%  Output:
%    summary - struct with .table and .passed.  The table is also written to
%              <OutDir>/helix_propulsive_matrix.csv.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-15

ip = inputParser;
addParameter(ip, 'OutDir', '');
parse(ip, varargin{:});

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(genpath(fullfile(root, 'src')));
addpath(fileparts(mfilename('fullpath')));

if isempty(ip.Results.OutDir)
    outdir = fullfile(root, 'output', 'cpc_checks');
else
    outdir = ip.Results.OutDir;
end
if ~isfolder(outdir)
    mkdir(outdir);
end

geom = local_geometries();
straight = straight_rod_check();
fprintf(['Johnson straight-rod perpendicular drag / Broersma = %.4f ', ...
    '(diagnostic; the helix pass line is the curved centreline)\n'], straight.ratio);

parts = cell(1, numel(geom));
prod_segs = zeros(1, numel(geom));
for g = 1:numel(geom)
    [prod_segs(g), parts{g}] = converge_geometry(geom(g));
    prod = parts{g}(parts{g}.is_production, :);
    fprintf(['%s  segs/turn = %g  RSS/Johnson  F %.4f  T %.4f  D %.4f  ', ...
        'max |rel| %.3f  conv %.3f\n'], geom(g).name, prod.segs_per_turn, ...
        prod.F_rss / prod.F_johnson, prod.T_rss / prod.T_johnson, ...
        prod.D_rss / prod.D_johnson, prod.rel_johnson, prod.rel_conv);
end
rows = vertcat(parts{:});

mid = geom(strcmp({geom.name}, 'pitch_mid'));
n_mid = prod_segs(strcmp({geom.name}, 'pitch_mid'));
rows = [rows; coefficient_row(mid, n_mid, 1, false, false, NaN)];
rows = [rows; coefficient_row(mid, n_mid, sqrt(exp(1)), false, false, NaN)];
rows = [rows; coefficient_row(mid, n_mid, sqrt(2), true, false, NaN)];

naked = rows(rows.with_body == 0 & abs(rows.beta - sqrt(2)) < 1e-12 & rows.is_production, :);
passed = all(naked.rel_conv <= 0.02) && all(naked.rel_johnson <= 0.05);

summary = struct('table', rows, 'passed', passed, 'straight_rod_ratio', straight.ratio);
writetable(rows, fullfile(outdir, 'helix_propulsive_matrix.csv'));
fprintf('Helix check %s.  Table: %s\n', pass_word(passed), fullfile(outdir, 'helix_propulsive_matrix.csv'));
end


function geom = local_geometries()
geom = struct( ...
    'name', {'pitch_short', 'pitch_mid', 'pitch_long', 'length_long'}, ...
    'lambda_over_R', {3, 6.5, 12, 6.5}, ...
    'L_over_R', {20, 20, 20, 40}, ...
    'a_over_R', {1/16, 1/16, 1/16, 1/16});
end


function word = pass_word(flag)
if flag
    word = 'PASSED';
else
    word = 'FAILED';
end
end


function check = straight_rod_check()
L = 20;
a = 1/16;
n_seg = 160;
x = linspace(-L / 2, L / 2, n_seg + 1).';
X = [zeros(n_seg + 1, 2), x];
U = [ones(n_seg + 1, 1), zeros(n_seg + 1, 2)];
F = slender_body_load(X, U, a);
p_ar = L / (2 * a);
broersma = 4 * pi * L / (log(p_ar) + 0.839 + 0.185 / p_ar + 0.233 / p_ar^2);
check = struct('ratio', F(1) / broersma, 'broersma', broersma);
end


function [n_prod, detail] = converge_geometry(geom)
c16 = coefficient_row(geom, 16, sqrt(2), false, false, NaN);
c32 = coefficient_row(geom, 32, sqrt(2), false, false, NaN);
c32.rel_conv = rel_coeff(c16, c32);
if c32.rel_conv > 0.02
    c64 = coefficient_row(geom, 64, sqrt(2), false, true, NaN);
    c64.rel_conv = rel_coeff(c32, c64);
    n_prod = 64;
    detail = [c16; c32; c64];
else
    c32.is_production = true;
    n_prod = 32;
    detail = [c16; c32];
end
end


function value = rel_coeff(row_a, row_b)
value = max([abs(row_b.F_rss - row_a.F_rss) / abs(row_a.F_rss), ...
    abs(row_b.T_rss - row_a.T_rss) / abs(row_a.T_rss), ...
    abs(row_b.D_rss - row_a.D_rss) / abs(row_a.D_rss)]);
end


function row = coefficient_row(geom, segs_per_turn, beta, with_body, is_production, rel_conv)
[X, ~] = helix_centerline(geom.lambda_over_R, geom.L_over_R, segs_per_turn);
eps_line = beta * geom.a_over_R;
[F_rot, T_rot] = project_motion(X, eps_line, 'rot', with_body);
D_trans = project_motion(X, eps_line, 'trans', with_body);
[F_johnson, T_johnson] = reference_motion(X, geom.a_over_R, 'rot', @slender_body_load);
D_johnson = reference_motion(X, geom.a_over_R, 'trans', @slender_body_load);
[F_rft, T_rft] = rft_motion(X, geom, 'rot');
D_rft = rft_motion(X, geom, 'trans');
rel_j = max([abs(F_rot - F_johnson) / abs(F_johnson), ...
    abs(T_rot - T_johnson) / abs(T_johnson), ...
    abs(D_trans - D_johnson) / abs(D_johnson)]);
row = table(string(geom.name), geom.lambda_over_R, geom.L_over_R, geom.a_over_R, ...
    segs_per_turn, beta, with_body, is_production, ...
    F_rot, T_rot, D_trans, F_johnson, T_johnson, D_johnson, F_rft, T_rft, D_rft, ...
    rel_conv, rel_j, ...
    'VariableNames', {'name', 'lambda_over_R', 'L_over_R', 'a_over_R', ...
    'segs_per_turn', 'beta', 'with_body', 'is_production', ...
    'F_rss', 'T_rss', 'D_rss', 'F_johnson', 'T_johnson', 'D_johnson', ...
    'F_rft', 'T_rft', 'D_rft', 'rel_conv', 'rel_johnson'});
end


function [F, T] = project_motion(X, eps_line, mode, with_body)
U = rigid_velocity(X, mode);
if ~with_body
    [force, torque] = rss_line_load(X, U, eps_line);
else
    radius = 0.2;
    X_body = fibonacci_sphere(200, radius);
    [force, torque] = rss_body_line_load(X_body, X, rigid_velocity(X_body, mode), U, ...
        mean_nearest(X_body) / 2.6667, eps_line);
end
F = force(3);
T = torque(3);
end


function [F, T] = reference_motion(X, a_over_R, mode, loader)
U = rigid_velocity(X, mode);
[force, torque] = loader(X, U, a_over_R);
F = force(3);
T = torque(3);
if strcmp(mode, 'trans')
    T = 0;
end
end


function [F, T] = rft_motion(X, geom, mode)
U = rigid_velocity(X, mode);
[force, torque] = rft_helix_load(X, U, geom.lambda_over_R, geom.a_over_R);
F = force(3);
T = torque(3);
if strcmp(mode, 'trans')
    T = 0;
end
end


function [X, n_seg] = helix_centerline(lambda_over_R, L_over_R, segs_per_turn)
n_turns = L_over_R / lambda_over_R;
n_seg = max(8, round(segs_per_turn * n_turns));
theta = linspace(0, 2 * pi * n_turns, n_seg + 1).';
X = [cos(theta), sin(theta), (lambda_over_R / (2 * pi)) * theta];
end


function U = rigid_velocity(X, mode)
n = size(X, 1);
if strcmp(mode, 'rot')
    U = [-X(:, 2), X(:, 1), zeros(n, 1)];
else
    U = [zeros(n, 2), ones(n, 1)];
end
end


function X = fibonacci_sphere(n_point, radius)
i = (0:n_point - 1).' + 0.5;
phi = acos(1 - 2 * i / n_point);
theta = pi * (1 + sqrt(5)) * i;
X = radius * [sin(phi) .* cos(theta), sin(phi) .* sin(theta), cos(phi)];
end


function h = mean_nearest(X)
n = size(X, 1);
nearest = inf(n, 1);
for i = 1:n
    d = sqrt(sum((X - X(i, :)).^2, 2));
    d(i) = inf;
    nearest(i) = min(d);
end
h = mean(nearest);
end
