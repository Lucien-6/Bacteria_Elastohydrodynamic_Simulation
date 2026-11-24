function tab = chlamy_inverse_moment(mat_file, opts)
%CHLAMY_INVERSE_MOMENT  Active internal moment q(s,t) that reproduces a measured gait
%
%  Description:
%    Inverse problem of the elastohydrodynamic model.  A measured flagellar
%    gait (Qin et al. 2015 / Gutierrez et al. 2025 data file) is fitted with
%    a mean tangent-angle profile plus NHARM Fourier harmonics, mapped to the
%    frame of the project geometry (X = swimming axis, flagellum 1 on +Y,
%    flagellum 2 its mirror image) and imposed as a prescribed state
%    Z(t) = [X_cm; q; theta = pi/2, phi = tangent angle, psi = 0] together
%    with its exact rate vector
%
%        zdot = [V; Omega; 0; dphi/dt + Omega_z; 0].
%
%    The unknown rigid-body velocities (V_x, V_y, Omega_z) follow from the
%    force- and torque-balance rows of the project's own linear system
%    (build_coupled_lin_sys), so the hydrodynamics, regularisation and
%    discretisation are exactly those of the forward solver.  The rod rows of
%    the forward system read  lin_sys * zdot = R,  R = kappa_dev - M_act;
%    hence the active moment that makes the prescribed motion an exact
%    solution of the discrete forward model is
%
%        M_act(:, t) = kappa_dev(Z(t)) - lin_sys(Z(t)) * zdot(t)      (rod rows).
%
%    This is the INTERNAL moment across the cross-section at node j (rows of
%    the 'custom' active moment, in units of EI/L, d2 component).  It contains
%    the elastic part EI*(kappa - kappa0) and the moment of the hydrodynamic
%    load beyond node j.  Because the clamped root forces the first segment
%    to be fixed in the body frame, the measured oscillation of the tangent
%    angle is ramped from zero at the first segment to its full measured
%    value at RAMP_S (relative arc length).
%
%    q(s,t) is stored as a Fourier series in time (NQ harmonics) per node
%    and flagellum, to be used by chlamy_table_moment_fun / chlamy_run_case
%    (override field .moment_table).
%
%  Syntax:
%    tab = chlamy_inverse_moment(mat_file)
%    tab = chlamy_inverse_moment(mat_file, opts)
%
%  Inputs:
%    mat_file - Path of a data .mat file (absolute or relative to the root)
%    opts     - (Optional) struct:
%                 .config    JSON (default configs/config_chlamy_biflagellate_gz.json)
%                 .f_Hz      beat frequency (default: measured value)
%                 .nharm     harmonics of the gait fit (default 2)
%                 .ramp_s    end of the proximal ramp, s/L (default 0.15)
%                 .n_t       time samples per period (default 48)
%                 .nq        Fourier harmonics kept for q (default 8)
%                 .verify    true: evaluate ode_coupled_rhs with the tabulated
%                            moment and compare with the prescribed rates
%                            (default true)
%                 .save      path of a .mat file to store TAB (default none)
%                 .quiet     suppress printing
%
%  Output:
%    tab - struct: .f_Hz, .T_0, .w_star, .nq, .N_seg, .coef{k} [N_seg x (2nq+1)]
%          (columns: mean, cos 1, sin 1, cos 2, ...; dimensionless moment,
%          units EI/L_k), .q_samples{k} [N_seg x n_t], .t_s, .U_t_um_s
%          (axial swimming velocity of the prescribed gait), .U_mean_um_s,
%          .Omega_z_t, .Q_peak_Nm (peak |q| of each flagellum), .fit_err
%          (relative Fourier-fit error), .res_rows (max residual of the body
%          and root rows, relative), .off_plane (max |rows 1,3| / max |row 2|),
%          .L_um, .th_mean_deg (mean tangent angle profile of flagellum 1),
%          .file.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-09

proj_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(proj_root);
addpath(genpath(fullfile(proj_root, 'src')));
addpath(fileparts(mfilename('fullpath')));
if nargin < 2, opts = struct(); end
cfg_file = getopt(opts, 'config', 'configs/config_chlamy_biflagellate_gz.json');
nharm    = getopt(opts, 'nharm', 2);
ramp_s   = getopt(opts, 'ramp_s', 0.15);
n_t      = getopt(opts, 'n_t', 48);
nq       = getopt(opts, 'nq', 8);
quiet    = getopt(opts, 'quiet', false);
do_ver   = getopt(opts, 'verify', true);

[params, Z0] = chlamy_build_params(cfg_file, @(varargin) []);
assert(params.M_pili == 2 && all(params.N_segments == params.N_segments(1)), ...
    'chlamy_inverse_moment:Geometry', 'Two flagella with equal segment numbers are required.');
N_seg = params.N_segments(1);
L_0 = params.scales.L_0;  T_0 = params.scales.T_0;
L_um = params.ell(1) * L_0 * 1e6;
n_state = 6 + 3 * sum(params.N_segments);
assert(all(abs(Z0(4:7) - [1; 0; 0; 0]) < 1e-12), 'chlamy_inverse_moment:Orientation', ...
    'Identity initial orientation expected.');
assert(params.Eg == 0, 'chlamy_inverse_moment:Gravity', 'Gravity must be off.');

%% ---- measured gait: mean + harmonics of the tangent angle -----------------
if ~isfile(mat_file), mat_file = fullfile(proj_root, mat_file); end
m   = load(mat_file);
phi = unwrap(m.phi, [], 1);                    % [31 x Nt] tangent angle (data frame)
[Nsd, Nt] = size(phi);
dt_d = m.time_scale(1);  t_d = (0:Nt-1) * dt_d;
x  = phi(round(Nsd/2), :) - mean(phi(round(Nsd/2), :));
fr = 30 : 0.05 : 90;
amp = arrayfun(@(v) abs(sum(x .* exp(-2i * pi * v * t_d))), fr);
[~, ib] = max(amp);  f_meas = fr(ib);
f_Hz = getopt(opts, 'f_Hz', f_meas);

B = zeros(Nt, 2*nharm + 1);
B(:, 1) = 1;
for n = 1:nharm
    B(:, 2*n)   = cos(2*pi*n*f_meas*t_d.');
    B(:, 2*n+1) = sin(2*pi*n*f_meas*t_d.');
end
coef = B \ phi.';                              % [(2n+1) x Nsd]
s_d  = ((1:Nsd) - 0.5) / Nsd;
s_m  = ((1:N_seg).' - 0.5) / N_seg;
coefN = zeros(size(coef, 1), N_seg);
for ic = 1:size(coef, 1)
    c = interp1(s_d, coef(ic, :), s_m, 'linear', 'extrap');
    c(s_m < s_d(1))   = coef(ic, 1);
    c(s_m > s_d(end)) = coef(ic, end);
    coefN(ic, :) = c;
end

% model frame: theta_t = pi/2 - phi; mean aligned to the rest root tangent
d3ref  = params.directors_ref{1}.d3(1, :);
th_rest1 = atan2(d3ref(2), d3ref(1));
th_mean  = pi/2 - coefN(1, :).';
th_mean  = th_mean + (th_rest1 - th_mean(1));
x_r = min(max((s_m - s_m(1)) / (ramp_s - s_m(1)), 0), 1);
ramp = x_r.^2 .* (3 - 2 * x_r);                % exactly 0 at the clamped segment

w = 2 * pi * f_Hz;
T_b = 1 / f_Hz;
tt  = (0:n_t-1) * T_b / n_t;

%% ---- inversion at each time sample ---------------------------------------
N_k = N_seg;
rows_pili = cell(1, 2);
off = 6;
for k = 1:2
    rows_pili{k} = off + (1:3*N_k);
    off = off + 3*N_k;
end
rc = 'clamped';
q2  = cell(1, 2);  q2{1} = zeros(N_k, n_t);  q2{2} = zeros(N_k, n_t);
q13 = 0;  res_rows = 0;  res_root = 0;
U_t = zeros(n_t, 1);  Om_t = zeros(n_t, 1);
Zs = cell(1, n_t);  zd_s = cell(1, n_t);
for it = 1:n_t
    th = th_mean;  dth = zeros(N_seg, 1);
    for n = 1:nharm
        c = cos(n * w * tt(it));  s = sin(n * w * tt(it));
        th  = th  - (coefN(2*n, :).' * c + coefN(2*n+1, :).' * s) .* ramp;
        dth = dth - n * w * (-coefN(2*n, :).' * s + coefN(2*n+1, :).' * c) .* ramp;
    end
    dth = dth * T_0;                           % dimensionless time derivative
    Z = Z0;
    z_a = zeros(n_state, 1);  z_w = zeros(n_state, 1);
    z_w(6) = 1;
    col = 6;  idx = 8;
    for k = 1:2
        sg = (k == 1) * 1 + (k == 2) * (-1);   % mirror image
        Z(idx           : idx + N_k - 1)     = pi/2;
        Z(idx + N_k     : idx + 2*N_k - 1)   = sg * th;
        Z(idx + 2*N_k   : idx + 3*N_k - 1)   = 0;
        z_a(col + N_k + (1:N_k)) = sg * dth;
        z_w(col + N_k + (1:N_k)) = 1;
        idx = idx + 3*N_k;  col = col + 3*N_k;
    end
    Z = enforce_clamped_geometry(Z, params);
    [lin_sys, geom] = build_coupled_lin_sys(Z, params);
    z_x = zeros(n_state, 1);  z_x(1) = 1;
    z_y = zeros(n_state, 1);  z_y(2) = 1;
    ri = [1 2 6];
    A3 = lin_sys(ri, :) * [z_x, z_y, z_w];
    u  = A3 \ (-lin_sys(ri, :) * z_a);
    zdot = z_a + [z_x, z_y, z_w] * u;
    Rrow = lin_sys * zdot;
    % body rows must vanish (no gravity, no Brownian force)
    res_rows = max(res_rows, max(abs(Rrow(1:6))) / max(abs(Rrow(rows_pili{1}(:))), [], 'all'));
    for k = 1:2
        [k1, k2, k3] = compute_pili_elastic_moments(geom.d1_basis{k}, geom.d2_basis{k}, ...
            geom.d3_basis{k}, N_k, params.sigma, params.kb0{k}, rc);
        kap = reshape([k1.'; k2.'; k3.'], [], 1);
        Mact = reshape(kap - Rrow(rows_pili{k}), 3, N_k);
        q2{k}(:, it) = Mact(2, :).';
        q13 = max(q13, max(abs(Mact([1 3], 2:end)), [], 'all'));
        res_root = max(res_root, max(abs(Rrow(rows_pili{k}(1:3)))));   % root BC rows
    end
    U_t(it) = zdot(1) * L_0 / T_0 * 1e6;
    Om_t(it) = zdot(6) / T_0;
    Zs{it} = Z;  zd_s{it} = zdot;
end

%% ---- Fourier series in time ------------------------------------------------
Bq = zeros(n_t, 2*nq + 1);
Bq(:, 1) = 1;
for n = 1:nq
    Bq(:, 2*n)   = cos(n * w * tt.');
    Bq(:, 2*n+1) = sin(n * w * tt.');
end
tab = struct('f_Hz', f_Hz, 'T_0', T_0, 'w_star', w * T_0, 'nq', nq, 'N_seg', N_seg);
tab.coef = cell(1, 2);  tab.q_samples = q2;
fit_err = 0;  Qpk = zeros(1, 2);
EI = params.scales.EI;
for k = 1:2
    cq = (Bq \ q2{k}.').';                     % [N_k x (2nq+1)]
    tab.coef{k} = cq;
    fit_err = max(fit_err, norm(q2{k} - cq * Bq.', 'fro') / norm(q2{k}, 'fro'));
    Qpk(k) = max(abs(q2{k}), [], 'all') * EI / (params.ell(k) * L_0);
end
tab.Q_peak_Nm = Qpk;
tab.fit_err = fit_err;
tab.t_s = tt;  tab.U_t_um_s = U_t;  tab.U_mean_um_s = mean(U_t);
tab.Omega_z_t = Om_t;
tab.res_rows = res_rows;  tab.res_root = res_root;
tab.off_plane = q13 / max(abs([q2{1}(:); q2{2}(:)]));
tab.L_um = L_um;  tab.th_mean_deg = rad2deg(th_mean);
tab.file = mat_file;  tab.f_meas_Hz = f_meas;  tab.ramp_s = ramp_s;  tab.nharm = nharm;
tab.config = cfg_file;

%% ---- closed-loop identity: forward RHS with the tabulated moment ----------
if do_ver
    pv = params;
    for k = 1:2
        pv.am_cfg{k} = struct('type', 'custom', 'fun', chlamy_table_moment_fun(tab, k));
    end
    err_V = 0;  err_phi = 0;
    for it = 1:n_t
        dZ = ode_coupled_rhs(tt(it) / T_0, Zs{it}, pv);
        zd = zd_s{it};
        err_V = max(err_V, abs(dZ(1) - zd(1)) / max(abs(U_t(it)) * T_0 / L_0 / 1e6, 1e-3));
        idx = 7;  col = 6;
        for k = 1:2
            ph_r = dZ(idx + N_k + (1:N_k));
            ph_t = zd(col + N_k + (1:N_k));
            err_phi = max(err_phi, norm(ph_r - ph_t) / norm(ph_t));
            idx = idx + 3*N_k;  col = col + 3*N_k;
        end
    end
    tab.verify_err_phi = err_phi;  tab.verify_err_V = err_V;
end

if ~isempty(getopt(opts, 'save', []))
    sp = opts.save;
    if ~isfolder(fileparts(sp)), mkdir(fileparts(sp)); end
    save(sp, 'tab');
end

if ~quiet
    [~, nm] = fileparts(mat_file);
    fprintf(['[INV] %-28s f = %.1f Hz | gait mean turning %.0f deg | prescribed-gait U = %.1f um/s ', ...
             '(U/(fL) = %.3f) | peak q = %.2e / %.2e N m\n'], nm, f_Hz, ...
        th_mean(end) * 180/pi - th_mean(1) * 180/pi, tab.U_mean_um_s, ...
        tab.U_mean_um_s / (f_Hz * L_um), Qpk(1), Qpk(2));
    fprintf('[INV] Fourier fit err %.2e | body-row residual %.2e | off-plane rows %.2e', ...
        fit_err, res_rows, tab.off_plane);
    if do_ver
        fprintf(' | forward-RHS check: dphi err %.2e, V err %.2e', tab.verify_err_phi, tab.verify_err_V);
    end
    fprintf('\n');
end
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
