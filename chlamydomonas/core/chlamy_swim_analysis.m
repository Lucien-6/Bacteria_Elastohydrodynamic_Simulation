function s = chlamy_swim_analysis(results, opts)
%CHLAMY_SWIM_ANALYSIS  Swimming kinematics and 3-D geometric post-checks
%
%  Description:
%    Evaluates, over an integer number of beat periods after a start-up
%    transient:
%      * body-frame mean translational velocity V_body (averaged step-by-step
%        in the instantaneous body frame, exact even over multi-turn helical
%        trajectories > 1 s) and lab-frame net displacement velocity V_lab,
%      * unwrapped body-frame angular velocity Omega_body accumulated from
%        frame-to-frame SO(3) logarithms log(R_m' * R_{m+1}) (exact for any
%        number of axial revolutions, avoiding the < pi folding error),
%      * axial spin frequency spin_Hz = Omega_body(1) / (2*pi) and total
%        revolutions over the analysis window,
%      * screw/helical trajectory parameters (pitch, radius, period, axis
%        tilt, handedness screw_sign),
%      * flagellum-body spheroid clearance (f_body_min > 1 required) and
%        minimum inter-flagellar distance,
%      * body-frame tip 3-D orbit diagnostics and instantaneous flagellar
%        centerline 3-D non-planarity.
%
%  Syntax:
%    s = chlamy_swim_analysis(results)
%    s = chlamy_swim_analysis(results, opts)
%
%  Inputs:
%    results - Struct from run_simulation
%    opts    - (Optional) struct:
%                .f_Hz    reference beat frequency [Hz]    (default 57.75)
%                .n_skip  beats discarded at the start     (default 1)
%                .c_um    spheroid semi-axis along X [um]  (default 5)
%                .a_um    spheroid semi-axes along Y,Z     (default 4)
%
%  Output:
%    s - Struct of swimming, helical, clearance, and 3-D waveform metrics
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-06
%  Updated: 2026-11-24  unwrapped step-by-step SO(3) rotation & body-frame
%                       velocity integration (exact for multi-turn 2 s runs)

arguments
    results (1,1) struct
    opts    (1,1) struct = struct()
end

proj_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(proj_root, 'src')));

f_Hz   = getopt(opts, 'f_Hz', 57.75);
n_skip = getopt(opts, 'n_skip', 1);
c_um   = getopt(opts, 'c_um', 5);
a_um   = getopt(opts, 'a_um', 4);
T_b    = 1 / f_Hz;

t  = results.t(:);
NT = numel(t);
n_beats = floor((t(end) - n_skip * T_b) / T_b + 1e-9);
if n_beats < 1
    error('chlamy_swim_analysis:TooShort', 'Run is shorter than %d beat(s).', n_skip + 1);
end
t0 = n_skip * T_b;
t1 = t0 + n_beats * T_b;
i0 = find(t >= t0 - 1e-9, 1, 'first');
i1 = find(t <= t1 + 1e-9, 1, 'last');
iw = i0:i1;
Tw = t(i1) - t(i0);

Rm = zeros(3, 3, NT);
for it = 1:NT
    Rm(:, :, it) = quat_to_rotmat(results.q(:, it));
end

%% Unwrapped body-frame translation and SO(3) rotation (exact for > 1 turn)
Xc = results.X_cm_um;
V_lab = (Xc(:, i1) - Xc(:, i0)) / Tw;

dX_body_acc = zeros(3, 1);
dTh_body_acc = zeros(3, 1);
for m = i0:(i1 - 1)
    R_step = Rm(:, :, m).' * Rm(:, :, m + 1);
    dth = rotmat_log(R_step);
    dTh_body_acc = dTh_body_acc + dth;
    % Mid-step body frame for unbiased trapezoidal displacement projection
    R_half = Rm(:, :, m) * rotmat_exp(0.5 * dth);
    dX_body_acc = dX_body_acc + R_half.' * (Xc(:, m + 1) - Xc(:, m));
end

V_body = dX_body_acc / Tw;                       % [um/s] mean body-frame velocity
Om_b   = dTh_body_acc / Tw;                      % [rad/s] unwrapped body-frame omega
Rmid   = Rm(:, :, iw(round(end/2)));

s.V_lab_um_s       = V_lab;
s.V_body_um_s      = V_body;
s.U_fwd_um_s       = V_body(1);                  % mean forward speed along body +X
s.U_net_um_s       = norm(V_lab);                % net lab-frame chord speed
s.speed_um_s       = norm(V_body);
s.Omega_body_rad_s = Om_b;
s.omega_lab_rad_s  = Rmid * Om_b;
s.spin_Hz          = Om_b(1) / (2 * pi);
s.total_spin_turns = dTh_body_acc(1) / (2 * pi);

%% Screw motion (steady helical path) from mean body-frame (V_body, Om_b)
nOm = norm(Om_b);
s.helix = struct('omega_rad_s', nOm, 'pitch_um', NaN, 'radius_um', NaN, ...
    'axis_tilt_deg', NaN, 'period_s', NaN, 'screw_sign', 0, 'spin_vs_fwd', 0);
if nOm > 1e-6
    e_ax  = Om_b / nOm;
    V_par = dot(V_body, e_ax);
    s.helix.pitch_um      = 2 * pi * V_par / nOm;
    s.helix.radius_um     = norm(V_body - V_par * e_ax) / nOm;
    s.helix.axis_tilt_deg = acosd(min(1, abs(e_ax(1))));
    s.helix.period_s      = 2 * pi / nOm;
    s.helix.screw_sign    = sign(V_par);
    s.helix.spin_vs_fwd   = sign(s.spin_Hz) * sign(s.U_fwd_um_s + (s.U_fwd_um_s == 0));
end

% Oscillation of the body long-axis direction over the window
e_all = squeeze(Rm(:, 1, iw));
e_bar = mean(e_all, 2);
e_bar = e_bar / max(norm(e_bar), 1e-12);
dev   = acosd(max(-1, min(1, e_bar.' * e_all)));
s.axis_wobble_deg = max(dev);

%% Geometric clearance post-checks over the entire trajectory
nf = numel(results.X_pili_um);
f_body_min = inf;
for k = 1:nf
    Xp = results.X_pili_um{k};                   % [(N+1) x 3 x NT]
    for it = 1:NT
        Xb = (Xp(:, :, it) - Xc(:, it).') * Rm(:, :, it);
        f = sqrt((Xb(:, 1) / c_um).^2 + (Xb(:, 2) / a_um).^2 + (Xb(:, 3) / a_um).^2);
        f_body_min = min(f_body_min, min(f(2:end)));
    end
end
s.f_body_min = f_body_min;

d_min = inf;
if nf >= 2
    skip = 3;
    for it = 1:NT
        P1 = results.X_pili_um{1}(skip+1:end, :, it);
        P2 = results.X_pili_um{2}(skip+1:end, :, it);
        D = sqrt((P1(:,1) - P2(:,1).').^2 + (P1(:,2) - P2(:,2).').^2 + ...
                 (P1(:,3) - P2(:,3).').^2);
        d_min = min(d_min, min(D(:)));
    end
end
s.d_flag_min_um = d_min;
s.n_beats       = n_beats;

%% Tip 3-D orbit and instantaneous centerline non-planarity over the window
tip_oop       = NaN(1, nf);
tip_np_ratio  = NaN(1, nf);
tip_plane_rms = NaN(1, nf);
tip_z_rms     = NaN(1, nf);
tip_z_pp      = NaN(1, nf);
tip_orbit     = cell(1, nf);
c_np_all      = NaN(1, nf);

for k = 1:nf
    tip = squeeze(results.X_pili_um{k}(end, :, iw));   % [3 x n_w]
    for j = 1:size(tip, 2)
        tip(:, j) = Rm(:, :, iw(j)).' * (tip(:, j) - Xc(:, iw(j)));
    end
    tip_orbit{k} = tip;
    ry = sqrt(mean(tip(2, :).^2));
    rz = sqrt(mean(tip(3, :).^2));
    tip_oop(k)   = rz / max(ry, 1e-12);
    tip_z_rms(k) = rz;
    tip_z_pp(k)  = max(tip(3, :)) - min(tip(3, :));
    tc = tip - mean(tip, 2);
    sv = svd(tc, 'econ');
    if numel(sv) >= 3 && sv(1) > 1e-12
        tip_np_ratio(k)  = sv(3) / sv(1);
        tip_plane_rms(k) = sv(3) / sqrt(size(tc, 2));
    else
        tip_np_ratio(k)  = 0;
        tip_plane_rms(k) = 0;
    end

    % Instantaneous centerline 3-D non-planarity (mean sv3/sv1 over window)
    Xp = results.X_pili_um{k};
    sv3_t = zeros(1, numel(iw));
    for j = 1:numel(iw)
        P  = Xp(:, :, iw(j));
        Pc = P - mean(P, 1);
        svc = svd(Pc, 'econ');
        sv3_t(j) = svc(3) / max(svc(1), 1e-12);
    end
    c_np_all(k) = mean(sv3_t);
end

s.tip_oop_ratio             = tip_oop;
s.tip_oop_mean              = mean(tip_oop);
s.tip_nonplanar_ratio       = tip_np_ratio;
s.tip_nonplanar_mean        = mean(tip_np_ratio);
s.tip_plane_rms_um          = tip_plane_rms;
s.tip_plane_rms_mean        = mean(tip_plane_rms);
s.tip_z_rms_um              = tip_z_rms;
s.tip_z_pp_um               = tip_z_pp;
s.tip_orbit                 = tip_orbit;
s.centerline_nonplanar_mean = mean(c_np_all);

fprintf(['[SWIM] window %.0f-%.0f ms (%d beats): U_fwd = %.1f um/s, |V_body| = %.1f um/s, ', ...
    'V_body = [%.1f %.1f %.1f] um/s\n'], 1e3*t(i0), 1e3*t(i1), n_beats, ...
    s.U_fwd_um_s, s.speed_um_s, V_body);
fprintf(['[SWIM] spin = %+.3f Hz (%.2f turns), long-axis wobble = %.2f deg | ', ...
    'helix: period %.2f s, pitch %.1f um, radius %.2f um, tilt %.1f deg\n'], ...
    s.spin_Hz, s.total_spin_turns, s.axis_wobble_deg, s.helix.period_s, ...
    s.helix.pitch_um, s.helix.radius_um, s.helix.axis_tilt_deg);
fprintf('[SWIM] min body clearance f = %.3f (>1 required), min flag-flag dist = %.2f um\n', ...
    s.f_body_min, s.d_flag_min_um);
fprintf(['[SWIM] tip oop = %.3f | tip non-planar sv3/sv1 = %.4f (%.3f um) | ', ...
    'centerline non-planar = %.4f\n'], ...
    s.tip_oop_mean, s.tip_nonplanar_mean, s.tip_plane_rms_mean, s.centerline_nonplanar_mean);
end


function w = rotmat_log(R)
%ROTMAT_LOG  SO(3) logarithm: rotation vector w [3x1] such that exp([w]_x) = R
tr = max(-1, min(1, 0.5 * (trace(R) - 1)));
ang = acos(tr);
v = [R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];
if ang < 1e-6
    w = 0.5 * v;
else
    w = (0.5 * ang / sin(ang)) * v;
end
end


function R = rotmat_exp(w)
%ROTMAT_EXP  Rodrigues SO(3) exponential map from rotation vector w [3x1]
ang = norm(w);
if ang < 1e-12
    R = eye(3);
    return
end
u = w / ang;
K = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
R = eye(3) + sin(ang) * K + (1 - cos(ang)) * (K * K);
end


function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
