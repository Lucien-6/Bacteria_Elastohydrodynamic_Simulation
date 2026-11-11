function s = chlamy_swim_analysis(results, opts)
%CHLAMY_SWIM_ANALYSIS  Swimming kinematics and geometric post-checks of a Chlamydomonas run
%
%  Description:
%    Evaluates, over an integer number of beat periods after a start-up
%    transient,
%      * net displacement and mean swimming velocity of the cell centre,
%        decomposed in the (time-averaged) body frame; the body long axis
%        is body-X and the flagella sit at its +X end (forward direction),
%      * mean angular velocity of the body (axis-angle of the net rotation)
%        and its component along the long axis (spin),
%      * peak-to-peak oscillation of the body's long-axis direction,
%      * minimum clearance between every flagellar node and the spheroid
%        surface (f = sqrt((x/c)^2 + (y/a)^2 + (z/a)^2), body frame; f < 1 is
%        penetration) and the minimum flagellum-flagellum node distance
%        (nodes closer than 3 segments to the root are skipped).
%    The code has no contact model, so the clearances are the post-check
%    that the run is physically admissible.
%
%  Syntax:
%    s = chlamy_swim_analysis(results)
%    s = chlamy_swim_analysis(results, opts)
%
%  Inputs:
%    results - Struct from run_simulation
%    opts    - (Optional) struct:
%                .f_Hz    reference beat frequency [Hz]    (default 50)
%                .n_skip  beats discarded at the start     (default 1)
%                .c_um    spheroid semi-axis along X [um]  (default 5)
%                .a_um    spheroid semi-axes along Y,Z     (default 4)
%
%  Output:
%    s - Struct with fields V_lab_um_s [3x1], V_body_um_s [3x1] (components
%        along body X, Y, Z at the window mid-point), U_fwd_um_s,
%        speed_um_s, omega_lab_rad_s [3x1], spin_Hz (about the body long
%        axis), axis_wobble_deg, f_body_min, d_flag_min_um, n_beats
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-06

proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(proj_root, 'src')));
if nargin < 2, opts = struct(); end
f_Hz   = getopt(opts, 'f_Hz', 50);
n_skip = getopt(opts, 'n_skip', 1);
c_um   = getopt(opts, 'c_um', 5);
a_um   = getopt(opts, 'a_um', 4);
T_b = 1 / f_Hz;

t  = results.t(:);
NT = numel(t);
n_beats = floor((t(end) - n_skip * T_b) / T_b + 1e-9);
if n_beats < 1
    error('chlamy_swim_analysis:TooShort', 'Run is shorter than %d beat(s).', n_skip + 1);
end
t0 = n_skip * T_b;  t1 = t0 + n_beats * T_b;
i0 = find(t >= t0 - 1e-9, 1, 'first');
i1 = find(t <= t1 + 1e-9, 1, 'last');
iw = i0:i1;
Tw = t(i1) - t(i0);

Rm = zeros(3, 3, NT);
for it = 1:NT, Rm(:, :, it) = quat_to_rotmat(results.q(:, it)); end

%% Translation
Xc = results.X_cm_um;
V_lab = (Xc(:, i1) - Xc(:, i0)) / Tw;                 % [um/s]
ex = squeeze(Rm(:, 1, iw));  ex = mean(ex, 2);  ex = ex / norm(ex);
Rmid = Rm(:, :, iw(round(end/2)));
V_body = Rmid.' * V_lab;
s.V_lab_um_s  = V_lab;
s.V_body_um_s = V_body;
s.U_fwd_um_s  = dot(V_lab, ex);
s.speed_um_s  = norm(V_lab);

%% Rotation (net rotation over the window; exact when it is < pi)
dR = Rm(:, :, i1) * Rm(:, :, i0).';
ang = acos(max(-1, min(1, (trace(dR) - 1) / 2)));
if ang < 1e-12
    ax = [0; 0; 0];
else
    ax = [dR(3,2) - dR(2,3); dR(1,3) - dR(3,1); dR(2,1) - dR(1,2)] / (2 * sin(ang));
end
s.omega_lab_rad_s = ax * ang / Tw;
s.spin_Hz = dot(s.omega_lab_rad_s, ex) / (2 * pi);

% Screw motion (steady helical path) from the window-mean velocities in the body frame
Om_b = Rmid.' * s.omega_lab_rad_s;
nOm = norm(Om_b);
s.Omega_body_rad_s = Om_b;
s.helix = struct('omega_rad_s', nOm, 'pitch_um', NaN, 'radius_um', NaN, ...
    'axis_tilt_deg', NaN, 'period_s', NaN);
if nOm > 1e-6
    e_ax = Om_b / nOm;
    V_par = dot(V_body, e_ax);
    s.helix.pitch_um  = 2 * pi * V_par / nOm;
    s.helix.radius_um = norm(V_body - V_par * e_ax) / nOm;
    s.helix.axis_tilt_deg = acosd(min(1, abs(e_ax(1))));        % screw axis vs body long axis (body +X)
    s.helix.period_s  = 2 * pi / nOm;
end

% oscillation of the long-axis direction inside the window
e_all = squeeze(Rm(:, 1, iw));
e_bar = mean(e_all, 2);  e_bar = e_bar / norm(e_bar);
dev = acosd(max(-1, min(1, e_bar.' * e_all)));
s.axis_wobble_deg = max(dev);

%% Geometric post-checks over the whole run
nf = numel(results.X_pili_um);
f_body_min = inf;
for k = 1:nf
    Xp = results.X_pili_um{k};                        % [(N+1) x 3 x NT]
    for it = 1:NT
        Xb = (Xp(:, :, it) - Xc(:, it).') * Rm(:, :, it);   % body frame
        f = sqrt((Xb(:,1)/c_um).^2 + (Xb(:,2)/a_um).^2 + (Xb(:,3)/a_um).^2);
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
s.n_beats = n_beats;

fprintf(['[SWIM] window %.0f-%.0f ms (%d beats): U_fwd = %.1f um/s, |V| = %.1f um/s, ', ...
    'V_body = [%.1f %.1f %.1f] um/s\n'], 1e3*t(i0), 1e3*t(i1), n_beats, ...
    s.U_fwd_um_s, s.speed_um_s, V_body);
fprintf(['[SWIM] omega_lab = [%.2f %.2f %.2f] rad/s, spin about long axis = %.3f Hz, ', ...
    'long-axis wobble = %.1f deg\n'], s.omega_lab_rad_s, s.spin_Hz, s.axis_wobble_deg);
fprintf(['[SWIM] omega_body = [%.2f %.2f %.2f] rad/s | helix: period %.2f s, pitch %.1f um, radius %.1f um, ', ...
    'axis tilt %.1f deg from the body axis\n'], Om_b, s.helix.period_s, s.helix.pitch_um, ...
    s.helix.radius_um, s.helix.axis_tilt_deg);
fprintf('[SWIM] min flagellum/spheroid clearance f = %.3f (>1 required), min flagellum-flagellum distance = %.2f um\n', ...
    s.f_body_min, s.d_flag_min_um);
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
