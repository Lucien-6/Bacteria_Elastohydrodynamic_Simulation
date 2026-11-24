function tab = chlamy_make_perp_envelope(tab, opts)
%CHLAMY_MAKE_PERP_ENVELOPE  Build 3-D out-of-plane active-moment table (1w + 2w)
%
%  Description:
%    Adds tab.coef_perp{1:2} for a proximal out-of-plane bending envelope
%    phase-locked to the first-harmonic phase phi_par(s) of tab.coef{1}:
%
%        q_perp(s, t) = W(s) * [ A_perp * cos(w t + phi_par(s) + delta_perp)
%                              + B_perp * cos(2 w t + 2 phi_par(s) + delta_perp2) ],
%        W(s)         = exp(-s / s_w),   W(root) = 0,
%        A_perp       = r_perp  * rms(q_par^osc over active nodes s, phase t),
%        B_perp       = r_perp2 * rms(q_par^osc over active nodes s, phase t),
%
%    with zero static bias (c0_perp == 0).  Default parameters correspond to
%    the calibrated 3-D drive (r_perp=0.25, r_perp2=0.05, s_w=0.20,
%    delta_perp=delta_perp2=pi/2).
%
%    Chirality convention (single multiplication):
%      Both flagella receive the identical envelope in tab.coef_perp = {c1, c1}.
%      Cis/trans chirality is applied ONCE at run time via ov.perp_sign in
%      chlamy_run_case (default [+1, -1], yielding exact C2 chiral symmetry
%      with tilt_deg = [+17, +17]).
%
%  Syntax:
%    tab = chlamy_make_perp_envelope(tab)
%    tab = chlamy_make_perp_envelope(tab, opts)
%
%  Inputs:
%    tab  - Moment table struct from chlamy_fit_moment / chlamy_inverse_moment
%    opts - (Optional) struct:
%             .r_perp       1w amplitude ratio A_perp / rms(q_par^osc) (0.25)
%             .r_perp2      2w amplitude ratio B_perp / rms(q_par^osc) (0.05)
%             .s_w          proximal decay length in s/L               (0.20)
%             .delta_perp   1w phase offset vs phi_par [rad]           (pi/2)
%             .delta_perp2  2w phase offset vs 2*phi_par [rad]         (pi/2)
%             .quiet        suppress console summary                   (false)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-23
%  Updated: 2026-11-24  arguments block & calibrated 3-D defaults

arguments
    tab  (1,1) struct
    opts (1,1) struct = struct()
end

r_perp  = getopt(opts, 'r_perp', 0.25);
r_perp2 = getopt(opts, 'r_perp2', 0.05);
s_w     = getopt(opts, 's_w', 0.20);
delta   = getopt(opts, 'delta_perp', pi/2);
delta2  = getopt(opts, 'delta_perp2', pi/2);
quiet   = logical(getopt(opts, 'quiet', false));

assert(s_w > 0, 'chlamy_make_perp_envelope:DecayLength', 's_w must be positive.');
assert(isfield(tab, 'coef') && iscell(tab.coef) && ~isempty(tab.coef), ...
    'chlamy_make_perp_envelope:InvalidTable', 'tab.coef{1} is required.');

cq = tab.coef{1};
N  = size(cq, 1);
nq = tab.nq;
s  = (0:N-1).' / max(N, 1);
W  = exp(-s / s_w);
W(1) = 0;   % root node is clamped reaction

% Oscillating RMS of q_par over active balance nodes (2..N) and one beat
n_ph = 96;
th = linspace(0, 2*pi, n_ph + 1);
th(end) = [];
B = zeros(n_ph, 2*nq);
for n = 1:nq
    B(:, 2*n-1) = cos(n * th.');
    B(:, 2*n)   = sin(n * th.');
end
q_osc   = cq(2:end, 2:end) * B.';    % [(N-1) x n_ph]
rms_par = sqrt(mean(q_osc(:).^2));
A_perp  = r_perp  * rms_par;
B_perp  = r_perp2 * rms_par;

% First-harmonic complex amplitude a1 = c_cos - i c_sin => Re[a1 exp(i w t)]
a1  = cq(:, 2) - 1i * cq(:, 3);
phi = angle(a1);
psi = phi + delta;

c1 = zeros(N, 2*nq + 1);
c1(:, 2) =  W * A_perp .* cos(psi);  % coeff of cos(w t)
c1(:, 3) = -W * A_perp .* sin(psi);  % coeff of sin(w t)
if r_perp2 ~= 0
    assert(nq >= 2, 'chlamy_make_perp_envelope:Harm2', ...
        'r_perp2 ~= 0 requires tab.nq >= 2.');
    psi2 = 2 * phi + delta2;
    c1(:, 4) =  W * B_perp .* cos(psi2);  % coeff of cos(2 w t)
    c1(:, 5) = -W * B_perp .* sin(psi2);  % coeff of sin(2 w t)
end
c1(1, :) = 0;

% Identical envelope for both flagella; chirality applied ONCE by ov.perp_sign
tab.coef_perp = {c1, c1};
tab.perp = struct('r_perp', r_perp, 'r_perp2', r_perp2, 's_w', s_w, ...
    'delta_perp', delta, 'delta_perp2', delta2, ...
    'A_perp', A_perp, 'B_perp', B_perp, 'rms_par', rms_par);

if ~quiet
    fprintf(['[PERP] r_perp=%.3f  r_perp2=%.3f  s_w=%.3f  delta=%.2f rad  ', ...
        'A_perp=%.4f  B_perp=%.4f (rms_par=%.4f)\n'], ...
        r_perp, r_perp2, s_w, delta, A_perp, B_perp, rms_par);
end
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
