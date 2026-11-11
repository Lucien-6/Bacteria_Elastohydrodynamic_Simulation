function [fit, P] = chlamy_fit_moment(tab, opts)
%CHLAMY_FIT_MOMENT  Compact analytic fit of a tabulated active moment q(s,t)
%
%  Description:
%    Fits the measured-gait moment table of chlamy_inverse_moment (flagellum 1;
%    flagellum 2 is its mirror image, q2 = -q1) with
%
%        q(s,t) = c0(s) + A(s) cos( w t + phi(s) )
%                       + B(s) cos( 2 w t + 2 phi(s) + delta(s) ),   s = (j-1)/N, j = 2..N
%
%    where
%      A(s)     = A0 (1 - dA s) + Apk exp(-((s - sp)/wp)^2)    plateau with a
%                                                              proximal peak
%      phi(s)   = phi0 - k1 s - (k2 - k1) tau log(1 + exp((s - sb)/tau))
%                                                              two-slope phase
%                                                              (tau = 0.03)
%      B(s)     = B0 + Bpk exp(-((s - sp)/wp)^2)               second harmonic,
%                                                              same peak as A(s)
%      delta(s) = d0 + (d1 - d0) / (1 + exp(-(s - sd)/wd))     its phase offset
%                                                              relative to 2 phi
%      c0(s)    = piecewise-linear through the knots            static moment
%
%    A and phi are fitted to the complex first-harmonic coefficient of the
%    table (9 parameters: A0, dA, Apk, sp, wp, phi0, k1, k2, sb); B and delta
%    to its second harmonic (6 parameters: B0, Bpk, d0, d1, sd, wd, with sp,
%    wp, phi(s) fixed); c0 is a linear least-squares fit on 8 knots.  The
%    result is returned as a table in the format of chlamy_table_moment_fun
%    (nq = 2, or 1 with opts.harm2 = false) so that it can be run
%    (chlamy_run_case 'moment_table') without further changes.  Harmonics
%    above the second are dropped.
%
%  Syntax:
%    fit = chlamy_fit_moment(tab)
%    [fit, P] = chlamy_fit_moment(tab, opts)
%
%  Inputs:
%    tab  - struct from chlamy_inverse_moment
%    opts - (Optional) struct: .knots  c0 knots on s in [0,1]
%           (default [0 0.03 0.06 0.10 0.15 0.30 0.60 1.0]),
%           .harm2 (default true), .quiet
%
%  Outputs:
%    fit - moment table struct (same fields as TAB used by
%          chlamy_table_moment_fun, plus .P and fit errors .err_a1, .err_a2,
%          .err_c0, .err_q)
%    P   - struct of fitted parameters (units of EI/L): .A0 .dA .Apk .sp .wp
%          .phi0 .k1 .k2 .sb .tau; .B0 .Bpk .d0 .d1 .sd .wd; .knots .c0
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-10

if nargin < 2, opts = struct(); end
knots = getopt(opts, 'knots', [0 0.03 0.06 0.10 0.15 0.30 0.60 1.0]);
harm2 = getopt(opts, 'harm2', true);
quiet = getopt(opts, 'quiet', false);
tau = 0.03;

N  = tab.N_seg;
s  = (0:N-1).' / N;
c  = tab.coef{1};
a1 = c(:, 2) - 1i * c(:, 3);               % q1 = Re[a1 exp(i w t)]
a2 = c(:, 4) - 1i * c(:, 5);               % q2 = Re[a2 exp(2i w t)]
j  = (2:N).';                              % root node excluded (zeroed by the solver)
opt = optimset('MaxFunEvals', 4e4, 'MaxIter', 4e4, 'TolX', 1e-9, 'TolFun', 1e-12);

%% first harmonic: amplitude / phase
x0 = [1.8, 0.5, 11, 0.10, 0.06, 1.2, 6.0, 2.7, 0.2];
x = fminsearch(@(x) resid1(x, s(j), a1(j), tau), x0, opt);
x = fminsearch(@(x) resid1(x, s(j), a1(j), tau), x,  opt);
P = struct('A0', x(1), 'dA', x(2), 'Apk', x(3), 'sp', x(4), 'wp', x(5), 'phi0', x(6), ...
    'k1', x(7), 'k2', x(8), 'sb', x(9), 'tau', tau, 'knots', knots(:));

%% second harmonic
P.B0 = 0;  P.Bpk = 0;  P.d0 = 0;  P.d1 = 0;  P.sd = 0.25;  P.wd = 0.05;
if harm2
    y0 = [1.2, 2.0, 2.7, 5.7, 0.25, 0.05];
    cost2 = @(y) resid2(y, P, s(j), a2(j));
    y = fminsearch(cost2, y0, opt);
    y = fminsearch(cost2, y,  opt);
    P.B0 = y(1);  P.Bpk = y(2);  P.d0 = y(3);  P.d1 = y(4);  P.sd = y(5);  P.wd = abs(y(6));
end

%% static moment: linear least squares on hat functions
H = zeros(numel(j), numel(knots));
for i = 1:numel(knots)
    e = zeros(numel(knots), 1);  e(i) = 1;
    H(:, i) = interp1(knots(:), e, s(j), 'linear');
end
P.c0 = H \ c(j, 1);

%% evaluate on the nodes and assemble the table
nq = 1 + harm2;
q1c = eval_coef(P, s, harm2);
fit = struct('f_Hz', tab.f_Hz, 'T_0', tab.T_0, 'w_star', tab.w_star, 'nq', nq, 'N_seg', N);
fit.coef = {q1c, -q1c};
fit.Q_peak_Nm = zeros(1, 2);
EI_L = tab.Q_peak_Nm(1) / max(abs(tab.q_samples{1}(:)));     % EI/L_k [N m] (units of the table)
th = linspace(0, 2*pi, 97);
Bt = ones(numel(th), 1);
for n = 1:nq, Bt = [Bt, cos(n*th.'), sin(n*th.')]; end %#ok<AGROW>
for k = 1:2
    fit.Q_peak_Nm(k) = max(abs(fit.coef{k} * Bt.'), [], 'all') * EI_L;
end
fit.P = P;
a1m = q1c(:, 2) - 1i * q1c(:, 3);
fit.err_a1 = norm(a1m(j) - a1(j)) / norm(a1(j));
fit.err_a2 = NaN;
if harm2
    a2m = q1c(:, 4) - 1i * q1c(:, 5);
    fit.err_a2 = norm(a2m(j) - a2(j)) / norm(a2(j));
end
fit.err_c0 = norm(q1c(j, 1) - c(j, 1)) / norm(c(j, 1));
ph = tab.w_star * tab.t_s.' / tab.T_0;
Bq = ones(numel(ph), 1);
for n = 1:nq, Bq = [Bq, cos(n*ph), sin(n*ph)]; end %#ok<AGROW>
qf = q1c * Bq.';
fit.err_q = norm(qf(j, :) - tab.q_samples{1}(j, :), 'fro') / norm(tab.q_samples{1}(j, :), 'fro');
if ~quiet
    fprintf(['[FIT] A(s) = %.2f (1 - %.2f s) + %.2f exp(-((s - %.3f)/%.3f)^2) | ', ...
             'phi(s): phi0 = %.2f, slopes %.2f -> %.2f rad/L at s_b = %.3f\n'], ...
        P.A0, P.dA, P.Apk, P.sp, P.wp, P.phi0, P.k1, P.k2, P.sb);
    if harm2
        fprintf(['[FIT] B(s) = %.2f + %.2f exp(...) | delta(s): %.2f -> %.2f rad ', ...
                 '(logistic at s_d = %.3f, width %.3f)\n'], P.B0, P.Bpk, P.d0, P.d1, P.sd, P.wd);
    end
    fprintf('[FIT] c0 knots (s: value) %s\n', sprintf('%.2f:%.2f  ', [knots(:), P.c0].'));
    fprintf('[FIT] rel. error: a1 %.3f | a2 %.3f | c0 %.3f | q(s,t) %.3f\n', ...
        fit.err_a1, fit.err_a2, fit.err_c0, fit.err_q);
end
end

function coef = eval_coef(P, s, harm2)
A  = P.A0 * (1 - P.dA * s) + P.Apk * exp(-((s - P.sp) / P.wp).^2);
ph = P.phi0 - P.k1 * s - (P.k2 - P.k1) * P.tau * log1pexp((s - P.sb) / P.tau);
c0 = interp1(P.knots, P.c0, s, 'linear');
a1 = A .* exp(1i * ph);
coef = [c0, real(a1), -imag(a1)];         % q = c0 + Re[a1 e^{iwt}] = c0 + Re(a1) cos - Im(a1) sin
if harm2
    Bm = P.B0 + P.Bpk * exp(-((s - P.sp) / P.wp).^2);
    dl = P.d0 + (P.d1 - P.d0) ./ (1 + exp(-(s - P.sd) / P.wd));
    a2 = Bm .* exp(1i * (2 * ph + dl));
    coef = [coef, real(a2), -imag(a2)];
end
end

function r = resid1(x, s, a1, tau)
P = struct('A0', x(1), 'dA', x(2), 'Apk', x(3), 'sp', x(4), 'wp', x(5), 'phi0', x(6), ...
    'k1', x(7), 'k2', x(8), 'sb', x(9), 'tau', tau, 'knots', [0; 1], 'c0', [0; 0]);
cf = eval_coef(P, s, false);
r = sum(abs((cf(:, 2) - 1i * cf(:, 3)) - a1).^2);
end

function r = resid2(y, P, s, a2)
P.B0 = y(1);  P.Bpk = y(2);  P.d0 = y(3);  P.d1 = y(4);  P.sd = y(5);  P.wd = abs(y(6)) + 1e-6;
P.knots = [0; 1];  P.c0 = [0; 0];
cf = eval_coef(P, s, true);
r = sum(abs((cf(:, 4) - 1i * cf(:, 5)) - a2).^2);
end

function y = log1pexp(z)
y = max(z, 0) + log1p(exp(-abs(z)));
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
