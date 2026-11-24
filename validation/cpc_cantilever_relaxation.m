function summary = cpc_cantilever_relaxation(varargin)
%CPC_CANTILEVER_RELAXATION  First bending time of a clamped straight filament
%
%  Description:
%    External elastic check for the CPC program paper.  The body rates are
%    held at zero and the slowest bending eigenvalue of the coupled operator
%    is compared with the local cantilever formula
%
%        tau_ref = zeta_perp * L^4 / (EI * (beta_1 * L)^4),
%        beta_1 * L = 1.8751040684,
%
%    where zeta_perp = (F/U)/L is the perpendicular drag per unit length of
%    the same straight filament translated as a rigid rod by the project RSS
%    kernel.  The free-body bending time is reported beside the fixed-body
%    time and is not the pass line: that mode moves the body, so it is shorter.
%
%    Cases: L = 3 um at 10 and 20 segments per um, and L = 6 um at 10 per um.
%    Material is the project default (E = 3.279e8 Pa, r = 3 nm, mu = 1e-3).
%    The filament is attached at site 2, the off-axis site used by the
%    publication cantilever residual, so the root Euler angle stays away
%    from the coordinate singularity.
%
%    Pass: the 10 and 20 segment times at 3 um differ by less than 2 percent;
%    every tau/tau_ref lies in [0.90, 1.10]; the 6 um / 3 um time ratio lies
%    within 10 percent of the ratio of the two reference times.
%
%  Syntax:
%    summary = cpc_cantilever_relaxation()
%    summary = cpc_cantilever_relaxation('OutDir', outdir)
%
%  Output:
%    summary - struct with .table and .passed.  The table is also written to
%              <OutDir>/cantilever_relaxation.csv.
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
addpath(fullfile(root, 'MATLAB Add-On', ...
    'Suite of functions to perform uniform sampling of a sphere'));

if isempty(ip.Results.OutDir)
    outdir = fullfile(root, 'output', 'cpc_checks');
else
    outdir = ip.Results.OutDir;
end
if ~isfolder(outdir)
    mkdir(outdir);
end

spec = struct( ...
    'length_um', {3, 3, 6}, ...
    'seg_per_um', {10, 20, 10});
nolog = @(varargin) [];
beta4 = 1.8751040684039087^4;
rows = cell(1, numel(spec));
for k = 1:numel(spec)
    rows{k} = one_length(spec(k).length_um, spec(k).seg_per_um, beta4, nolog);
    fprintf(['L = %.0f um, %g /um: tau_fixed = %.4f s, tau_ref = %.4f s, ', ...
        'ratio = %.3f, tau_free = %.4f s\n'], ...
        rows{k}.length_um, rows{k}.seg_per_um, rows{k}.tau_fixed_s, ...
        rows{k}.tau_ref_s, rows{k}.ratio, rows{k}.tau_free_s);
end
T = vertcat(rows{:});

i10 = T.length_um == 3 & T.seg_per_um == 10;
i20 = T.length_um == 3 & T.seg_per_um == 20;
i6 = T.length_um == 6;
seg_rel = abs(T.tau_fixed_s(i10) - T.tau_fixed_s(i20)) / T.tau_fixed_s(i20);
scale = (T.tau_fixed_s(i6) / T.tau_fixed_s(i10)) / (T.tau_ref_s(i6) / T.tau_ref_s(i10));
passed = seg_rel < 0.02 && all(T.ratio >= 0.90) && all(T.ratio <= 1.10) ...
    && abs(scale - 1) <= 0.10;

summary = struct('table', T, 'passed', passed, ...
    'segment_relative_change', seg_rel, 'length_scale_over_reference', scale);
writetable(T, fullfile(outdir, 'cantilever_relaxation.csv'));
fprintf(['Cantilever check %s.  10 vs 20 seg: %.3f relative; ', ...
    'length ratio / reference ratio = %.3f\n'], pass_word(passed), seg_rel, scale);
fprintf('Table: %s\n', fullfile(outdir, 'cantilever_relaxation.csv'));
end


function word = pass_word(flag)
if flag
    word = 'PASSED';
else
    word = 'FAILED';
end
end


function row = one_length(length_um, seg_per_um, beta4, nolog)
E = 3.279e8;
r_pili = 3e-9;
mu = 1e-3;
EI = E * pi * r_pili^4 / 4;
N = max(8, round(seg_per_um * length_um));
pc = struct('type', 'linear', 'length_um', length_um, 'site', 2);
params = init_coupled_params('N_body', 120, 'N_segments', N, ...
    'Pili_Matrix', parse_pili_config(pc), 'E', E, 'r_pili', r_pili, ...
    'mu', mu, 'sigma', 0.2, 'hydro_radius_factor', sqrt(2), ...
    'minor_axis', 0.4e-6, 'major_axis', 1.25e-6, ...
    'rho_bac', 998.232, 'rho_fluid', 998.232, ...
    'logf', nolog, 'logf_file', nolog);
Z0 = init_coupled_state(params, nolog);
[Z0, params] = adaptive_basis_check(Z0, params);

tau_fixed = slowest_bending_time(Z0, params, true);
tau_free = slowest_bending_time(Z0, params, false);
zeta = perpendicular_drag(length_um * 1e-6, r_pili, sqrt(2), N, mu);
tau_ref = zeta * (length_um * 1e-6)^4 / (EI * beta4);
row = table(length_um, seg_per_um, N, EI, zeta, tau_fixed, tau_free, tau_ref, ...
    tau_fixed / tau_ref, ...
    'VariableNames', {'length_um', 'seg_per_um', 'N_seg', 'EI', 'zeta_perp', ...
    'tau_fixed_s', 'tau_free_s', 'tau_ref_s', 'ratio'});
end


function tau = slowest_bending_time(Z, params, body_fixed)
[lin_sys, ~] = build_coupled_lin_sys(Z, params);
[~, layout] = pili_bending_coords(Z, params);
nx = layout.nx;
if body_fixed
    nf = size(lin_sys, 1) - 6;
    block = lin_sys(7:end, 7:end);
    rhs = zeros(nf, nx);
    rows = layout.rows - 6;
    rhs(sub2ind([nf, nx], rows, 1:nx)) = -1 ./ layout.wE;
    z_f = block \ rhs;
    z_dot = [zeros(6, nx); z_f];
else
    rhs = zeros(size(lin_sys, 1), nx);
    rhs(sub2ind(size(rhs), layout.rows, 1:nx)) = -1 ./ layout.wE;
    z_dot = lin_sys \ rhs;
end
Mx = pili_bending_rates(Z, params, layout, z_dot);
Ms = 0.5 * (Mx + Mx');
K_half = sqrt((layout.wK .* layout.wE).');
S = 0.5 * ((K_half .* Ms) .* K_half.' + ((K_half .* Ms) .* K_half.').');
[U, D] = eig(S);
[lam, order] = sort(real(diag(D)), 'ascend');
U = U(:, order);
bend = sum(U(layout.comp < 3, :).^2, 1).' >= 0.5;
if lam(1) <= 0
    error('cpc_cantilever_relaxation:NotDissipative', ...
        'The curvature mobility has a non-positive eigenvalue (%.3e).', lam(1));
end
if ~any(bend)
    error('cpc_cantilever_relaxation:NoBendingMode', ...
        'No bending-dominated mode was found.');
end
tau = params.scales.T_0 / lam(find(bend, 1, 'first'));
end


function zeta = perpendicular_drag(L, r_pili, beta, n_seg, mu)
x = linspace(-L / 2, L / 2, n_seg + 1).';
X = [x, zeros(n_seg + 1, 2)];
U = [zeros(n_seg + 1, 1), ones(n_seg + 1, 1), zeros(n_seg + 1, 1)];
force_over_mu_U = rss_line_load(X, U, beta * r_pili);
zeta = mu * force_over_mu_U(2) / L;
end
