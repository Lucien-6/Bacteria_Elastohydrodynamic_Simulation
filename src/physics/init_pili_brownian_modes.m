function bp = init_pili_brownian_modes(Z_ref, params, n_orders, opts)
%INIT_PILI_BROWNIAN_MODES  Slow bending modes of the pili that receive thermal noise
%
%  Description:
%    Builds the mode set used by the pilus thermal-fluctuation scheme
%    (project v9.0.0+).  Thermal forcing is applied to the cell body and to
%    the slowest elastohydrodynamic bending modes of the pili only; every
%    faster pilus mode stays deterministic.
%
%    The modes are computed once, at the reference state Z_ref, from the
%    discrete coupled operator itself:
%
%      1. Curvature coordinates x (pili_bending_coords): material curvature
%         and twist deviations at the interior joints.  The elastic energy is
%         E = (1/2) x' K x with the diagonal stiffness K = diag(wK .* wE).
%      2. Curvature mobility Mx: x_dot = Mx * f, where f is the generalised
%         force conjugate to x.  A unit f(r) is the moment-balance right-hand
%         side R(rows(r)) = -1 / wE(r); the response is read back through
%         the curvature kinematics (pili_bending_rates).  Body and all other
%         joints are force-free in each probe, so Mx carries the full
%         hydrodynamic coupling between filaments and with the moving body.
%      3. Relaxation modes: x_dot = -Mx K x.  With Ms = (Mx + Mx')/2 the
%         symmetric problem  sqrt(K) Ms sqrt(K) u = lambda u  gives the
%         relaxation rates lambda and K-orthonormal shapes w = K^(-1/2) u.
%      4. Mode selection:
%         - Default ('orders'): keep the 2 * n_orders * M_pili slowest
%           bending-dominated modes.  With opts.protect_pairs = true, if
%           the cut splits a near-degenerate cluster (gap_ratio < 1.15),
%           the cluster is expanded to the next clean spectral gap.
%         - Time-scale cutoff (opts.tau_cut_s > 0): keep every
%           bending-dominated mode with relaxation time tau_s >= tau_cut_s
%           (recommended for unequal-length pili; protects near-degenerate
%           pairs by default).
%
%    The remaining modes (all twist modes and the faster bending modes) are
%    stored as hold patterns: over the correlation time of a thermal mode
%    they are slaved by their own elasticity, so the joint mobility of the
%    thermal coordinates is probed with those modes held fixed.
%
%    Each kept mode j has an amplitude a_j = M_meas(:, j)' * x with elastic
%    energy a_j^2 / 2 (in F_0 L_0), and a generalised force g_j conjugate to
%    a_j whose moment pattern in the row space of the coupled system is
%    R(rows) = -Q_rows(:, j) * g_j.  The patterns are stored per joint in
%    the MATERIAL frame, so they rotate with the pili and are unaffected by
%    Euler-basis changes.  The moment pattern defines the modal force; in
%    the time integration it is exerted by the probe's centreline force
%    density and balancing axial torque density
%    (compute_effective_coupled_brownian, pili_thermal_rhs).
%
%  Syntax:
%    bp = init_pili_brownian_modes(Z_ref, params)
%    bp = init_pili_brownian_modes(Z_ref, params, n_orders)
%    bp = init_pili_brownian_modes(Z_ref, params, n_orders, opts)
%
%  Inputs:
%    Z_ref    - [(7 + 3*sum(N_segments)) x 1] Reference state (normally the
%               stress-free initial state Z0)
%    params   - Parameter structure (from init_coupled_params)
%    n_orders - Number of bending orders kept per pilus and polarisation
%               (default 2: first and second bending mode; ignored when
%               opts.tau_cut_s is set)
%    opts     - (Optional) Struct with fields:
%                 .tau_cut_s     positive scalar [s]; keep all bending
%                                modes with tau_s >= tau_cut_s (default [])
%                 .protect_pairs logical; expand cut across near-degenerate
%                                clusters with gap < 1.15 (default false
%                                for 'orders', true for 'tau_cut')
%
%  Outputs:
%    bp - Struct stored as params.brown_pili:
%      .enabled        true
%      .selection_mode 'orders' | 'tau_cut'
%      .tau_cut_s      cutoff time [s] (NaN for 'orders')
%      .protect_pairs  logical
%      .n_orders       bending orders requested (or max per-pilus order)
%      .n_modes        number of thermal modes kept
%      .layout         curvature-coordinate layout (pili_bending_coords)
%      .rows           [1 x nx] rows of R_global of the coordinates
%      .U_sel          [nx x n_modes] K-orthonormal eigenvectors kept
%      .Q_rows         [nx x n_modes] right-hand-side pattern per unit force
%      .M_meas         [nx x n_modes] amplitude read-out, a = M_meas' * x
%      .Q_hold         [nx x (nx - n_modes)] right-hand-side patterns of the
%                      non-thermal modes (twist and faster bending modes)
%      .A_hold         [nx x (nx - n_modes)] their amplitude read-out
%      .rate           [n_modes x 1] relaxation rate [1/T_0]
%      .tau_s          [n_modes x 1] relaxation time [s]
%      .order          [n_modes x 1] nominal bending order
%      .dom_pilus      [n_modes x 1] dominant pilus index (argmax pilus_share)
%      .pilus_order    [n_modes x 1] per-pilus bending order estimate
%      .bend_share     [n_modes x 1] bending fraction of the mode energy
%      .pilus_share    [n_modes x M_pili] energy share per pilus
%      .next_rate      slowest rate that is NOT thermalised [1/T_0]
%      .gap_ratio      next_rate / max(rate)
%      .asym_curv      relative asymmetry of Mx before symmetrisation
%      .joint_asym     relative asymmetry of the joint body-mode mobility
%      .pattern_resid, .hold_pattern_resid, .twist_row_resid
%
%  Errors:
%    init_pili_brownian_modes:NoPili, :BadOrders, :BadTauCut, :FreeTorsion,
%    :TooManyModes, :InsufficientBendingModes, :NotDissipative.
%  Warnings:
%    init_pili_brownian_modes:CurvedReference, :SmallGap, :PatternResidual.
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] Doi & Edwards (1986), The Theory of Polymer Dynamics, ch. 3.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11
%  Updated: 2026-11-24 — optional tau_cut_s time-scale selection and
%                        near-degenerate pair protection for unequal-length
%                        pili; dom_pilus / pilus_order metadata.

if nargin < 3 || isempty(n_orders)
    n_orders = 2;
end
if nargin < 4 || isempty(opts)
    opts = struct();
end
M = params.M_pili;
if M < 1
    error('init_pili_brownian_modes:NoPili', ...
        'Pilus thermal fluctuations need at least one pilus.');
end

use_tau_cut = isfield(opts, 'tau_cut_s') && ~isempty(opts.tau_cut_s);
tau_cut_s   = NaN;
if use_tau_cut
    tau_cut_s = double(opts.tau_cut_s);
    if ~(isscalar(tau_cut_s) && isnumeric(tau_cut_s) && isfinite(tau_cut_s) && tau_cut_s > 0)
        error('init_pili_brownian_modes:BadTauCut', ...
            'opts.tau_cut_s must be a positive finite scalar [s].');
    end
    selection_mode = 'tau_cut';
    default_protect = true;
else
    if ~(isscalar(n_orders) && isnumeric(n_orders) && isfinite(n_orders) && ...
            n_orders >= 1 && n_orders == round(n_orders))
        error('init_pili_brownian_modes:BadOrders', ...
            'The number of thermal bending orders must be a positive integer.');
    end
    selection_mode = 'orders';
    default_protect = false;
end
if isfield(opts, 'protect_pairs') && ~isempty(opts.protect_pairs)
    protect_pairs = logical(opts.protect_pairs);
else
    protect_pairs = default_protect;
end
n_modes_req = 2 * n_orders * M;

%  A free-torsion root adds the spin of the filament about its root axis, a
%  coordinate without elastic restoring moment that is neither thermalised
%  nor held by this scheme; on a bent filament the joint mobility then
%  loses its symmetry and definiteness.  The scheme requires clamped roots.
if isfield(params, 'root_coupling') && iscell(params.root_coupling) && ...
        any(strcmp(params.root_coupling, 'free_torsion'))
    error('init_pili_brownian_modes:FreeTorsion', ...
        ['Pilus thermal fluctuations require clamped roots; %d filament(s) use ', ...
         'root_coupling = free_torsion.'], ...
        sum(strcmp(params.root_coupling, 'free_torsion')));
end

logf_h = @(varargin) [];
if isfield(params, 'logf') && ~isempty(params.logf)
    logf_h = params.logf;
end

%% Reference state in a regular Euler basis
%  The initial state may sit on an Euler-angle singularity (run_simulation
%  repairs it before the first segment).  The modes are material-frame
%  quantities, so a local basis change does not alter them.
[Z_r, params_r] = adaptive_basis_check(Z_ref(:), params);

%  Warn when any pilus has non-zero intrinsic curvature (curved parametric
%  types 2-5, curved custom centrelines, or curvature overrides).  The
%  thermal centreline-load scheme is validated for naturally straight pili.
max_kb0 = 0;
for k_kb = 1:M
    if params_r.N_segments(k_kb) >= 2
        kb_int = params_r.kb0{k_kb}(2:params_r.N_segments(k_kb), :);
        max_kb0 = max(max_kb0, max(abs(kb_int(:))));
    end
end
if max_kb0 > 1e-6
    logf_h(['[BROWN] WARNING: non-zero reference curvature detected (max |kb0| = %.2e); ', ...
            'the pilus thermal scheme is validated for naturally straight filaments only.\n'], ...
            max_kb0);
    warning('init_pili_brownian_modes:CurvedReference', ...
        ['Non-zero reference curvature detected (max |kb0| = %.2e); ', ...
         'the pilus thermal scheme is validated for naturally straight filaments only.'], ...
        max_kb0);
end

[~, layout] = pili_bending_coords(Z_r, params_r);
nx      = layout.nx;
n_bend  = sum(layout.comp < 3);
if ~use_tau_cut && n_modes_req > n_bend
    error('init_pili_brownian_modes:TooManyModes', ...
        '%d thermal modes requested but the pili have only %d bending joints.', ...
        n_modes_req, n_bend);
end

%% Curvature mobility: unit generalised force on every joint coordinate
[lin_sys, ~] = build_coupled_lin_sys(Z_r, params_r);
R_probe = zeros(layout.n_state, nx);
R_probe(sub2ind(size(R_probe), layout.rows, 1:nx)) = -1 ./ layout.wE;
z_probe = lin_sys \ R_probe;
check_linear_solve(lin_sys, z_probe, R_probe, 'pilus curvature mobility probe solve');

Mx = pili_bending_rates(Z_r, params_r, layout, z_probe);          % x_dot = Mx * f
asym_curv = norm(Mx - Mx', 'fro') / max(norm(Mx, 'fro'), eps);
Ms = 0.5 * (Mx + Mx');

%% Relaxation modes of  x_dot = -Mx K x
K_diag = (layout.wK .* layout.wE)';                               % [nx x 1]
K_half = sqrt(K_diag);
S      = (K_half .* Ms) .* K_half';
S      = 0.5 * (S + S');
[U, D] = eig(S);
[lam, order_idx] = sort(real(diag(D)), 'ascend');
U = U(:, order_idx);

if lam(1) <= 0
    error('init_pili_brownian_modes:NotDissipative', ...
        ['The symmetrised curvature mobility has a non-positive eigenvalue (%.3e, ', ...
         'asym_curv = %.2e); either the hydrodynamic operator is not dissipative or ', ...
         'strongly curved 3D geometry breaks discrete curvature-moment conjugacy.'], ...
        lam(1), asym_curv);
end

bend_share_all = sum(U(layout.comp < 3, :).^2, 1)';
is_bend = bend_share_all >= 0.5;
bend_idx_all = find(is_bend);
tau_all_s = params.scales.T_0 ./ lam;

if use_tau_cut
    sel = bend_idx_all(tau_all_s(bend_idx_all) >= tau_cut_s);
    if isempty(sel)
        max_tau_bend = max(tau_all_s(bend_idx_all));
        error('init_pili_brownian_modes:InsufficientBendingModes', ...
            ['No bending-dominated mode has relaxation time >= tau_cut_s = %.4g s ', ...
             '(slowest bending mode is %.4g s). Reduce pili_tau_cut_s.'], ...
            tau_cut_s, max_tau_bend);
    end
else
    sel = bend_idx_all(1 : min(n_modes_req, numel(bend_idx_all)));
    if numel(sel) < n_modes_req
        error('init_pili_brownian_modes:InsufficientBendingModes', ...
            ['Requested %d bending-dominated thermal modes (bend_share >= 0.5), ', ...
             'but found only %d; strong bend-twist coupling on curved filaments ', ...
             'prevents clean separation of bending and twist modes.'], ...
            n_modes_req, numel(sel));
    end
end

% Optional near-degenerate pair/cluster protection: avoid splitting a
% polarisation pair or tightly clustered unequal-length multiplet (gap < 1.15).
n_expanded = 0;
if protect_pairs
    pair_gap_tol = 1.15;
    pos_end = numel(sel);
    while pos_end < numel(bend_idx_all)
        r_last = lam(bend_idx_all(pos_end));
        r_next = lam(bend_idx_all(pos_end + 1));
        if r_next / r_last < pair_gap_tol
            pos_end = pos_end + 1;
            n_expanded = n_expanded + 1;
        else
            break;
        end
    end
    sel = bend_idx_all(1:pos_end);
end

n_modes = numel(sel);
rest = bend_idx_all(~ismember(bend_idx_all, sel));

U_sel = U(:, sel);
%  Fix the sign of every mode (largest component positive) so the mode set
%  does not depend on the arbitrary sign returned by eig.
[~, i_max] = max(abs(U_sel), [], 1);
sgn = sign(U_sel(sub2ind(size(U_sel), i_max, 1:n_modes)));
sgn(sgn == 0) = 1;
U_sel = U_sel .* sgn;

rate = lam(sel);
if isempty(rest)
    next_rate = Inf;
else
    next_rate = lam(rest(1));
end

pilus_share = zeros(n_modes, M);
for k = 1:M
    pilus_share(:, k) = sum(U_sel(layout.pilus == k, :).^2, 1)';
end
[~, dom_pilus] = max(pilus_share, [], 2);
pilus_order = zeros(n_modes, 1);
for k = 1:M
    in_k = find(dom_pilus == k);
    if ~isempty(in_k)
        pilus_order(in_k) = ceil((1:numel(in_k))' / 2);
    end
end
if use_tau_cut || n_expanded > 0
    order_vec    = pilus_order;
    n_orders_eff = max(pilus_order);
else
    order_vec    = ceil((1:n_modes)' / (2 * M));
    n_orders_eff = n_orders;
end

%% Pack
bp = struct();
bp.enabled        = true;
bp.selection_mode = selection_mode;
bp.tau_cut_s      = tau_cut_s;
bp.protect_pairs  = protect_pairs;
bp.n_orders       = n_orders_eff;
bp.n_modes        = n_modes;
bp.layout         = layout;
bp.rows           = layout.rows;
bp.U_sel          = U_sel;
bp.Q_rows         = sqrt(layout.wK ./ layout.wE)' .* U_sel;
bp.M_meas         = K_half .* U_sel;
%  Non-thermal modes (all twist modes and the faster bending modes): the
%  K-orthogonal complement of the thermal modes.  They are held fixed when
%  the joint mobility is probed (compute_effective_coupled_brownian).
U_hold         = U(:, setdiff(1:nx, sel, 'stable'));
bp.Q_hold      = sqrt(layout.wK ./ layout.wE)' .* U_hold;
bp.A_hold      = K_half .* U_hold;
bp.rate        = rate;
bp.tau_s       = params.scales.T_0 ./ rate;
bp.order       = order_vec;
bp.dom_pilus   = dom_pilus;
bp.pilus_order = pilus_order;
bp.bend_share  = bend_share_all(sel);
bp.pilus_share = pilus_share;
bp.next_rate   = next_rate;
bp.gap_ratio   = next_rate / max(rate);
bp.asym_curv   = asym_curv;

%% Joint mobility and centreline loads at the reference state
params_r.brown_pili = bp;
[~, joint_info]       = compute_effective_coupled_brownian(Z_r, params_r);
bp.joint_asym         = joint_info.rel_asym;
bp.pattern_resid      = joint_info.pattern_resid;
bp.hold_pattern_resid = joint_info.hold_pattern_resid;
bp.twist_row_resid    = joint_info.twist_row_resid;

%% Report
if use_tau_cut
    logf_h('[BROWN] Pilus thermal modes: %d modes selected by tau_cut_s = %.4g s (%d pili, pair-protected = %d)\n', ...
        n_modes, tau_cut_s, M, protect_pairs);
else
    logf_h('[BROWN] Pilus thermal modes: %d modes (%d bending order(s) x 2 polarisations x %d pili)\n', ...
        n_modes, n_orders, M);
    if n_expanded > 0
        logf_h('[BROWN]   pair protection expanded cut by +%d mode(s) to avoid splitting a cluster\n', ...
            n_expanded);
    end
end
logf_h('[BROWN]   curvature-mobility asymmetry = %.2e (symmetrised for the mode analysis)\n', asym_curv);
unequal_lengths = (M > 1) && (max(params.ell) - min(params.ell) > 1e-6 * max(params.ell));
if unequal_lengths || use_tau_cut || n_expanded > 0
    for k = 1:M
        in_k = (dom_pilus == k);
        if any(in_k)
            logf_h('[BROWN]   pilus %d (ell=%.3f L0): %d dominant mode(s), tau %.4g .. %.4g s\n', ...
                k, params.ell(k), sum(in_k), min(bp.tau_s(in_k)), max(bp.tau_s(in_k)));
        else
            logf_h('[BROWN]   pilus %d (ell=%.3f L0): 0 dominant modes (all faster than cut)\n', ...
                k, params.ell(k));
        end
    end
else
    for o = 1:n_orders_eff
        in_o = bp.order == o;
        logf_h('[BROWN]   order %d: relaxation time %.4g .. %.4g s (rate %.4g .. %.4g 1/s)\n', o, ...
            min(bp.tau_s(in_o)), max(bp.tau_s(in_o)), ...
            min(1 ./ bp.tau_s(in_o)), max(1 ./ bp.tau_s(in_o)));
    end
end
logf_h('[BROWN]   fastest thermal mode: %.4g s | slowest deterministic bending mode: %.4g s (gap ratio %.2f)\n', ...
    min(bp.tau_s), params.scales.T_0 / next_rate, bp.gap_ratio);
if bp.gap_ratio < 1.5
    if bp.gap_ratio < 1.15
        gap_note = ' (likely splits a near-degenerate polarisation pair or unequal-length cluster)';
    else
        gap_note = '';
    end
    logf_h(['[BROWN] WARNING: the thermal/deterministic mode cut passes through a cluster of ', ...
            'relaxation rates (gap ratio %.2f < 1.5)%s; modes of similar speed are treated differently.\n'], ...
            bp.gap_ratio, gap_note);
    warning('init_pili_brownian_modes:SmallGap', ...
        ['Thermal/deterministic mode cut passes through a cluster of relaxation rates ', ...
         '(gap ratio %.2f < 1.5)%s.'], bp.gap_ratio, gap_note);
end
logf_h(['[BROWN]   joint body-mode mobility asymmetry = %.2e | centreline-load ', ...
        'pattern residual = %.2e (hold = %.2e, twist rows = %.2e)\n'], ...
    bp.joint_asym, bp.pattern_resid, bp.hold_pattern_resid, bp.twist_row_resid);
if bp.pattern_resid > 1e-2
    logf_h(['[BROWN] WARNING: the centreline loads reproduce the modal forces only to ', ...
            '%.1e; the thermal modes carry twist that a centreline force cannot drive.\n'], ...
            bp.pattern_resid);
    warning('init_pili_brownian_modes:PatternResidual', ...
        'Centreline loads reproduce the modal forces only to %.1e.', bp.pattern_resid);
end

end
