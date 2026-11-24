function [M_joint, info] = compute_effective_coupled_brownian(Z, params)
%COMPUTE_EFFECTIVE_COUPLED_BROWNIAN  Joint mobility of the body and the thermal pilus modes
%
%  Description:
%    Extends compute_effective_body_brownian to the pilus thermal-fluctuation
%    scheme (v9.0.0).  The thermal coordinates are
%
%        y = [ X_cm (3); body rotation (3); a_1 .. a_n ] ,
%
%    with a_j the amplitudes of the slow pilus bending modes selected by
%    init_pili_brownian_modes (params.brown_pili).  The joint mobility
%
%        y_dot = M_joint * [ f_body (6); g_1 .. g_n ]
%
%    is obtained by probing the fully coupled linear system with one unit
%    generalised force per thermal coordinate (one factorisation): the body
%    rows of the response are read directly and the modal rates through the
%    curvature kinematics (pili_bending_rates).
%
%    Every pilus mode that is NOT thermalised - all twist modes and all
%    faster bending modes of the rest-state mode analysis - is held fixed in
%    every probe (nx - n additional right-hand sides, params.brown_pili
%    .Q_hold / .A_hold).  Those modes relax within the deterministic
%    relaxation times (tens of microseconds for twist, at most a few
%    milliseconds for the third bending order), so over the correlation
%    time of a thermal mode the filament moves within the thermal-mode
%    subspace and the fast modes act as constraints: the mobility that
%    governs the thermal coordinates is the one with the fast modes held,
%    i.e. the inverse of the resistance restricted to the thermal
%    coordinates.  Leaving them force-free instead returns the
%    instantaneous mobility, which on a bent filament is inflated - by one
%    to two orders of magnitude for twist (a material-frame moment spins
%    the segments of a bent filament about their own axes), and by up to
%    several times for the first bending order (the rest-state mode
%    patterns mix with fast bending modes) - by motion that the elastic
%    moments stop long before a Brownian interval ends.  Noise sampled
%    from the inflated mobility is too weak for the response the
%    integration actually produces, and the modes fall short of
%    equipartition.  For a straight filament the hold is inactive on the
%    modal block (the modes are eigenvectors there).
%
%    The body block is therefore the mobility of the body with the
%    non-thermal pilus modes held and the thermal modes free.  It lies
%    between the fully force-free mobility of compute_effective_body_
%    brownian (valid below the fast relaxation times) and the rigid-shape
%    mobility (reached beyond the slowest thermal relaxation time).
%
%    The probe also yields the load that realises the thermal forces in the
%    time integration: the hydrodynamic force density carried by the pili
%    in the constrained probe motion (info.P_load, one column per thermal
%    coordinate).  With that load on the centrelines - plus the body
%    force / torque itself for the six body coordinates - the coupled
%    system responds with exactly the probe motion, so the thermal forces
%    do not kick the deterministic modes and the response of the
%    integration matches M_joint (fluctuation-dissipation consistency).  A
%    force on the centreline does not twist a naturally straight filament.
%
%    A load held in the body frame does not exert a constant modal force
%    while the filament moves.  The routine therefore also returns the
%    resulting drift coefficient (info.load_drift), which the sampler
%    compensates so that the thermal modes relax to equipartition.
%
%    The routine is deterministic and draws no random numbers.  The matrix
%    is symmetrised and, if round-off requires it, projected to SPD in the
%    single coherent dimensionless scale of the coupled system (forces in
%    F_0, torques and modal forces in F_0 L_0, rates in 1/T_0).
%
%  Syntax:
%    [M_joint, info] = compute_effective_coupled_brownian(Z, params)
%
%  Inputs:
%    Z      - [(7 + 3*sum(N_segments)) x 1] Global state vector
%    params - Parameter structure with params.brown_pili set by
%             init_pili_brownian_modes
%
%  Outputs:
%    M_joint - [(6+n) x (6+n)] Symmetric positive-definite dimensionless
%              joint mobility
%    info    - Struct: .rel_asym (before symmetrisation), .was_projected,
%              .min_eig_before, .rcond,
%              .P_load [3*N_nodes x (6+n)] centreline force density at the
%              pilus nodes per unit thermal force (body frame, node-major;
%              see pili_thermal_rhs).  Columns 1:6 accompany a unit body
%              force / torque, columns 7:end carry a unit modal force,
%              .pattern_resid mismatch between the modal forces exerted by
%              the loads and the target [0, I] (thermal-mode projection of
%              the generalised force; the axial torque density of the
%              probe motion is not part of the load),
%              .hold_pattern_resid relative mismatch of the loads in the
%              non-thermal (Q_hold) subspace,
%              .twist_row_resid relative mismatch of the loads on the twist
%              moment rows (comp == 3),
%              .hold_force Frobenius norm of the holding forces of the
%              modal probes (zero for straight pili),
%              .load_drift [n x 1] c_j = sum_k d g_j(P_k) / d a_k: per unit
%              kT, the spurious modal force left by holding the loads in
%              the body frame over an interval (at the straight rest
%              state zero for a pilus on the body axis, non-zero for the
%              symmetric modes of off-axis pili);
%              sample_brownian_coupled_force subtracts kT * c
%
%  Warnings:
%    compute_effective_coupled_brownian:AsymmetryHigh  (rel. asymmetry > 1e-2)
%    compute_effective_coupled_brownian:NotPosDef      (round-off projection)
%    compute_effective_coupled_brownian:IllConditioned (rcond < 1e-10)
%  Errors:
%    compute_effective_coupled_brownian:NotInitialised, :Indefinite,
%    :NonFiniteHold, :NonFiniteDrift
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11
%  Updated: 2026-11-24 — diagonal preconditioning & solve check for T_hh,
%                        symmetrise/SPD-repair M_joint before load_drift W,
%                        adaptive angle-capped h_fd, hold/twist row residuals.

if ~isfield(params, 'brown_pili') || ~isstruct(params.brown_pili) || ...
        ~isfield(params.brown_pili, 'enabled') || ~params.brown_pili.enabled
    error('compute_effective_coupled_brownian:NotInitialised', ...
        'params.brown_pili is not initialised (call init_pili_brownian_modes).');
end
bp = params.brown_pili;
if ~isfield(bp, 'Q_hold') || ~isfield(bp, 'A_hold')
    error('compute_effective_coupled_brownian:NotInitialised', ...
        ['params.brown_pili has no hold patterns (.Q_hold / .A_hold); ', ...
         'rebuild it with init_pili_brownian_modes.']);
end
n  = bp.n_modes;

logf_h = @(varargin) fprintf(varargin{:});
if isfield(params, 'logf') && ~isempty(params.logf)
    logf_h = params.logf;
end

%% Probe the coupled system: unit body force/torque and unit modal force
%  lin_sys = -Xi_gen, and an applied generalised force enters R with a minus
%  sign (R(1:3) = -Eg*g_hat for gravity, R(rows) = -Q_rows*g for a mode).
%  The last n_h columns are unit forces on the non-thermal modes; they are
%  used to hold those modes in the 6 + n thermal probes.
[lin_sys, geom, hydro] = build_coupled_lin_sys(Z, params, 'rhs_light');
n_state = size(lin_sys, 1);
layout  = bp.layout;
n_h     = size(bp.Q_hold, 2);

R_probe = zeros(n_state, 6 + n + n_h);
R_probe(1:6, 1:6)                 = -eye(6);
R_probe(bp.rows, 6 + (1:n))       = -bp.Q_rows;
R_probe(bp.rows, 6 + n + (1:n_h)) = -bp.Q_hold;

z_all = lin_sys \ R_probe;
check_linear_solve(lin_sys, z_all, R_probe, 'joint body-pilus mobility probe solve');
x_dot_all = pili_bending_rates(Z, params, layout, z_all);

%% Hold the non-thermal modes
%  The deterministic pilus modes (twist and the faster bending orders) are
%  stopped by their own elasticity long before a Brownian interval ends, so
%  the thermal probes are taken with those modes fixed: the holding forces
%  lam cancel the rates of the non-thermal amplitudes, A_hold' * x_dot = 0,
%  in each of the 6 + n probes.  At the straight rest state the hold does
%  not act on the modal block; on a bent filament it removes the fast
%  motion (segment spin about the local axis, admixture of fast bending
%  modes) that would otherwise be returned as "mobility" of a thermal mode.
if n_h > 0
    T_all = bp.A_hold' * x_dot_all;
    T_hh  = T_all(:, 6 + n + (1:n_h));
    T_hy  = T_all(:, 1:6 + n);
    % Symmetric diagonal preconditioning: T_hh spans O(N_k^4) across fast
    % bending and twist modes; scaling by |diag(T_hh)|^(-1/2) reduces the
    % condition number to ~1 at rest and ~O(10^2) on bent filaments.
    d_scale  = 1 ./ sqrt(max(abs(diag(T_hh)), realmin('double')));
    T_hh_pre = (d_scale .* T_hh) .* d_scale';
    T_hy_pre = d_scale .* T_hy;
    lam_pre  = T_hh_pre \ T_hy_pre;
    check_linear_solve(T_hh_pre, lam_pre, T_hy_pre, 'pilus non-thermal hold solve');
    lam = d_scale .* lam_pre;
    if any(~isfinite(lam(:)))
        error('compute_effective_coupled_brownian:NonFiniteHold', ...
            'Holding forces for non-thermal pilus modes are non-finite.');
    end
    z_probe = z_all(:, 1:6 + n)     - z_all(:, 6 + n + (1:n_h))     * lam;
    x_dot   = x_dot_all(:, 1:6 + n) - x_dot_all(:, 6 + n + (1:n_h)) * lam;
    hold_force = norm(lam(:, 7:end), 'fro');
else
    lam     = zeros(0, 6 + n);
    z_probe = z_all;
    x_dot   = x_dot_all;
    hold_force = 0;
end

a_dot   = bp.M_meas' * x_dot;
M_joint = [z_probe(1:6, :); a_dot];

%% Symmetrise (fluctuation-dissipation path only) and check BEFORE load_drift
%  Symmetrisation and SPD repair must precede W = z_probe / M_joint so that
%  the kinematic directions used for held-load drift match the exact
%  symmetrised mobility used by sample_brownian_coupled_force.
rel_asym = norm(M_joint - M_joint', 'fro') / max(norm(M_joint, 'fro'), eps);
if rel_asym > 1e-2
    logf_h(['[BROWN] WARNING: joint body-pilus mobility relative asymmetry = %.2e ', ...
            '(> 1e-2); symmetrisation is still applied. The deterministic/', ...
            'stochastic operator-splitting bias may be non-negligible.\n'], rel_asym);
    warning('compute_effective_coupled_brownian:AsymmetryHigh', ...
        'Joint body-pilus mobility relative asymmetry = %.2e (> 1e-2).', rel_asym);
end
M_joint = 0.5 * (M_joint + M_joint');

[M_joint, was_projected, min_eig_before] = project_spd(M_joint);
if was_projected
    logf_h(['[BROWN] WARNING: joint body-pilus mobility required an SPD projection ', ...
            '(min eigenvalue = %.2e).\n'], min_eig_before);
    warning('compute_effective_coupled_brownian:NotPosDef', ...
        ['Joint body-pilus mobility was not strictly positive definite ', ...
         '(min eigenvalue = %.2e); a round-off-scale SPD projection was applied.'], ...
        min_eig_before);
end

rc_mobility = rcond(M_joint);
if ~isfinite(rc_mobility) || rc_mobility < 1e-10
    logf_h('[BROWN] WARNING: joint body-pilus mobility rcond = %.2e (< 1e-10).\n', rc_mobility);
    warning('compute_effective_coupled_brownian:IllConditioned', ...
        'Joint body-pilus mobility rcond = %.2e (< 1e-10).', rc_mobility);
end

%% Centreline force and axial torque patterns that realise the thermal forces
%  In probe column c the pili carry the hydrodynamic force density
%  f = AQ_trans * z_probe(:, c) and the scalar axial torque density
%  m = AQ_twist_scalar * z_probe(:, c) along each segment tangent d3.
%  Applied together as an internal load, B_global * [f; m .* d3] reproduces
%  on the pilus rows the EXACT generalised force of the constrained probe
%  (to machine precision): the modal moment pattern (columns 7:end) together
%  with the holding forces of every non-thermal mode.  Including the probe's
%  balancing axial torque cancels viscous spin drag on bent filaments
%  without exciting twist (reducing twist energy by ~10^2 and eliminating
%  the ~3 % bent-state instantaneous response mismatch).  Force patterns are
%  stored in the body frame; the axial torque densities are frame-invariant
%  material scalars applied along the current segment tangents d3.
n_segs_all = sum(params.N_segments(:));
n_nodes    = n_segs_all + numel(params.N_segments);
f_rows     = 3 * params.N_body + (1 : 3 * n_nodes);
f_lab      = hydro.AQ_trans(f_rows, :) * z_probe;                % [3*n_nodes x (6+n)]
P_load     = reshape(geom.R_body' * reshape(f_lab, 3, []), 3 * n_nodes, 6 + n);

has_tau = isfield(hydro, 'AQ_twist_scalar') && ~isempty(hydro.AQ_twist_scalar);
if has_tau
    P_tau = hydro.AQ_twist_scalar * z_probe;                     % [n_segs_all x (6+n)]
else
    P_tau = zeros(0, 6 + n);
end

proj   = (bp.M_meas ./ layout.wK')';             % row-space load -> modal forces
Bf_all = hydro.B_global(bp.rows, f_rows) * f_lab;
if has_tau
    t_cols = 3 * (params.N_body + n_nodes) + (1 : 3 * n_segs_all);
    d3_all = zeros(n_segs_all, 3);
    seg_off = 0;
    for kp = 1:params.M_pili
        Nk = params.N_segments(kp);
        d3_all(seg_off + (1:Nk), :) = geom.d3{kp};
        seg_off = seg_off + Nk;
    end
    B_tau_scalar = hydro.B_global(bp.rows, t_cols(1:3:end)) .* d3_all(:, 1)' ...
                 + hydro.B_global(bp.rows, t_cols(2:3:end)) .* d3_all(:, 2)' ...
                 + hydro.B_global(bp.rows, t_cols(3:3:end)) .* d3_all(:, 3)';
    Bf_all = Bf_all + B_tau_scalar * P_tau;
end
g_check       = proj * Bf_all;
pattern_resid = norm(g_check - [zeros(n, 6), eye(n)], 'fro') / sqrt(n);

%  Diagnostics across the full row space (thermal + non-thermal / twist rows):
if n_h > 0
    proj_hold = (bp.A_hold ./ layout.wK')';
    hold_pattern_resid = norm(proj_hold * Bf_all + lam, 'fro') / max(norm(lam, 'fro'), 1);
else
    hold_pattern_resid = 0;
end
rhs_full_rows   = -lin_sys(bp.rows, :) * z_probe;
diff_rows       = Bf_all - rhs_full_rows;
is_tw_row       = (layout.comp == 3);
twist_row_resid = norm(diff_rows(is_tw_row, :), 'fro') / max(norm(rhs_full_rows, 'fro'), eps);

%% Drift of the held loads
%  The loads are held over the interval (centreline force in the body frame,
%  axial torque density as a material scalar along current d3), but the
%  modal force they exert changes as the filament moves (the arms and
%  directors of the moment balance change).  With g_j(P_k; y) the force on
%  mode j of the load pattern of mode k at configuration y, averaging the
%  second-order displacement of an interval over the thermal forces leaves,
%  next to the Ito drift kT div(M_joint), the term kT * M_joint * c with
%
%      c_j = sum_k  d g_j(P_k; y) / d a_k ,
%
%  i.e. a spurious generalised force kT * c.  The sampler subtracts it.
%  The derivative is taken along the kinematic direction of each thermal
%  mode, W = z_probe / M_joint (unit rate of a_k, other thermal coordinates
%  and all non-thermal modes fixed, body at rest), by central differences
%  via the matrix-free O(N_k) helper apply_pili_B_force_rows.
load_drift = zeros(n, 1);
if n > 0
    W = z_probe / M_joint;
    if any(~isfinite(W(:)))
        error('compute_effective_coupled_brownian:NonFiniteDrift', ...
            'Kinematic mode-direction matrix W = z_probe / M_joint is non-finite.');
    end
    h_fd_nom = 2e-3;                             % nominal step in mode amplitude (dimensionless)
    max_ang_step = 1e-2;                         % safety cap on largest Euler-angle step [rad]
    Z_c  = Z(:);
    for k = 1:n
        w_ang = W(7:end, 6 + k);
        s_ang = max(abs(w_ang));
        if s_ang > 0
            h_fd = min(h_fd_nom, max_ang_step / s_ang);
        else
            h_fd = h_fd_nom;
        end
        dz_k   = h_fd * w_ang;
        fk_lab = f_lab(:, 6 + k);
        if has_tau
            mk_tau = P_tau(:, 6 + k);
        else
            mk_tau = [];
        end
        g_pm = zeros(n, 2);
        for i_s = 1:2
            Z_d = Z_c;
            Z_d(8:end) = Z_d(8:end) + (3 - 2 * i_s) * dz_k;
            % Matrix-free O(N_k) evaluation of B_d(bp.rows, :) * [fk_lab; mk_tau .* d3_d]
            Bf_rows = apply_pili_B_force_rows(Z_d, fk_lab, params, layout, mk_tau);
            g_pm(:, i_s) = proj * Bf_rows;
        end
        load_drift = load_drift + (g_pm(:, 1) - g_pm(:, 2)) / (2 * h_fd);
    end
end

info = struct('rel_asym', rel_asym, 'was_projected', was_projected, ...
              'min_eig_before', min_eig_before, 'rcond', rc_mobility, ...
              'P_load', P_load, 'P_tau', P_tau, 'pattern_resid', pattern_resid, ...
              'hold_pattern_resid', hold_pattern_resid, ...
              'twist_row_resid', twist_row_resid, ...
              'hold_force', hold_force, 'load_drift', load_drift);

end


%% ==================== Local Helper Functions ====================

function Bf_rows = apply_pili_B_force_rows(Z, f_lab_vec, params, layout, m_twist_vec)
%APPLY_PILI_B_FORCE_ROWS  Matrix-free O(N_k) evaluation of B(bp.rows, :) * [f_lab; m_twist .* d3]
%  Matches build_force_balance exactly by using the scalar triple product
%  identity (dN x arg) . f = dN . (arg x f) and reverse cumulative sums.
    M      = params.M_pili;
    N_segs = params.N_segments(:)';
    [~, ~, pili_angles] = unpack_state(Z, N_segs);

    has_tau = (nargin >= 5) && ~isempty(m_twist_vec);
    Bf_rows = zeros(layout.nx, 1);
    out_pos = 0;
    f_pos   = 0;
    seg_pos = 0;

    for k = 1:M
        N_k     = N_segs(k);
        rot_k   = params.rot{k};
        ell_k   = params.ell(k);
        Nk2     = N_k^2;
        n_f_k   = 3 * (N_k + 1);

        if N_k < 2
            f_pos   = f_pos + n_f_k;
            seg_pos = seg_pos + N_k;
            continue;
        end

        % Reconstruct pilus centreline in rot_k basis (scaled by 1/ell_k, root = 0;
        % translation cancels in all arm differences inter(j) - X_k(i)/N_k).
        [d1_k, d2_k, d3_k] = pili_directors(pili_angles{k}.theta, ...
                                            pili_angles{k}.phi, ...
                                            pili_angles{k}.psi);
        [d1N, d2N, d3N] = compute_node_directors(d1_k, d2_k, d3_k, N_k);
        X_k = [0, 0, 0; cumsum(d3_k, 1) / N_k];          % [(N_k+1) x 3]

        % Rotate nodal force densities into rot_k basis and apply ell_k^3 factor
        f_k_lab = reshape(f_lab_vec(f_pos + (1:n_f_k)), 3, N_k + 1)';
        f_k     = (f_k_lab * rot_k') * (ell_k^3);        % [(N_k+1) x 3]

        % Accumulated hydrodynamic moment vector M_vec(i, :) for i = 1..N_k
        % Block 2: Diagonal (j = i)
        arg_diag = d3_k / (6 * Nk2);                     % [N_k x 3]
        M_vec    = cross(arg_diag, f_k(1:N_k, :), 2);

        % Block 3: Tip (j = N_k + 1)
        arg_tip = (X_k(N_k, :) - X_k(1:N_k, :)) / (2 * N_k) + d3_k(N_k, :) / (3 * Nk2);
        M_vec   = M_vec + cross(arg_tip, repmat(f_k(N_k + 1, :), N_k, 1), 2);

        % Block 1: Off-diagonal (j = i+1 .. N_k) via O(N_k) reverse cumsums
        %   sum_{j=i+1}^{N_k} (inter(j-1) - X_k(i)/N_k) x f_k(j)
        inter = (d3_k(2:end, :) + 2 * d3_k(1:end-1, :)) / (6 * Nk2) ...
              + (X_k(2:end-1, :) + X_k(1:end-2, :)) / (2 * N_k);   % [(N_k-1) x 3]
        f_off = f_k(2:N_k, :);                                     % [(N_k-1) x 3]
        cs_ixf = cumsum(cross(inter, f_off, 2), 1, 'reverse');     % [(N_k-1) x 3]
        cs_f   = cumsum(f_off, 1, 'reverse');                      % [(N_k-1) x 3]
        M_vec(1:N_k-1, :) = M_vec(1:N_k-1, :) + cs_ixf ...
                          - cross(X_k(1:N_k-1, :) / N_k, cs_f, 2);

        % Block 4: Axial torque density (upper triangular j >= i, ell_k^2 / N_k)
        if has_tau
            m_k   = m_twist_vec(seg_pos + (1:N_k));
            tau_k = d3_k .* m_k(:);                                % [N_k x 3] in rot_k basis
            M_vec = M_vec + (ell_k^2 / N_k) * cumsum(tau_k, 1, 'reverse');
        end

        % Project onto nodal directors at interior joints i = 2..N_k
        M_int = M_vec(2:N_k, :);
        proj_k = [sum(d1N(2:N_k, :) .* M_int, 2), ...
                  sum(d2N(2:N_k, :) .* M_int, 2), ...
                  sum(d3N(2:N_k, :) .* M_int, 2)]';                % [3 x (N_k-1)]
        n_out_k = 3 * (N_k - 1);
        Bf_rows(out_pos + (1:n_out_k)) = proj_k(:);

        out_pos = out_pos + n_out_k;
        f_pos   = f_pos   + n_f_k;
        seg_pos = seg_pos + N_k;
    end
end


function [A_spd, was_projected, min_eig_before] = project_spd(A)
%PROJECT_SPD  SPD repair in dimensionless coordinates (round-off scale only)
%  Same rule as compute_effective_body_brownian: a materially negative
%  eigenvalue signals a non-dissipative operator and is never hidden.
    A_sym = 0.5 * (A + A');
    [V, D] = eig(A_sym);
    lam = real(diag(D));
    min_eig_before = min(lam);
    spectral_scale = max(max(abs(lam)), realmin('double'));
    if min_eig_before < -1e-10 * spectral_scale
        error('compute_effective_coupled_brownian:Indefinite', ...
            ['Joint body-pilus mobility is materially indefinite: ', ...
             'lambda_min = %.3e, spectral scale = %.3e.'], ...
            min_eig_before, spectral_scale);
    end
    floor_eig = max(100 * size(A_sym, 1) * eps(spectral_scale), realmin('double'));
    if min_eig_before > floor_eig
        A_spd = A_sym;
        was_projected = false;
        return;
    end
    lam = max(lam, floor_eig);
    A_spd = V * diag(lam) * V';
    A_spd = 0.5 * (A_spd + A_spd');
    was_projected = true;
end
