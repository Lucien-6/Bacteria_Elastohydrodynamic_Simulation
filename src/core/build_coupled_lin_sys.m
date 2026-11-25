function [lin_sys, geom, hydro] = build_coupled_lin_sys(Z, params, export_mode)
%BUILD_COUPLED_LIN_SYS  Assemble the coupled body-pili linear system
%
%  Description:
%    Builds the full coupled linear system for the body-pili
%    elastohydrodynamic problem:
%      lin_sys * z_dot = R
%    where lin_sys encodes the hydrodynamic resistance (mobility
%    inversion + force/moment balance) and root boundary conditions
%    (clamped or free-torsion independently for each pilus).
%
%    The input ODE state Z has length 7 + 3*sum(N_segments), whereas
%    z_dot, R, and both dimensions of lin_sys are
%    6 + 3*sum(N_segments).  The difference is the four-component
%    quaternion state versus the three-component angular velocity.
%    With N_total = size(geom.X_all, 1), the principal assembly dimensions are
%      n_state = 6 + 3*sum(N_segments),
%      n_vel   = 3*N_total + sum(N_segments), and
%      n_force = 3*N_total + 3*sum(N_segments).
%    Accordingly Q_global is [n_vel x n_state], while the hydrodynamic
%    force map uses n_force rows before reduction to lin_sys.
%
%    This function extracts Steps 1-7 of ode_coupled_rhs.m into a
%    reusable subroutine, shared by ode_coupled_rhs (deterministic
%    RHS) and compute_effective_body_brownian (Brownian diffusion).
%
%    The elastic restoring moments (RHS) and state-rate packing are
%    NOT included here — those remain in ode_coupled_rhs.m.
%
%  Syntax:
%    [lin_sys, geom] = build_coupled_lin_sys(Z, params)
%    [lin_sys, geom, hydro] = build_coupled_lin_sys(Z, params)
%
%  Inputs:
%    Z      - [(7 + 3*sum(N_segments)) x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    lin_sys - [(6 + 3*sum(N_segments)) x
%               (6 + 3*sum(N_segments))] Coupled linear system matrix
%              with each pilus's configured root BC applied.  A clamped
%              root replaces the theta_1, phi_1, and psi_1 rows; a
%              free-torsion root replaces only theta_1 and phi_1 and
%              retains the psi_1 natural moment-balance row.  Represents
%              -Xi_gen in dimensionless units.
%    geom    - Geometry struct from reconstruct_geometry (needed by
%              caller for elastic moments, director frames, etc.)
%    hydro   - (Optional, 3rd output) Hydrodynamic intermediate struct
%              for background flow computation.  Fields:
%                .B_global, .M_weighted, .dM_weighted,
%                .M_twist, .D_twist, .C_twist, .C_twist_reverse,
%                .W_trans_diag, .W_twist_diag, .Minv_C_twist,
%                .S_twist, .dS_twist, .M_rotlet, .n_trans, .AQ_trans
%              .dM_weighted is a reusable LU decomposition object for
%              M_weighted \ rhs solves (lets apply_Ainv_u skip a second
%              full LU factorisation in the background-flow branch of
%              ode_coupled_rhs).
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] Cortez et al. (2005), Phys. Fluids 17, 031504.
%    [3] Cortez (2018), J. Comput. Phys. 375, 783-796.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-31
%  Updated: 2026-07-13
%  Updated: 2026-11-04 — v8.7.2: comment and message-text corrections only; no behavioural change.
%  Updated: 2026-11-12 — v9.0.0: hydro.AQ_trans exported (third output only).
%  Updated: 2026-11-12 — v9.0.0: B_global assembly moved unchanged to build_force_balance (results bit-identical).

if nargin < 3 || isempty(export_mode)
    export_mode = 'full';
end

N_body = params.N_body;
M      = params.M_pili;
N_segs = params.N_segments;

%% Step 1: Reconstruct geometry
geom    = reconstruct_geometry(Z, params);
X_all   = geom.X_all;
X_cm    = geom.X_cm;
N_total = size(X_all, 1);

%% Fast Path: Body-Frame Block Schur Complement using precomputed M_BB Cholesky
%  Exact algebraic reformulation of Steps 2-6 when precomputed rigid-body
%  factors are available and a full dense M_weighted is not required.
use_block_schur = isfield(params, 'body_precomp') && isstruct(params.body_precomp) && ...
    ~isempty(params.body_precomp) && ...
    isfield(params.body_precomp, 'dM_BB') && ...
    size(params.body_precomp.dM_BB, 1) == 3 * N_body && ...
    params.body_precomp.N_body == N_body && ...
    abs(params.body_precomp.eps_body - params.eps_body) <= 1e-15 * max(params.eps_body, 1) && ...
    (nargout < 3 || strcmp(export_mode, 'rhs_light'));

if use_block_schur
    [lin_sys, hydro] = assemble_block_schur_lin_sys(geom, params, nargout >= 3);
    % Proceed directly to Step 7 (root boundary conditions)
    if M > 0
        [lin_sys, geom] = apply_root_boundary_conditions(lin_sys, geom, params);
    end
    if any(~isfinite(lin_sys(:)))
        error('numerical_diagnostics:NonFiniteMatrix', ...
            'Coupled linear system contains non-finite values.');
    end
    return;
end

%% Step 2: Build global mobility matrices (full dense path)
[M_weighted, M_twist, M_rotlet, D_twist] = build_global_resistance(X_all, params, geom);
if any(~isfinite(M_weighted(:)))
    error('numerical_diagnostics:NonFiniteMobility', ...
        'Translational mobility matrix contains non-finite values.');
end
if ~isempty(M_twist) && any(~isfinite(M_twist(:)))
    error('numerical_diagnostics:NonFiniteMobility', ...
        'Twist mobility matrix contains non-finite values.');
end

%% Step 3-4: Build B, Q matrices
%  NOTE: Step 5 (elastic restoring moments RHS) is intentionally omitted
%  from this function — it remains in ode_coupled_rhs.m where it is
%  assembled after the linear system solve.  Hence the step numbering
%  jumps from 4 to 6 here.
n_state = 6 + 3 * sum(N_segs);
n_vel   = 3 * N_total + sum(N_segs);

%% Build Q_global: [n_vel x n_state]
Q_global = zeros(n_vel, n_state);

%  Vectorised body block of Q_global.  For i = 1..N_body:
%    rows (3i-2, 3i-1, 3i), cols 1..3 are identity 3x3
%    rows (3i-2, 3i-1, 3i), cols 4..6 are -skew(r_i)
%  where r_i = X_all(i,:) - X_cm.  -skew(r) = [0 z -y; -z 0 x; y -x 0].
r_body = X_all(1:N_body, :) - X_cm';              % [N_body x 3]
rbx    = r_body(:,1);  rby = r_body(:,2);  rbz = r_body(:,3);

idx_bx = (1 : 3 : 3*N_body)';
idx_by = (2 : 3 : 3*N_body)';
idx_bz = (3 : 3 : 3*N_body)';

Q_global(idx_bx, 1) = 1;  Q_global(idx_by, 2) = 1;  Q_global(idx_bz, 3) = 1;
Q_global(idx_bx, 5) =  rbz;   Q_global(idx_bx, 6) = -rby;
Q_global(idx_by, 4) = -rbz;   Q_global(idx_by, 6) =  rbx;
Q_global(idx_bz, 4) =  rby;   Q_global(idx_bz, 5) = -rbx;

pt_offset  = N_body;
col_offset = 6;
twist_row  = 3 * N_total;

for k = 1:M
    N_k   = N_segs(k);
    rot_k = params.rot{k};
    ell_k = params.ell(k);
    ds_k  = ell_k / N_k;

    theta_k = geom.pili_angles{k}.theta;
    phi_k   = geom.pili_angles{k}.phi;

    r_root = X_all(pt_offset + 1, :)' - X_cm;

    %  Body rigid-body velocity columns (1:6): identical for every node of pilus k
    J_base = [eye(3), -skew_matrix(r_root)];
    trans_rows = 3*pt_offset + (1 : 3*(N_k+1));
    Q_global(trans_rows, 1:6) = repmat(J_base, N_k + 1, 1);

    %  Pili angle-rate columns without dense kron(eye, rot_k') multiplication:
    %  rotate the 3 x (2*N_k) tangent derivatives directly by rot_k', then
    %  broadcast with the strict lower-triangular cumulative-sum mask.
    cT = cos(theta_k(:))';     sT = sin(theta_k(:))';     % 1 x N_k
    cP = cos(phi_k(:))';       sP = sin(phi_k(:))';
    d_lab = rot_k' * [cT .* cP, -sT .* sP; ...
                      cT .* sP,  sT .* cP; ...
                      -sT,       zeros(1, N_k)] * ds_k;   % [3 x 2*N_k]

    maskL   = tril(true(N_k + 1, N_k), -1);               % (N_k+1) x N_k
    cols_th = col_offset + (1 : N_k);
    cols_ph = col_offset + N_k + (1 : N_k);
    cols_ps = col_offset + 2*N_k + (1 : N_k);

    rx_k = 3*pt_offset + (1 : 3 : 3*(N_k+1));
    ry_k = rx_k + 1;
    rz_k = rx_k + 2;

    Q_global(rx_k, cols_th) = maskL .* d_lab(1, 1:N_k);
    Q_global(rx_k, cols_ph) = maskL .* d_lab(1, N_k+1:end);
    Q_global(ry_k, cols_th) = maskL .* d_lab(2, 1:N_k);
    Q_global(ry_k, cols_ph) = maskL .* d_lab(2, N_k+1:end);
    Q_global(rz_k, cols_th) = maskL .* d_lab(3, 1:N_k);
    Q_global(rz_k, cols_ph) = maskL .* d_lab(3, N_k+1:end);

    %  Twist rows: omega_axial = cos(theta_k) * dphi_k + dpsi_k
    tw_r = twist_row + (1 : N_k);
    Q_global(sub2ind([n_vel, n_state], tw_r, cols_ph)) = cT;
    Q_global(sub2ind([n_vel, n_state], tw_r, cols_ps)) = 1;

    pt_offset  = pt_offset + N_k + 1;
    col_offset = col_offset + 3 * N_k;
    twist_row  = twist_row + N_k;
end

%% Build B_global: [n_state x n_force]
%  Force / moment-balance operator.  It depends on the configuration only
%  and is assembled by build_force_balance (the assembly formerly inlined
%  here, moved unchanged), which the pilus thermal scheme also evaluates
%  on its own (v9.0.0, held-load drift).
[B_global, twist_seg_lengths, B_blocks] = build_force_balance(geom, params);

%% Step 6: Assemble lin_sys = -B * AQ
n_trans  = 3 * N_total;
Q_trans  = Q_global(1:n_trans, :);

% Determine solve-check cadence for the large multi-RHS mobility solve.
% Always run the full backward residual check when hydro is exported
% (health checks, Brownian probes, flow-field reconstruction) or when
% params.solve_check_mode == 'strict'.  On the hot ODE RHS path, periodic
% sampling (default stride 50) verifies the full residual regularly while
% checking finiteness on every call.
do_full_mob_check = should_run_full_solve_check(params, nargout >= 3);

% ---- Direction 1: factor M_weighted ONCE, reuse across solves.            ----
dM_weighted = decomposition(M_weighted, 'lu', 'CheckCondition', false);

C_twist      = [];
C_twist_reverse = [];
Minv_C_twist = [];
S_twist      = [];
dS_twist     = [];
W_trans_diag = repelem(params.quad_weights(:), 3);
W_twist_diag = twist_seg_lengths(:);

if ~isempty(M_twist)
    Q_twist = Q_global(n_trans+1:end, :);

    % Weighted-reciprocal mobility inverse:
    %   [u_t] = [Mtt       C_forward] [f]
    %   [u_w]   [C_reverse Mll      ] [m],
    % where C_reverse = W_m \ (C_forward' * W_f).
    [C_twist, C_twist_reverse, W_trans_diag, W_twist_diag] = ...
        assemble_twist_coupling(M_rotlet, D_twist, ...
                                params.quad_weights, twist_seg_lengths);
    % Single combined multi-RHS solve against the shared LU factors
    RHS_comb = [Q_trans, C_twist];
    Sol_comb = dM_weighted \ RHS_comb;
    if do_full_mob_check
        check_linear_solve(M_weighted, Sol_comb, RHS_comb, ...
            'coupled translational mobility solve');
    elseif any(~isfinite(Sol_comb(:)))
        error('numerical_diagnostics:NonFiniteValue', ...
            'Non-finite value detected in coupled translational mobility solve.');
    end
    Minv_Q_trans = Sol_comb(:, 1:n_state);
    Minv_C_twist = Sol_comb(:, n_state+1:end);

    S_twist      = M_twist - C_twist_reverse * Minv_C_twist;
    dS_twist     = decomposition(S_twist, 'lu', 'CheckCondition', false);

    rhs_schur       = Q_twist - C_twist_reverse * Minv_Q_trans;
    AQ_twist_scalar = dS_twist \ rhs_schur;
    if do_full_mob_check
        check_linear_solve(S_twist, AQ_twist_scalar, rhs_schur, ...
            'coupled twist Schur solve');
    elseif any(~isfinite(AQ_twist_scalar(:)))
        error('numerical_diagnostics:NonFiniteValue', ...
            'Non-finite value detected in coupled twist Schur solve.');
    end
    AQ_trans        = Minv_Q_trans - Minv_C_twist * AQ_twist_scalar;

    %  Block-sparse assembly of lin_sys = -B_global * [AQ_trans; D_twist * AQ_twist_scalar]
    %  B_global(1:6, :) spans all forces/couples; rows 7:end are strictly
    %  block-diagonal per pilus with zero coupling to the body columns.
    lin_sys = zeros(n_state, n_state);
    lin_sys(1:6, :) = -B_global(1:6, 1:n_trans) * AQ_trans;

    %  Add net axial couple to body torque rows 4:6 without forming dense AQ_twist:
    %    sum_s L_s * d3_s * AQ_twist_scalar(s, :)
    N_seg_all = sum(N_segs);
    d3_all = zeros(3, N_seg_all);
    seg_off = 0;
    for k = 1:M
        N_k = N_segs(k);
        d3_all(:, seg_off + (1:N_k)) = geom.d3{k}';
        seg_off = seg_off + N_k;
    end
    lin_sys(4:6, :) = lin_sys(4:6, :) - (d3_all .* twist_seg_lengths) * AQ_twist_scalar;

    %  Per-pilus moment-balance rows (block-sparse; contracts B_tau * D_twist_k directly)
    seg_off = 0;
    for k = 1:M
        N_k = N_segs(k);
        bk  = B_blocks(k);
        d3_k = geom.d3{k};                           % [N_k x 3]
        B_tau_scalar = bk.B_tau(:, 1:3:end) .* d3_k(:, 1)' ...
                     + bk.B_tau(:, 2:3:end) .* d3_k(:, 2)' ...
                     + bk.B_tau(:, 3:3:end) .* d3_k(:, 3)';   % [3*N_k x N_k]
        lin_sys(bk.row_start:bk.row_end, :) = ...
            -(bk.B_force * AQ_trans(bk.f_start:bk.f_end, :) ...
            + B_tau_scalar * AQ_twist_scalar(seg_off + (1:N_k), :));
        seg_off = seg_off + N_k;
    end
else
    AQ_trans = dM_weighted \ Q_trans;
    if do_full_mob_check
        check_linear_solve(M_weighted, AQ_trans, Q_trans, ...
            'coupled translational mobility solve');
    elseif any(~isfinite(AQ_trans(:)))
        error('numerical_diagnostics:NonFiniteValue', ...
            'Non-finite value detected in coupled translational mobility solve.');
    end
    lin_sys = -B_global * AQ_trans;
end

%% Step 6b: Optionally export hydrodynamic intermediates (for background flow)
if nargout >= 3
    hydro.B_global      = B_global;
    hydro.M_weighted    = M_weighted;
    hydro.dM_weighted   = dM_weighted;
    hydro.M_twist       = M_twist;
    hydro.D_twist       = D_twist;
    hydro.C_twist       = C_twist;
    hydro.C_twist_reverse = C_twist_reverse;
    hydro.W_trans_diag  = W_trans_diag;
    hydro.W_twist_diag  = W_twist_diag;
    hydro.Minv_C_twist  = Minv_C_twist;
    hydro.S_twist       = S_twist;
    hydro.dS_twist      = dS_twist;
    hydro.M_rotlet      = M_rotlet;
    hydro.n_trans        = n_trans;
    % Force map of the generalised rates (v9.0.0): f_trans = AQ_trans * z_dot
    % is the hydrodynamic force (body points) / force density (pilus nodes)
    % exerted on the fluid; m_twist = AQ_twist_scalar * z_dot is the scalar
    % axial torque density on the pilus segments.  Both are used to build the
    % constrained thermal load patterns.
    hydro.AQ_trans        = AQ_trans;
    if ~isempty(M_twist)
        hydro.AQ_twist_scalar = AQ_twist_scalar;
    else
        hydro.AQ_twist_scalar = zeros(0, n_state);
    end
end

%% Step 7: Apply root boundary conditions (clamped or free-torsion)
if M > 0
    [lin_sys, geom] = apply_root_boundary_conditions(lin_sys, geom, params);
end

if any(~isfinite(lin_sys(:)))
    error('numerical_diagnostics:NonFiniteMatrix', ...
        'Coupled linear system contains non-finite values.');
end

end


%% ==================== Local Helper Functions ====================

function [lin_sys, hydro] = assemble_block_schur_lin_sys(geom, params, want_hydro)
%ASSEMBLE_BLOCK_SCHUR_LIN_SYS  Exact body-frame Block Schur linear system
    pre    = params.body_precomp;
    N_body = params.N_body;
    nB     = pre.nB;
    M      = params.M_pili;
    N_segs = params.N_segments;
    R      = geom.R_body;
    Rt     = R';
    X_cm   = geom.X_cm;
    N_tot  = size(geom.X_all, 1);
    n_trans = 3 * N_tot;
    n_state = 6 + 3 * sum(N_segs);

    hydro = struct();

    % ---- Bare body (M_pili == 0): rotate constant precomputed 6x6 resistance ----
    if M == 0
        T_rot = blkdiag(R, R);
        lin_sys = T_rot * pre.lin_sys_bare_loc * T_rot';
        if want_hydro
            [B_global, ~] = build_force_balance(geom, params);
            AQ_bf = pre.U_QB_loc;
            AQ_lab = reshape(R * reshape(AQ_bf, 3, []), nB, 6);
            AQ_lab(:, 1:3) = AQ_lab(:, 1:3) * Rt;
            AQ_lab(:, 4:6) = AQ_lab(:, 4:6) * Rt;
            hydro.B_global    = B_global;
            hydro.M_twist     = [];
            hydro.n_trans     = n_trans;
            hydro.AQ_trans    = AQ_lab;
            hydro.AQ_twist_scalar = zeros(0, 6);
            hydro.block_schur = struct('R_body', R, 'dM_BB', pre.dM_BB, ...
                'nB', nB, 'nP', 0, 'U_BP', [], 'M_PB', [], 'dS_PP', []);
        end
        return;
    end

    % ---- Body + Pili: transform pili points & tangents to body frame ----
    N_pili = N_tot - N_body;
    nP     = 3 * N_pili;
    Ns     = sum(N_segs);

    X_all_bf = zeros(N_tot, 3);
    X_all_bf(1:N_body, :)   = params.body_local';
    X_all_bf(N_body+1:end, :) = (geom.X_all(N_body+1:end, :) - X_cm') * R;

    d3_all_lab = zeros(Ns, 3);
    seg_off = 0;
    for k = 1:M
        N_k = N_segs(k);
        d3_all_lab(seg_off + (1:N_k), :) = geom.d3{k};
        seg_off = seg_off + N_k;
    end
    d3_all_bf = d3_all_lab * R;

    % Build only non-BB mobility blocks in body frame (skip M_BB, M_rotlet, D_twist)
    [M_PB, M_BP, M_PP, C_twist_bf, M_twist, L_all] = ...
        build_pili_mobility_blocks(X_all_bf, d3_all_bf, params);
    if any(~isfinite(M_PB(:))) || any(~isfinite(M_BP(:))) || ...
            any(~isfinite(M_PP(:))) || any(~isfinite(M_twist(:)))
        error('numerical_diagnostics:NonFiniteMobility', ...
            'Hydrodynamic mobility blocks contain non-finite values.');
    end

    % Reciprocal reverse twist coupling (frame-invariant weights)
    W_trans_diag = repelem(params.quad_weights(:), 3);
    W_twist_diag = L_all(:);
    C_rev_bf = (C_twist_bf' .* W_trans_diag') ./ W_twist_diag;

    % Assemble Q_P [nP x n_state] and Q_twist [Ns x n_state] in body frame
    Q_P  = zeros(nP, n_state);
    Q_tw = zeros(Ns, n_state);
    po_p = 0;
    co   = 6;
    tr   = 0;
    for k = 1:M
        N_k    = N_segs(k);
        rk_bf  = params.rot{k} * R;       % maps body-frame vectors -> pilus basis
        rk_bfT = rk_bf';
        ds_k   = params.ell(k) / N_k;
        th_k   = geom.pili_angles{k}.theta;
        ph_k   = geom.pili_angles{k}.phi;

        r_root_bf = X_all_bf(N_body + po_p + 1, :)';
        J_base    = [eye(3), -skew_matrix(r_root_bf)];
        Q_P(3*po_p + (1:3*(N_k+1)), 1:6) = repmat(J_base, N_k + 1, 1);

        cT = cos(th_k(:))';  sT = sin(th_k(:))';
        cP = cos(ph_k(:))';  sP = sin(ph_k(:))';
        d_bf = rk_bfT * [cT .* cP, -sT .* sP; ...
                         cT .* sP,  sT .* cP; ...
                         -sT,       zeros(1, N_k)] * ds_k;

        maskL   = tril(true(N_k + 1, N_k), -1);
        cols_th = co + (1:N_k);
        cols_ph = co + N_k + (1:N_k);
        cols_ps = co + 2*N_k + (1:N_k);
        rx_k    = 3*po_p + (1 : 3 : 3*(N_k+1));
        ry_k    = rx_k + 1;
        rz_k    = rx_k + 2;

        Q_P(rx_k, cols_th) = maskL .* d_bf(1, 1:N_k);
        Q_P(rx_k, cols_ph) = maskL .* d_bf(1, N_k+1:end);
        Q_P(ry_k, cols_th) = maskL .* d_bf(2, 1:N_k);
        Q_P(ry_k, cols_ph) = maskL .* d_bf(2, N_k+1:end);
        Q_P(rz_k, cols_th) = maskL .* d_bf(3, 1:N_k);
        Q_P(rz_k, cols_ph) = maskL .* d_bf(3, N_k+1:end);

        tw_r = tr + (1:N_k);
        Q_tw(sub2ind([Ns, n_state], tw_r, cols_ph)) = cT;
        Q_tw(sub2ind([Ns, n_state], tw_r, cols_ps)) = 1;

        po_p = po_p + N_k + 1;
        co   = co   + 3 * N_k;
        tr   = tr   + N_k;
    end

    % ---- Block Schur Complement eliminating all 3*N_body body unknowns ----
    %  Uses precomputed Cholesky factor dM_BB and U_QB_loc = M_BB \ Q_B6_loc.
    C_B = C_twist_bf(1:nB, :);
    C_P = C_twist_bf(nB+1:end, :);

    X_B  = pre.dM_BB \ [M_BP, C_B];
    U_BP = X_B(:, 1:nP);
    U_CB = X_B(:, nP+1:end);

    S_PP  = M_PP - M_PB * U_BP;
    dS_PP = decomposition(S_PP, 'lu', 'CheckCondition', false);

    RHS_P_Q = Q_P;
    RHS_P_Q(:, 1:6) = RHS_P_Q(:, 1:6) - M_PB * pre.U_QB_loc;
    RHS_P_C = C_P - M_PB * U_CB;

    RHS_P_comb = [RHS_P_Q, RHS_P_C];
    Sol_P      = dS_PP \ RHS_P_comb;
    do_full_chk = should_run_full_solve_check(params, want_hydro);
    if do_full_chk
        check_linear_solve(S_PP, Sol_P, RHS_P_comb, ...
            'coupled block-Schur pili mobility solve');
    elseif any(~isfinite(Sol_P(:)))
        error('numerical_diagnostics:NonFiniteValue', ...
            'Non-finite value detected in coupled block-Schur pili mobility solve.');
    end

    fP_Q = Sol_P(:, 1:n_state);
    fP_C = Sol_P(:, n_state+1:end);

    fB_Q = -U_BP * fP_Q;
    fB_Q(:, 1:6) = fB_Q(:, 1:6) + pre.U_QB_loc;
    fB_C = U_CB - U_BP * fP_C;

    Minv_Q_bf = [fB_Q; fP_Q];
    Minv_C_bf = [fB_C; fP_C];

    % ---- Twist Schur Complement (size Ns x Ns) ----
    S_twist  = M_twist - C_rev_bf * Minv_C_bf;
    dS_twist = decomposition(S_twist, 'lu', 'CheckCondition', false);

    rhs_tw = Q_tw - C_rev_bf * Minv_Q_bf;
    AQ_twist_scalar_bf = dS_twist \ rhs_tw;
    if do_full_chk
        check_linear_solve(S_twist, AQ_twist_scalar_bf, rhs_tw, ...
            'coupled twist Schur solve');
    elseif any(~isfinite(AQ_twist_scalar_bf(:)))
        error('numerical_diagnostics:NonFiniteValue', ...
            'Non-finite value detected in coupled twist Schur solve.');
    end

    AQ_trans_bf = Minv_Q_bf - Minv_C_bf * AQ_twist_scalar_bf;

    % ---- Map hydrodynamic force operator back to lab frame ----
    %  Force vectors transform by R on rows; body-velocity columns 1:3 and 4:6
    %  transform by R' on the right (V_bf = R' * V_lab, Omega_bf = R' * Omega_lab).
    AQ_trans = reshape(R * reshape(AQ_trans_bf, 3, []), n_trans, n_state);
    AQ_trans(:, 1:3) = AQ_trans(:, 1:3) * Rt;
    AQ_trans(:, 4:6) = AQ_trans(:, 4:6) * Rt;

    AQ_twist_scalar = AQ_twist_scalar_bf;
    AQ_twist_scalar(:, 1:3) = AQ_twist_scalar(:, 1:3) * Rt;
    AQ_twist_scalar(:, 4:6) = AQ_twist_scalar(:, 4:6) * Rt;

    % ---- Block-sparse assembly of lin_sys = -B_global * AQ ----
    [B_global, twist_seg_lengths, B_blocks] = build_force_balance(geom, params);

    lin_sys = zeros(n_state, n_state);
    lin_sys(1:6, :) = -B_global(1:6, 1:n_trans) * AQ_trans;
    lin_sys(4:6, :) = lin_sys(4:6, :) - (d3_all_lab' .* twist_seg_lengths) * AQ_twist_scalar;

    seg_off = 0;
    for k = 1:M
        N_k = N_segs(k);
        bk  = B_blocks(k);
        d3_k = geom.d3{k};
        B_tau_scalar = bk.B_tau(:, 1:3:end) .* d3_k(:, 1)' ...
                     + bk.B_tau(:, 2:3:end) .* d3_k(:, 2)' ...
                     + bk.B_tau(:, 3:3:end) .* d3_k(:, 3)';
        lin_sys(bk.row_start:bk.row_end, :) = ...
            -(bk.B_force * AQ_trans(bk.f_start:bk.f_end, :) ...
            + B_tau_scalar * AQ_twist_scalar(seg_off + (1:N_k), :));
        seg_off = seg_off + N_k;
    end

    if want_hydro
        hydro.B_global        = B_global;
        hydro.M_twist         = M_twist;
        hydro.C_twist_reverse = C_rev_bf;
        hydro.Minv_C_twist    = Minv_C_bf;
        hydro.S_twist         = S_twist;
        hydro.dS_twist        = dS_twist;
        hydro.d3_all_lab      = d3_all_lab;
        hydro.n_trans         = n_trans;
        hydro.AQ_trans        = AQ_trans;
        hydro.AQ_twist_scalar = AQ_twist_scalar;
        hydro.block_schur     = struct('R_body', R, 'dM_BB', pre.dM_BB, ...
            'nB', nB, 'nP', nP, 'U_BP', U_BP, 'M_PB', M_PB, 'dS_PP', dS_PP);
    end
end


function [lin_sys, geom] = apply_root_boundary_conditions(lin_sys, geom, params)
%APPLY_ROOT_BOUNDARY_CONDITIONS  Enforce clamped or free-torsion root BCs (Step 7)
    M      = params.M_pili;
    N_segs = params.N_segments;
    col_offset = 6;
    row_offset = 6;

    for k = 1:M
        N_k    = N_segs(k);
        rot_k  = params.rot{k};
        theta1 = geom.pili_angles{k}.theta(1);
        phi1   = geom.pili_angles{k}.phi(1);

        sT = sin(theta1);  cT = cos(theta1);
        sP = sin(phi1);    cP = cos(phi1);

        if abs(sT) < 1e-12
            logf_h = get_params_logf(params);
            logf_h('[ODE] WARNING: sin(theta_1) = %.2e for pili %d -- Euler angle singularity.\n', sT, k);
            warning('build_coupled_lin_sys:SingularTheta', ...
                'sin(theta_1) = %.2e for pili %d -- Euler angle singularity.', sT, k);
            sT = sign(sT + 1e-20) * 1e-12;
        end

        r1 = row_offset + 1;
        r2 = row_offset + 2;
        r3 = row_offset + 3;

        omega_cols = 4:6;
        theta1_col = col_offset + 1;
        phi1_col   = col_offset + N_k + 1;
        psi1_col   = col_offset + 2*N_k + 1;

        J_euler_basis = [-sP,         cP,        0;
                         -cT*cP/sT,  -cT*sP/sT, 1;
                          cP/sT,      sP/sT,     0];
        J_euler = J_euler_basis * rot_k;

        rc_k = 'clamped';
        if isfield(params, 'root_coupling') && numel(params.root_coupling) >= k
            rc_k = params.root_coupling{k};
        end

        if strcmp(rc_k, 'clamped')
            lin_sys(r1, :) = 0;
            lin_sys(r2, :) = 0;
            lin_sys(r3, :) = 0;

            lin_sys(r1, theta1_col) = 1;   lin_sys(r1, omega_cols) = -J_euler(1,:);
            lin_sys(r2, phi1_col)   = 1;   lin_sys(r2, omega_cols) = -J_euler(2,:);
            lin_sys(r3, psi1_col)   = 1;   lin_sys(r3, omega_cols) = -J_euler(3,:);

        elseif strcmp(rc_k, 'free_torsion')
            lin_sys(r1, :) = 0;
            lin_sys(r2, :) = 0;

            lin_sys(r1, theta1_col) = 1;   lin_sys(r1, omega_cols) = -J_euler(1,:);
            lin_sys(r2, phi1_col)   = 1;   lin_sys(r2, omega_cols) = -J_euler(2,:);

        else
            error('build_coupled_lin_sys:UnknownRootCoupling', ...
                ['Pili %d: root_coupling = ''%s''. ', ...
                 'Supported: clamped, free_torsion.'], k, rc_k);
        end

        geom.J_euler3{k} = J_euler(3,:);

        col_offset = col_offset + 3 * N_k;
        row_offset = row_offset + 3 * N_k;
    end
end


function logf_h = get_params_logf(params)
%GET_PARAMS_LOGF  Retrieve logf handle from params, or fallback to fprintf
    if isstruct(params) && isfield(params, 'logf') && ~isempty(params.logf)
        logf_h = params.logf;
    else
        logf_h = @(varargin) fprintf(varargin{:});
    end
end


function tf = should_run_full_solve_check(params, force_strict)
%SHOULD_RUN_FULL_SOLVE_CHECK  Cadence gate for expensive multi-RHS residual checks
%  Default ('periodic', stride 50) checks the first call and every 50th call
%  in full while keeping finiteness checks on every call.  Passing
%  params.solve_check_mode = 'strict' or force_strict = true checks every call.
%  Tracks params.run_token so a new simulation in the same MATLAB session
%  always resets the counter and checks its first call.
    persistent call_counter last_run_token
    if force_strict
        tf = true;
        return;
    end
    check_mode = 'periodic';
    if isstruct(params) && isfield(params, 'solve_check_mode') && ~isempty(params.solve_check_mode)
        check_mode = params.solve_check_mode;
    end
    if strcmpi(check_mode, 'strict')
        tf = true;
        return;
    end
    stride = 50;
    if isstruct(params) && isfield(params, 'solve_check_stride') && ~isempty(params.solve_check_stride)
        stride = max(1, round(double(params.solve_check_stride)));
    end
    if isstruct(params) && isfield(params, 'run_token') && ~isempty(params.run_token)
        if isempty(last_run_token) || ~isequal(last_run_token, params.run_token)
            call_counter   = 0;
            last_run_token = params.run_token;
        end
    end
    if isempty(call_counter)
        call_counter = 0;
    end
    call_counter = call_counter + 1;
    tf = (call_counter == 1) || (mod(call_counter, stride) == 0);
end
