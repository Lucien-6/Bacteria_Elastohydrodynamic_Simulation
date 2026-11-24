function [f_body, g_modes, F_pili, info, tau_pili] = sample_brownian_coupled_force(Z, params, dt_B_s, rng_stream)
%SAMPLE_BROWNIAN_COUPLED_FORCE  Thermal forces on the body and the slow pilus modes
%
%  Description:
%    Pilus thermal-fluctuation scheme (v9.0.0): the counterpart of
%    sample_brownian_body_force when params.brown_pili is enabled.  One
%    Brownian interval of length dt_B is realised as a CONSTANT generalised
%    force on the thermal coordinates
%
%        y = [ X_cm (3); body rotation (3); a_1 .. a_n ] ,
%
%    where a_j are the amplitudes of the 2 * n_orders * M_pili slowest
%    bending modes of the pili (init_pili_brownian_modes).  With the joint
%    mobility M_joint frozen at the interval start and probed with every
%    non-thermal pilus mode held (compute_effective_coupled_brownian),
%
%        dY = L * eta * sqrt(dt_B),   L L' = 2 kT M_joint,   eta ~ N(0, I),
%        [f_body; g] = M_joint^{-1} * dY / dt_B ,
%
%    so that <[f;g][f;g]'> dt_B = 2 kT M_joint^{-1} (fluctuation-dissipation
%    for the thermal coordinates, including the body-mode and mode-mode
%    hydrodynamic cross-correlations).  The forces are held for the whole
%    interval by ode_coupled_rhs:
%
%        R_global(1:6) = R_global(1:6) - f_body
%        R_global      = R_global      - pili_thermal_rhs(F_pili, ...) .
%
%    F_pili is the force density on the pilus centrelines that realises
%    these forces at the interval start: the viscous load the pili carry in
%    the constrained motion driven by [f_body; g]
%    (compute_effective_coupled_brownian, info.P_load).  It exerts the modal
%    forces g and, for the body part, keeps the non-thermal pilus modes at
%    rest while the body is driven, so the response of the integration at
%    the interval start is M_joint * [f_body; g] and no deterministic mode
%    is kicked by the noise.  It is held in the body frame; its net force
%    and torque on the cell are not applied twice (pili_thermal_rhs zeroes
%    rows 1:6).  The modal force is deliberately not applied as a
%    material-frame moment: over an interval much longer than the twist
%    relaxation time (tens of microseconds) a held moment of white-noise
%    amplitude winds the filament into a helix, whereas a force on the
%    centreline of a naturally straight isotropic filament cannot twist it.
%
%    Faster pilus modes and twist receive no noise.  The elastic restoring
%    moments are unchanged, so each thermal mode relaxes to equipartition,
%    <a_j^2 / 2> = kT / 2, through the same coupled dynamics as before.  As
%    in the body-only scheme the Ito drift kT div(M_joint) is generated
%    implicitly by the mobility sampled along the trajectory.
%
%    That mechanism presumes a constant generalised force.  A centreline
%    load held in the body frame exerts a modal force that changes as the
%    filament moves, which leaves, averaged over the thermal forces, the
%    spurious modal force kT * c, c_j = sum_k d g_j(P_k) / d a_k (a few per
%    cent of the elastic force at 293 K; at the straight rest state it is
%    zero for a pilus on the body axis, not for off-axis pili).
%    The deterministic modal force -kT * c is therefore added to g for the
%    interval; without it the second bending order of a single pilus
%    settles about 5 % below equipartition at 293 K.
%
%    Holding the force over dt_B lowers the stationary variance of a mode of
%    relaxation rate lambda by the factor (2/x) tanh(x/2), x = lambda dt_B
%    (0.5 % at x = 0.25, 2 % at x = 0.5); keep dt_B below a quarter of the
%    fastest thermal relaxation time.
%
%    This routine advances the Brownian stream by exactly 6 + n normal
%    deviates per call.  If params.brown_pili_streams is a cell of one
%    RandStream per pilus, the 6 body deviates still come from rng_stream
%    and each pilus stream supplies the curvature-coordinate noise of that
%    filament only.  Those draws are projected with the K-orthonormal mode
%    matrix, whose blocks sum to the identity, so the modal vector stays
%    standard normal and the joint covariance is unchanged.  In that mode
%    rng_stream advances by exactly 6 deviates, the same call
%    sample_brownian_body_force makes, and each pilus stream advances by
%    its own number of curvature coordinates.
%
%  Syntax:
%    [f_body, g_modes, F_pili, info] = ...
%        sample_brownian_coupled_force(Z, params, dt_B_s, rng_stream)
%
%  Inputs:
%    Z          - [(7 + 3*sum(N_segments)) x 1] State at the interval start
%    params     - Parameter structure with params.brown_pili initialised
%    dt_B_s     - Interval length [s] (> 0)
%    rng_stream - RandStream used for the normal deviates
%
%  Outputs:
%    f_body   - [6x1] Body force (F_0) and torque (F_0 L_0), dimensionless
%    g_modes  - [n x 1] Modal generalised forces applied over the interval
%               (F_0 L_0), dimensionless: thermal part plus the drift
%               compensation info.g_drift
%    F_pili   - [3 x N_nodes] Dimensionless force density at the pilus nodes
%               (all pili, order of X_all) in the body frame; goes into
%               params.brown_pili_force
%    info     - Diagnostics: .eta, .dU_target_dim [m; rad], .f_th_dim
%               [N; N*m], .da_target (dimensionless), .dt_B_s, .M_joint,
%               .kT_dimless, .rel_asym, .pattern_resid,
%               .hold_pattern_resid, .twist_row_resid, .hold_force,
%               .g_drift [n x 1] (deterministic modal force -kT * c)
%
%  References:
%    [1] Ermak & McCammon (1978), J. Chem. Phys. 69, 1352.
%    [2] Fixman (1978), J. Chem. Phys. 69, 1527.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11

if ~(isscalar(dt_B_s) && isfinite(dt_B_s) && dt_B_s > 0)
    error('sample_brownian_coupled_force:InvalidInterval', ...
        'Brownian interval length must be a finite positive scalar, got %.4g s.', dt_B_s);
end

L_0 = params.scales.L_0;
T_0 = params.scales.T_0;
F_0 = params.scales.F_0;

kB     = 1.380649e-23;
Temper = params.dim.Temper;
if ~(isfinite(Temper) && Temper > 0)
    error('sample_brownian_coupled_force:InvalidTemperature', ...
        'Brownian temperature must be finite and positive, got %.4g K.', Temper);
end
kT_dimless = kB * Temper / (F_0 * L_0);

%% Joint mobility at the interval start
[M_joint, m_info] = compute_effective_coupled_brownian(Z, params);
n_y = size(M_joint, 1);

%% Target increment of the thermal coordinates (dimensionless)
%  Default: one joint stream.  With per-pilus streams, the body block is
%  drawn from rng_stream (6 deviates, paired with the body-only scheme)
%  and each filament contributes only through its own stream.
eta        = draw_thermal_eta(rng_stream, n_y, params);
dt_dimless = dt_B_s / T_0;
L_joint    = chol(2 * kT_dimless * M_joint, 'lower');
dY         = L_joint * eta * sqrt(dt_dimless);

%% Constant generalised force that realises it on the frozen system
force = M_joint \ (dY / dt_dimless);
if any(~isfinite(force))
    error('sample_brownian_coupled_force:NonFiniteForce', ...
        'Thermal generalised force is non-finite; the joint mobility may be singular.');
end

%% Compensation of the held-load drift
%  Holding the centreline loads in the body frame leaves, averaged over the
%  thermal forces, the spurious modal force kT * c
%  (compute_effective_coupled_brownian, info.load_drift).  It is removed by
%  a deterministic modal force -kT * c, applied through the same loads.
g_drift = -kT_dimless * m_info.load_drift;

f_body  = force(1:6);
g_modes = force(7:end) + g_drift;

%% Centreline force density and balancing axial torque density
%  Modal forces and, for the body force / torque, the load that keeps the
%  non-thermal pilus modes at rest while the body is driven.
fg_total = [f_body; g_modes];
F_pili   = reshape(m_info.P_load * fg_total, 3, []);
if any(~isfinite(F_pili(:)))
    error('sample_brownian_coupled_force:NonFiniteForce', ...
        'Thermal pilus force density is non-finite.');
end
if isfield(m_info, 'P_tau') && ~isempty(m_info.P_tau)
    tau_pili = m_info.P_tau * fg_total;
    if any(~isfinite(tau_pili(:)))
        error('sample_brownian_coupled_force:NonFiniteForce', ...
            'Thermal pilus axial torque density is non-finite.');
    end
else
    tau_pili = zeros(0, 1);
end

%% Diagnostics
info = struct();
info.eta                = eta;
info.dU_target_dim      = [dY(1:3) * L_0; dY(4:6)];
info.f_th_dim           = [f_body(1:3) * F_0; f_body(4:6) * F_0 * L_0];
info.da_target          = dY(7:end);
info.dt_B_s             = dt_B_s;
info.M_joint            = M_joint;
info.kT_dimless         = kT_dimless;
info.rel_asym           = m_info.rel_asym;
info.pattern_resid      = m_info.pattern_resid;
info.hold_pattern_resid = m_info.hold_pattern_resid;
info.twist_row_resid    = m_info.twist_row_resid;
info.hold_force         = m_info.hold_force;
info.g_drift            = g_drift;
info.tau_pili           = tau_pili;

end


function eta = draw_thermal_eta(rng_stream, n_y, params)
%DRAW_THERMAL_ETA  Standard normal vector for the joint thermal coordinates
if ~isfield(params, 'brown_pili_streams') || isempty(params.brown_pili_streams)
    eta = randn(rng_stream, n_y, 1);
    return;
end
bp = params.brown_pili;
n_modes = n_y - 6;
streams = params.brown_pili_streams;
if ~iscell(streams) || numel(streams) ~= params.M_pili || n_modes ~= bp.n_modes
    error('sample_brownian_coupled_force:BadPilusStreams', ...
        ['params.brown_pili_streams must hold one RandStream per pilus ', ...
         '(%d expected).'], params.M_pili);
end
eta_body  = randn(rng_stream, 6, 1);
eta_modes = zeros(n_modes, 1);
if isfield(bp, 'U_sel') && ~isempty(bp.U_sel)
    U = bp.U_sel;
else
    U = bp.M_meas ./ sqrt((bp.layout.wK .* bp.layout.wE)');
end
for k = 1:numel(streams)
    rows = (bp.layout.pilus == k);
    n_rows_k = sum(rows);
    if n_rows_k == 0
        continue;
    end
    v = randn(streams{k}, n_rows_k, 1);
    eta_modes = eta_modes + U(rows, :)' * v;
end
eta = [eta_body; eta_modes];
end
