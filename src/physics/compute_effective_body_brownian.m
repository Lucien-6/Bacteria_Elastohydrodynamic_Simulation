function [L_eff, D_eff, M_eff_dim] = compute_effective_body_brownian(Z, params)
%COMPUTE_EFFECTIVE_BODY_BROWNIAN  Effective 6x6 body-only Brownian parameters
%
%  Description:
%    Computes the effective 6x6 body mobility, diffusion tensor, and
%    Cholesky factor by probing the fully coupled body-pili linear
%    system.  The resulting mobility implicitly incorporates the
%    hydrodynamic contribution of ALL pili via the Schur complement
%    of the pili DOFs, but stochastic forcing is subsequently applied only
%    to the body pose; pili internal modes remain deterministic:
%
%      M_eff = [Xi_gen^{-1}]_{BB}
%            = (Xi_BB - Xi_BP * Xi_PP^{-1} * Xi_PB)^{-1}
%
%    This is obtained without explicit block inversion by solving
%    the full lin_sys with 6 unit probe vectors and extracting the
%    body rows of the solution.  This routine is deterministic: it
%    constructs the Brownian covariance and factor but does not draw
%    random numbers or advance any random-number generator.
%
%  Algorithm:
%    1. Assemble the coupled linear system via build_coupled_lin_sys.
%    2. Probe with 6 unit force/torque vectors on the body DOFs.
%    3. Extract the dimensionless 6x6 effective mobility, measure its
%       relative asymmetry before modification, and emit diagnostics
%       when the mixed-discretisation asymmetry envelope is exceeded.
%    4. Symmetrise only the mobility used for the Brownian
%       fluctuation-dissipation path:
%       M_eff <- (M_eff + M_eff')/2.
%       The deterministic ode_coupled_rhs solve continues to use the
%       original, unsymmetrised hydrodynamic operator.
%    5. Project the dimensionless mobility to SPD when round-off requires
%       it.  A materially negative eigenvalue is treated as an assembly error.
%    6. Dimensionalise using a congruence transformation.
%    7. Form D_eff = kB*T*M_eff_dim and map the dimensionless Cholesky
%       factor to physical coordinates.  No mixed-unit eigendecomposition
%       is performed and no stochastic increment is sampled here.
%
%  Diagnostics and Warnings:
%    Emits compute_effective_body_brownian:AsymmetryHigh when relative
%    asymmetry exceeds 1e-3; a further log note is written above 1e-2.
%    Emits compute_effective_body_brownian:IllConditioned when the
%    dimensionless effective mobility is nonfinite or has rcond < 1e-10.
%    Emits compute_effective_body_brownian:NotPosDef when the dimensionless
%    mobility requires a round-off-scale SPD projection.  A materially
%    indefinite mobility raises compute_effective_body_brownian:Indefinite.
%
%  Syntax:
%    [L_eff, D_eff, M_eff_dim] = compute_effective_body_brownian(Z, params)
%
%  Inputs:
%    Z      - [(7 + 3*sum(N_segments)) x 1] Global state vector
%             (dimensionless)
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    L_eff     - [6x6] Lower Cholesky factor of 2*D_eff
%                Units: [m/sqrt(s), m/sqrt(s), m/sqrt(s),
%                        rad/sqrt(s), rad/sqrt(s), rad/sqrt(s)]
%    D_eff     - [6x6] Full effective diffusion tensor
%                Units: [m^2/s (TT), m*rad/s (TR), rad^2/s (RR)]
%    M_eff_dim - [6x6] Dimensional effective mobility matrix
%                Units: [m/(s*N) (TT), 1/(s*N) (TR), rad/(s*N*m) (RR)]
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] Cortez et al. (2005), Phys. Fluids 17, 031504.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-31
%  History: 2026-05-06 — rcond diagnostic moved to M_eff_dimless
%           (project v7.4.1; see CHANGELOG.md).
%  Updated: 2026-07-13

%% Step 1: Assemble coupled linear system
[lin_sys, ~] = build_coupled_lin_sys(Z, params);

n_state = size(lin_sys, 1);

% Retrieve logger handle (if available in params)
logf_h = get_brown_logf(params);

%% Step 2: Probe with unit body forces/torques
%  lin_sys = -Xi_gen, so to get Xi_gen^{-1} * e_i we solve:
%    lin_sys * z_dot = -e_i
%  Then z_dot = -Xi_gen^{-1} * (-e_i) = Xi_gen^{-1} * e_i = M_gen * e_i
R_probe = zeros(n_state, 6);
R_probe(1:6, 1:6) = -eye(6);

z_dot_probe = lin_sys \ R_probe;
check_linear_solve(lin_sys, z_dot_probe, R_probe, ...
    'effective body mobility probe solve');

%% Step 3: Extract dimensionless 6x6 effective mobility
M_eff_dimless = z_dot_probe(1:6, :);

%  Retrieve regularisation lengths for the asymmetry diagnostic below.
%  eps_body may be absent in edge-case params structs; guard with isfield.
eps_body_val = NaN;
eps_pili_val = NaN;
if isfield(params, 'eps_body')
    eps_body_val = params.eps_body;
end
if isfield(params, 'epsilon')
    eps_pili_val = params.epsilon;
end

%  Diagnostic: measure asymmetry BEFORE symmetrisation.
%  For pure-RSM (NoPili) systems, LU round-off on a ~90-dim matrix yields
%  rel_asym ~ 1e-12..1e-10.  For RSM+RSS mixed configurations (eps_body !=
%  eps_pili), M_weighted is structurally non-symmetric; the Schur-complement
%  extraction attenuates this to rel_asym ~ 1e-5..1e-3, which is an expected
%  numerical consequence of the mixed-epsilon regularisation, not a code bug.
%  A value >= 1e-3 now triggers the warning (truly structural assembly bugs
%  or catastrophic conditioning produce rel_asym >> 1e-3).  The numerical
%  path is unchanged; this block only emits a warning if something is
%  clearly wrong beyond the expected mixed-epsilon asymmetry envelope.
asym_abs  = norm(M_eff_dimless - M_eff_dimless', 'fro');
sym_scale = max(norm(M_eff_dimless, 'fro'), eps);
rel_asym  = asym_abs / sym_scale;
if rel_asym > 1e-3
    logf_h(['[BROWN] WARNING: M_eff_dimless relative asymmetry = %.2e (expected <1e-3 for ' ...
         'RSM+RSS mixed epsilon; pure RSM yields <1e-10). Possible lin_sys assembly bug or ' ...
         'extreme conditioning; silent symmetrisation will still be applied.\n'], ...
        rel_asym);
    warning('compute_effective_body_brownian:AsymmetryHigh', ...
        ['M_eff_dimless relative asymmetry = %.2e (expected <1e-3 for ', ...
         'RSM+RSS mixed epsilon; pure RSM yields <1e-10). Possible lin_sys assembly bug or ', ...
         'extreme conditioning; silent symmetrisation will still be applied.'], ...
        rel_asym);
end
if rel_asym > 1e-2
    %  Quantitative guidance: at this level the deterministic/stochastic operator
    %  splitting (see docs/theory.md §10 "Known approximation") may carry a
    %  non-negligible bias in the equilibrium distribution.  Suggest the user
    %  consider tightening eps_body toward eps_pili (or vice versa), provided
    %  N_body is adjusted to keep the RSM leak within tolerance (§5.1).
    logf_h(['[BROWN] NOTE: rel_asym = %.2e > 1e-2. The deterministic/stochastic ' ...
         'operator-splitting bias may be non-negligible. Consider tightening eps_body ' ...
         'toward eps_pili (currently eps_body=%.4e, eps_pili=%.4e), provided N_body ' ...
         'is adjusted to keep the RSM leak within tolerance (see theory.md sec 5.1).\n'], ...
        rel_asym, eps_body_val, eps_pili_val);
end

% Symmetrise the mobility used by the Brownian covariance.  The deterministic
% coupled solve still uses the original hydrodynamic operator; this path
% constructs the fluctuation-dissipation tensor for the body-only stochastic
% increment.
M_eff_dimless = 0.5 * (M_eff_dimless + M_eff_dimless');

% Perform all spectral repair in dimensionless generalized coordinates.
% Translational and rotational entries of the dimensional covariance have
% different units, so an eigendecomposition of D_eff is not invariant under
% a change of length unit.  The dimensionless mobility has one coherent
% numerical scale and is the correct object on which to enforce SPD.
[M_eff_dimless, was_projected, min_eig_before] = ...
    project_spd_dimensionless(M_eff_dimless);
if was_projected
    logf_h(['[BROWN] WARNING: dimensionless effective mobility required an ', ...
        'SPD projection (min eigenvalue = %.2e).\n'], min_eig_before);
    warning('compute_effective_body_brownian:NotPosDef', ...
        ['Dimensionless effective mobility was not strictly positive ', ...
         'definite (min eigenvalue = %.2e); a round-off-scale SPD ', ...
         'projection was applied.'], min_eig_before);
end

%  Diagnostic: condition number of the DIMENSIONLESS effective mobility.
%  We deliberately probe `M_eff_dimless` rather than the dimensional
%  `M_eff_dim` because the dimensional matrix mixes TT/TR/RR blocks at
%  vastly different scales -- the (S * S') outer product divides the TT
%  block by 8*pi*mu*L_0 (~1e-8) and the RR block by 8*pi*mu*L_0^3
%  (~1e-20), so even a perfectly healthy mobility carries a 1-norm
%  condition number of order L_0^{-2} (~1e12 for L_0 ~ 1e-6 m), which
%  triggered spurious "near-singular" warnings on every Brownian step.
%  The dimensionless form is the physically meaningful diagnostic:
%  healthy body-pili configurations yield cond(M_eff_dimless) ~ 1..1e3,
%  and rcond < 1e-10 here genuinely signals near-degenerate TT/RR blocks
%  (extreme aspect ratio, singular clamped BC, or catastrophic mobility
%  cancellation) that would compromise Cholesky stability downstream.
%  Cholesky itself remains numerically stable on the dimensional matrix
%  because it depends on eigenvalue spread, not on the row-norm-based
%  1-norm rcond used by MATLAB.  rcond cost is O(1) for 6x6.
rc_mobility = rcond(M_eff_dimless);
if ~isfinite(rc_mobility) || rc_mobility < 1e-10
    logf_h(['[BROWN] WARNING: M_eff_dimless rcond = %.2e (< 1e-10): effective mobility is ' ...
         'nearly singular, Cholesky factor may be unreliable.\n'], rc_mobility);
    warning('compute_effective_body_brownian:IllConditioned', ...
        ['M_eff_dimless rcond = %.2e (< 1e-10): effective mobility is ', ...
         'nearly singular, Cholesky factor may be unreliable.'], rc_mobility);
end

%% Step 4: Dimensionalise
%  Scaling: the dimensionless system uses EI-based nondimensionalisation.
%  M_eff_dim(j,i) = M_eff_dimless(j,i) / (S_j * S_i)
%  where S encodes the force/velocity scaling per DOF component.
%
%  Derivation:
%    V_dim = V~ * U_0,     F_dim = F~ * F_0     => M_TT_dim = M_TT~ * U_0/F_0 = M_TT~ / (8*pi*mu*L_0)
%    Omega_dim = Omega~/T_0, tau_dim = tau~*F_0*L_0 => M_RR_dim = M_RR~ / (8*pi*mu*L_0^3)
%    V_dim / tau_dim => M_TR_dim = M_TR~ / (8*pi*mu*L_0^2)
mu  = params.dim.mu;
L_0 = params.scales.L_0;

coeff = 8 * pi * mu;
S = [sqrt(coeff * L_0) * ones(3,1); sqrt(coeff * L_0^3) * ones(3,1)];
M_eff_dim = M_eff_dimless ./ (S * S');

%  Note: no separate symmetrisation or spectral projection of M_eff_dim is
%  needed here.
%  M_eff_dimless was symmetrised at Step 3 (line above the rcond check),
%  and the Hadamard division by the symmetric outer product (S * S')
%  preserves symmetry:  (A ./ B)_ij = A_ij / B_ij, so if A_ij = A_ji and
%  B_ij = B_ji then (A ./ B)_ij = (A ./ B)_ji.  The downstream
%  M_eff_dim is already symmetric positive definite by congruence.

%% Step 5: Diffusion tensor and Cholesky decomposition
kB    = 1.380649e-23;
Temper = params.dim.Temper;
if ~(isfinite(Temper) && Temper > 0)
    error('compute_effective_body_brownian:InvalidTemperature', ...
        'Brownian temperature must be finite and positive, got %.4g K.', Temper);
end

D_eff = kB * Temper * M_eff_dim;

% Factor in dimensionless coordinates, then apply the physical row scaling.
% This is algebraically equivalent to chol(2*D_eff) but avoids asking a
% mixed-unit matrix to determine its own numerical spectral scale.
L_dimless = chol(2 * kB * Temper * M_eff_dimless, 'lower');
L_eff = L_dimless ./ S;

end


%% ==================== Local Helper Function ====================

function logf_h = get_brown_logf(params)
%GET_BROWN_LOGF  Retrieve logf handle from params, or fallback to fprintf
    if isstruct(params) && isfield(params, 'logf') && ~isempty(params.logf)
        logf_h = params.logf;
    else
        logf_h = @(varargin) fprintf(varargin{:});
    end
end


function [A_spd, was_projected, min_eig_before] = project_spd_dimensionless(A)
%PROJECT_SPD_DIMENSIONLESS  SPD repair in dimensionless coordinates only.
%  Round-off-scale non-positive eigenvalues are lifted to a relative
%  machine-precision floor.  A materially negative eigenvalue indicates a
%  non-dissipative hydrodynamic operator and is never silently hidden.
    A_sym = 0.5 * (A + A');
    [V, D] = eig(A_sym);
    lam = real(diag(D));
    min_eig_before = min(lam);
    spectral_scale = max(max(abs(lam)), realmin('double'));
    material_tol = 1e-10 * spectral_scale;
    if min_eig_before < -material_tol
        error('compute_effective_body_brownian:Indefinite', ...
            ['Dimensionless effective mobility is materially indefinite: ', ...
             'lambda_min = %.3e, spectral scale = %.3e.'], ...
            min_eig_before, spectral_scale);
    end

    floor_eig = max(100 * size(A_sym, 1) * eps(spectral_scale), ...
                    realmin('double'));

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
