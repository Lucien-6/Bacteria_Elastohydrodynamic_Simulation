function [L_eff, D_eff, M_eff_dim] = compute_effective_body_brownian(Z, params)
%COMPUTE_EFFECTIVE_BODY_BROWNIAN  Full-system effective 6x6 body Brownian parameters
%
%  Description:
%    Computes the effective 6x6 body mobility, diffusion tensor, and
%    Cholesky factor by probing the fully coupled body-pili linear
%    system.  The resulting mobility implicitly incorporates the
%    hydrodynamic contribution of ALL pili via the Schur complement
%    of the pili DOFs:
%
%      M_eff = [Xi_gen^{-1}]_{BB}
%            = (Xi_BB - Xi_BP * Xi_PP^{-1} * Xi_PB)^{-1}
%
%    This is obtained without explicit block inversion by solving
%    the full lin_sys with 6 unit probe vectors and extracting the
%    body rows of the solution.
%
%  Algorithm:
%    1. Assemble the coupled linear system via build_coupled_lin_sys.
%    2. Probe with 6 unit force/torque vectors on the body DOFs.
%    3. Extract dimensionless 6x6 effective mobility.
%    4. Dimensionalise using characteristic scales.
%    5. Compute diffusion tensor D = kB*T * M and Cholesky factor L.
%
%  Syntax:
%    [L_eff, D_eff, M_eff_dim] = compute_effective_body_brownian(Z, params)
%
%  Inputs:
%    Z      - [n x 1] Global state vector (dimensionless)
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
%  Updated: 2026-05-06 — Lucien — rcond diagnostic on M_eff_dimless (project v7.4.1; see CHANGELOG)

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

%% Step 3: Extract dimensionless 6x6 effective mobility
M_eff_dimless = z_dot_probe(1:6, :);

%  Diagnostic: measure asymmetry BEFORE symmetrisation.
%  LU round-off on a 96..186-dim system should yield rel_asym at 1e-12..1e-10.
%  A value >= 1e-6 indicates a structural issue (lin_sys assembly bug,
%  catastrophic conditioning, or mis-scaled quadrature weights) that
%  silent symmetrisation would otherwise mask.  The numerical path is
%  unchanged; this block only emits a warning if something is clearly wrong.
asym_abs  = norm(M_eff_dimless - M_eff_dimless', 'fro');
sym_scale = max(norm(M_eff_dimless, 'fro'), eps);
rel_asym  = asym_abs / sym_scale;
if rel_asym > 1e-6
    logf_h(['[BROWN] WARNING: M_eff_dimless relative asymmetry = %.2e (expected <1e-8 from ' ...
         'LU round-off alone). Possible lin_sys assembly bug or ' ...
         'extreme conditioning; silent symmetrisation will still be applied.\n'], ...
        rel_asym);
    warning('compute_effective_body_brownian:AsymmetryHigh', ...
        ['M_eff_dimless relative asymmetry = %.2e (expected <1e-8 from ', ...
         'LU round-off alone). Possible lin_sys assembly bug or ', ...
         'extreme conditioning; silent symmetrisation will still be applied.'], ...
        rel_asym);
end

% Symmetrise (numerical round-off may break exact symmetry)
M_eff_dimless = 0.5 * (M_eff_dimless + M_eff_dimless');

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

% Symmetrise dimensional mobility
M_eff_dim = 0.5 * (M_eff_dim + M_eff_dim');

%  NOTE: the rcond conditioning diagnostic was moved up to Step 3 (on the
%  dimensionless mobility) -- see the comment block above.  Probing the
%  dimensional `M_eff_dim` here is statistically meaningless because its
%  1-norm condition number is dominated by the TT/RR unit-scale ratio
%  (~L_0^{-2}) rather than by the underlying physical conditioning.

%% Step 5: Diffusion tensor and Cholesky decomposition
kB    = 1.380649e-23;
Temper = params.dim.Temper;

D_eff = kB * Temper * M_eff_dim;
D_eff = 0.5 * (D_eff + D_eff');

% Positive-definiteness check with diagonal regularisation fallback
[~, flag] = chol(2 * D_eff);
if flag ~= 0
    eigvals = eig(D_eff);
    min_eig = min(eigvals);
    logf_h('[BROWN] WARNING: 2*D_eff not positive definite (min eigenvalue = %.2e). Regularising...\n', min_eig);
    warning('compute_effective_body_brownian:NotPosDef', ...
        '2*D_eff not positive definite (min eigenvalue = %.2e). Regularising...', min_eig);
    reg = max(eps(norm(D_eff, 'fro')), abs(min_eig) * (1 + 1e-6));
    D_eff = D_eff + eye(6) * reg;
end

L_eff = chol(2 * D_eff, 'lower');

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
