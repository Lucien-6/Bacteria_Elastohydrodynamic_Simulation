function fun = chlamy_moment_fun(Q_Nm, f_Hz, lambda_rel, phi, s_taper, EI, L_k, T_0, standing, root_env)
%CHLAMY_MOMENT_FUN  Custom active-moment function for a Chlamydomonas flagellum
%
%  Description:
%    Returns the function handle required by active_moment.type = 'custom'.
%    The INTERNAL active bending moment carried across the cross-section at
%    arc length s (the quantity that the 'traveling_wave' mode builds as an
%    integral from s to the tip) is prescribed directly as
%
%        q(s,t) = Q * env(s) * cos( 2*pi*s/lambda_rel - w*t + phi )
%
%    with a cosine taper env(s) that is 1 for s <= s_taper and falls
%    smoothly to 0 at the tip (the internal moment vanishes at a free end):
%
%        env(s) = 1                                      , s <= s_taper
%        env(s) = 0.5*(1 + cos(pi*(s-s_taper)/(1-s_taper))), s >  s_taper
%
%    Compared with the built-in 'traveling_wave' mode, the envelope has no
%    spatially uniform, time-dependent offset (which produced a hump-shaped
%    amplitude and tied the wavelength to 4*pi/kw), and the wavelength is
%    set directly through lambda_rel = lambda / L.
%
%    The moment acts about d2 (in-plane bending of a planar centreline with
%    psi = 0).  Sign convention is the same as 'traveling_wave' with
%    d2_amp = 1: positive Q and positive phase velocity (base -> tip).
%
%  Syntax:
%    fun = chlamy_moment_fun(Q_Nm, f_Hz, lambda_rel, phi, s_taper, EI, L_k, T_0)
%    fun = chlamy_moment_fun(Q_Nm, f_Hz, lambda_rel, phi, s_taper, EI, L_k, T_0, standing)
%    fun = chlamy_moment_fun(..., standing, root_env)
%
%  Inputs:
%    Q_Nm       - Peak internal active moment [N m]
%    f_Hz       - Beat frequency [Hz]
%    lambda_rel - Phase wavelength / flagellar length [-]
%    phi        - Phase offset [rad] (pi for the mirror-image flagellum)
%    s_taper    - Start of the tip taper on s in [0,1]
%    EI         - Bending stiffness [N m^2]
%    L_k        - Flagellar length [m]
%    T_0        - Time scale of the code [s]
%    standing   - (Optional) true: standing-wave control
%                 q = Q env(s) cos(phi) cos(kw s) cos(w t), a reciprocal
%                 (time-reversible in the Stokes limit) drive; default false
%    root_env   - (Optional) struct describing a proximal ramp that is
%                 multiplied onto env(s) (the basal region carries no dynein);
%                 default [] = no ramp (env = 1 at the root).  Not a root
%                 spring: the root boundary condition is unchanged (clamped).
%                   .type = 'gauss'  : 1 - exp(-(kE*s)^2), kE = .kE (2*pi)
%                   .type = 'linear' : root_frac + (1-root_frac)*s/s_p for
%                                      s < s_p, else 1 (.root_frac, .s_p)
%                   .type = 'table'  : pchip through the points (.s, .val)
%                                      (e.g. an envelope of the measured
%                                      curvature amplitude, with a proximal
%                                      boost)
%                 s is the relative arc length in [0,1].
%
%  Output:
%    fun - Handle fun(t, k, N_k, d1, d2, d3, d1N, d2N, d3N) -> [3 x N_k]
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-06
%  Updated: 2026-11-11 (anonymous handle + chlamy_moment_eval, so saved
%                       results survive later edits of this file)

if nargin < 9 || isempty(standing), standing = false; end
if nargin < 10, root_env = []; end
q_star = Q_Nm * L_k / EI;            % same non-dimensionalisation as amplitude_dimless
w_star = 2 * pi * f_Hz * T_0;
kw     = 2 * pi / lambda_rel;

% All drive data are stored by value in PAR.  The handle is an anonymous
% function that calls the stand-alone chlamy_moment_eval, NOT a nested
% function: a nested handle saved in simulation_results.mat keeps only a
% snapshot of this workspace and breaks as soon as this file is edited.
par = struct('q_star', q_star, 'w_star', w_star, 'kw', kw, 'phi', phi, ...
    's_taper', s_taper, 'standing', logical(standing), 'root_env', {root_env});
fun = @(t, ~, N_k, ~, ~, ~, ~, ~, ~) chlamy_moment_eval(par, t, N_k);

    function M = moment_fun(t, ~, N_k, ~, ~, ~, ~, ~, ~)
        % LEGACY PATH - keep the name and the q_star/... closure so that
        % results saved before 2026-11-11 can still be evaluated (or have
        % their workspace read by restore_saved_moment_handles).  Newer runs
        % never call this path.
        par_legacy = struct('q_star', q_star, 'w_star', w_star, 'kw', kw, ...
            'phi', phi, 's_taper', s_taper, 'standing', logical(standing), ...
            'root_env', {root_env});
        M = chlamy_moment_eval(par_legacy, t, N_k);
    end
end
