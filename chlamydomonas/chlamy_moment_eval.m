function M = chlamy_moment_eval(par, t, N_k)
%CHLAMY_MOMENT_EVAL  Evaluate the uniform-wave Chlamydomonas active moment
%
%  Description:
%    Stand-alone evaluator behind chlamy_moment_fun.  The drive data travel
%    in the struct PAR, which is stored by value inside the anonymous handle
%    returned by chlamy_moment_fun, so a saved simulation_results.mat keeps
%    working after this file is edited.  See chlamy_moment_fun for the
%    definition of q(s,t) and of the envelopes.
%
%  Syntax:
%    M = chlamy_moment_eval(par, t, N_k)
%
%  Inputs:
%    par - struct with .q_star, .w_star, .kw, .phi, .s_taper, .standing,
%          .root_env ([] or struct)
%    t   - Dimensionless time (t_phys / T_0)
%    N_k - Number of balance nodes of the flagellum
%
%  Output:
%    M - [3 x N_k] internal active moment (bending about d2 only)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11

q_star   = par.q_star;
w_star   = par.w_star;
kw       = par.kw;
phi      = par.phi;
s_taper  = par.s_taper;
standing = par.standing;
root_env = par.root_env;

s   = (0 : N_k - 1)' / N_k;                 % balance nodes i = 1..N_k
env = ones(size(s));
tail = s > s_taper;
env(tail) = 0.5 * (1 + cos(pi * (s(tail) - s_taper) / (1 - s_taper)));
if ~isempty(root_env)
    switch root_env.type
        case 'gauss'
            env = env .* (1 - exp(-(root_env.kE * s).^2));
        case 'linear'
            ramp = ones(size(s));
            up = s < root_env.s_p;
            ramp(up) = root_env.root_frac + (1 - root_env.root_frac) * s(up) / root_env.s_p;
            env = env .* ramp;
        case 'table'
            % measured envelope: pchip through (.s, .val), end values held;
            % values may exceed 1 (proximal boost)
            env = env .* interp1([0; root_env.s(:); 1], ...
                [root_env.val(1); root_env.val(:); root_env.val(end)], s, 'pchip');
        otherwise
            error('chlamy_moment_fun:RootEnv', 'Unknown root_env.type ''%s''.', root_env.type);
    end
end
if standing
    % reciprocal control: fixed spatial shape, oscillating in time
    q = q_star * cos(phi) * env .* cos(kw * s) * cos(w_star * t);
else
    q = q_star * env .* cos(kw * s - w_star * t + phi);
end
M = [zeros(1, N_k); q.'; zeros(1, N_k)];    % bending about d2 only
end
