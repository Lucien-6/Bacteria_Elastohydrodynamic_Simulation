function [force_over_mu_U, torque_over_mu_U] = rss_body_line_load(X_body, X_line, U_body, U_line, eps_body, eps_line)
%RSS_BODY_LINE_LOAD  Rigid load on a regularised-Stokeslet body plus a filament
%
%  Description:
%    Same kernel and 8*pi normalisation as rss_line_load, with a surface
%    cloud of point forces on the body (regularised Stokeslets, Cortez 2005).
%    Body unknowns are point forces.  Filament unknowns are nodal force
%    densities, integrated with the trapezoidal rule.  The torque is about
%    the origin of the coordinates passed in.
%
%  Syntax:
%    [force_over_mu_U, torque_over_mu_U] = rss_body_line_load( ...
%        X_body, X_line, U_body, U_line, eps_body, eps_line)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-15

n_body = size(X_body, 1);
n_node = size(X_line, 1);
n_seg = n_node - 1;
dX = diff(X_line, 1, 1);
seg_len = sqrt(sum(dX.^2, 2));
d3 = dX ./ seg_len;
X = [X_body; X_line];
geom = struct('d3', {{d3}});
params = struct('N_body', n_body, 'M_pili', 1, 'N_segments', n_seg, ...
    'eps_body', eps_body, 'epsilon', eps_line);
M = build_global_resistance(X, params, geom);
u = [reshape(U_body.', [], 1); reshape(U_line.', [], 1)];
f = reshape(M \ u, 3, []).';
f_body = f(1:n_body, :);
f_line = f(n_body + (1:n_node), :);

w = zeros(n_node, 1);
w(1) = seg_len(1) / 2;
w(end) = seg_len(end) / 2;
if n_node > 2
    w(2:end-1) = (seg_len(1:end-1) + seg_len(2:end)) / 2;
end
force_over_mu_U = (8 * pi) * (sum(f_body, 1) + (w.' * f_line));
moment = cross(X_body, f_body, 2);
moment_line = cross(X_line, f_line, 2);
torque_over_mu_U = (8 * pi) * (sum(moment, 1) + (w.' * moment_line));
end
