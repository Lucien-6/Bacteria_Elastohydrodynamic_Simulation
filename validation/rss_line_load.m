function [force_over_mu_U, torque_over_mu_U] = rss_line_load(X, U, eps_len, origin)
%RSS_LINE_LOAD  Force and torque on a rigid centreline from the project RSS kernel
%
%  Description:
%    Prescribes a nodal velocity and returns the hydrodynamic load through
%    the same regularised-Stokeslet-segment matrix as build_global_resistance
%    (Cortez 2018).  No cell body is included.
%
%    X and eps_len are in one length unit L_u.  U is the velocity divided by
%    a velocity scale U_s.  The outputs are
%
%        force_over_mu_U  = F / (mu * U_s)          [L_u]
%        torque_over_mu_U = T / (mu * U_s)          [L_u^2]
%
%    with the torque about origin.  The 8*pi factor is the one absorbed by
%    the project's nondimensional kernel (see paper_static_checks, check S2).
%    Nodal unknowns are force densities; the line integral is trapezoidal.
%
%  Syntax:
%    [force_over_mu_U, torque_over_mu_U] = rss_line_load(X, U, eps_len)
%    [force_over_mu_U, torque_over_mu_U] = rss_line_load(X, U, eps_len, origin)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-15

if nargin < 4 || isempty(origin)
    origin = zeros(1, 3);
end
X = X - origin;
n_node = size(X, 1);
n_seg = n_node - 1;
if n_seg < 1
    error('rss_line_load:TooFewNodes', 'A centreline needs at least two nodes.');
end
dX = diff(X, 1, 1);
seg_len = sqrt(sum(dX.^2, 2));
if any(seg_len <= 0)
    error('rss_line_load:ZeroSegment', 'A centreline segment has zero length.');
end
d3 = dX ./ seg_len;

geom = struct('d3', {{d3}});
params = struct('N_body', 0, 'M_pili', 1, 'N_segments', n_seg, ...
    'eps_body', 1, 'epsilon', eps_len);
M = build_global_resistance(X, params, geom);
u = reshape(U.', [], 1);
f_node = reshape(M \ u, 3, []).';
w = trap_weights(seg_len);
force_over_mu_U = (8 * pi) * (w.' * f_node).';
moment_density = cross(X, f_node, 2);
torque_over_mu_U = (8 * pi) * (w.' * moment_density).';
end


function w = trap_weights(seg_len)
%TRAP_WEIGHTS  Trapezoidal line weights on nodes from segment lengths
n_node = numel(seg_len) + 1;
w = zeros(n_node, 1);
w(1) = seg_len(1) / 2;
w(end) = seg_len(end) / 2;
if n_node > 2
    w(2:end-1) = (seg_len(1:end-1) + seg_len(2:end)) / 2;
end
end
