function [force_over_mu_U, torque_over_mu_U] = rft_helix_load(X, U, lambda_over_R, a_over_R)
%RFT_HELIX_LOAD  Gray-Hancock resistive-force load on a helical centreline
%
%  Description:
%    Local drag per unit length with the Gray and Hancock coefficients
%
%        c_t = 2*pi*mu / (log(2*lambda/a) - 1/2)
%        c_n = 4*pi*mu / (log(2*lambda/a) + 1/2)
%
%    evaluated at segment midpoints.  This is the curve Rodenborn et al.
%    show departing from their measurements at short pitch.  It is a
%    negative control, not the pass line of the helical check.
%
%    lambda_over_R and a_over_R use the same length unit as X.  Outputs
%    match rss_line_load.
%
%  Syntax:
%    [force_over_mu_U, torque_over_mu_U] = rft_helix_load(X, U, lambda_over_R, a_over_R)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-15

dX = diff(X, 1, 1);
ds = sqrt(sum(dX.^2, 2));
t = dX ./ ds;
mid = 0.5 * (X(1:end-1, :) + X(2:end, :));
Um = 0.5 * (U(1:end-1, :) + U(2:end, :));
log_arg = log(2 * lambda_over_R / a_over_R);
c_t = 2 * pi / (log_arg - 0.5);
c_n = 4 * pi / (log_arg + 0.5);

vt = sum(Um .* t, 2);
f = c_n * Um + (c_t - c_n) * (vt .* t);
force_over_mu_U = (ds.' * f).';
moment = cross(mid, f, 2);
torque_over_mu_U = (ds.' * moment).';
end
