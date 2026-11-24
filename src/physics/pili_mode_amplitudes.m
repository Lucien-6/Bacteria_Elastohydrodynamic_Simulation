function [a, energy_kT, elastic_kT] = pili_mode_amplitudes(Z, params)
%PILI_MODE_AMPLITUDES  Amplitudes of the thermal pilus bending modes
%
%  Description:
%    Projects the current pilus curvature onto the thermal modes selected by
%    init_pili_brownian_modes:
%
%        a = M_meas' * x ,     x = pili_bending_coords(Z, params) .
%
%    The elastic energy stored in mode j is a_j^2 / 2 in units of F_0 L_0.
%    In thermal equilibrium every mode holds kT / 2, i.e.
%    <a_j^2> = kT / (F_0 L_0).
%
%  Syntax:
%    a              = pili_mode_amplitudes(Z, params)
%    [a, energy_kT] = pili_mode_amplitudes(Z, params)
%    [a, energy_kT, elastic_kT] = pili_mode_amplitudes(Z, params)
%
%  Inputs:
%    Z      - [(7 + 3*sum(N_segments)) x 1] Global state vector
%    params - Parameter structure with params.brown_pili initialised
%
%  Outputs:
%    a         - [n_modes x 1] Mode amplitudes (dimensionless)
%    energy_kT - [n_modes x 1] Elastic energy per mode in units of kT
%                (equilibrium mean 1/2)
%    elastic_kT - [2 x 1] Total elastic energy of all pili in units of kT:
%                bending; twist.  The thermal modes are part of it; the
%                rest is held by the deterministic modes.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11

bp = params.brown_pili;
x  = pili_bending_coords(Z, params, bp.layout);
a  = bp.M_meas' * x;

if nargout >= 2
    kB = 1.380649e-23;
    kT_dimless = kB * params.dim.Temper / (params.scales.F_0 * params.scales.L_0);
    energy_kT  = 0.5 * a.^2 / kT_dimless;
end
if nargout >= 3
    e_joint    = 0.5 * (bp.layout.wK .* bp.layout.wE)' .* x.^2;
    is_twist   = (bp.layout.comp == 3)';
    elastic_kT = [sum(e_joint(~is_twist)); sum(e_joint(is_twist))] / kT_dimless;
end

end
