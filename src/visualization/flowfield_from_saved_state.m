function hydro = flowfield_from_saved_state(Z, params, t_dimless, logf)
%FLOWFIELD_FROM_SAVED_STATE  Recover Stokeslet and rotlet data from one state
%
%  Description:
%    Reconstructs the coupled geometry at one dimensionless time, evaluates
%    the existing ODE right-hand side, and solves the same regularised
%    Stokeslet / rotlet system used by the final-state flow field.  The
%    returned structure is the input of flowfield_grid_velocity.
%
%  Syntax:
%    hydro = flowfield_from_saved_state(Z, params, t_dimless, logf)
%
%  Inputs:
%    Z          - Complete dimensionless state vector
%    params     - Parameter structure synchronised with Z, including rot
%    t_dimless  - Dimensionless time passed to ode_coupled_rhs
%    logf       - Logging function handle
%
%  Outputs:
%    hydro - Structure with geometry, kernel data, background-flow handle,
%            and unit-conversion factors
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28

if nargin < 4 || isempty(logf)
    logf = @(varargin) [];
end

if ~iscell(params.rot) || numel(params.rot) ~= params.M_pili
    error('flowfield_from_saved_state:RotShapeMismatch', ...
        ['params.rot must be a cell array of length M_pili = %d, got %d. ', ...
         'Z and params.rot must use the same pilus basis.'], ...
        params.M_pili, numel(params.rot));
end
expected_Z_len = 7 + 3 * sum(params.N_segments);
if numel(Z) ~= expected_Z_len
    error('flowfield_from_saved_state:ZShapeMismatch', ...
        'numel(Z) = %d does not match 7 + 3*sum(N_segs) = %d.', ...
        numel(Z), expected_Z_len);
end
Z = Z(:);

L_0     = params.scales.L_0;
U_0     = params.scales.U_0;
L_to_um = L_0 * 1e6;
um_to_L = 1 / L_to_um;

N_body = params.N_body;
M_pili = params.M_pili;
N_segs = params.N_segments;

logf('[FLOW] Computing ODE right-hand side for force extraction...\n');
[~, ode_hydro] = ode_coupled_rhs(t_dimless, Z, params);

geom     = ode_hydro.geom;
X_cm     = geom.X_cm;
X_all    = geom.X_all;
R_body   = geom.R_body;
f_coeffs = ode_hydro.f_coeffs;
tau      = ode_hydro.tau;

has_bg_flow = isfield(params, 'bg_flow') && isfield(params.bg_flow, 'type') ...
    && ~strcmp(params.bg_flow.type, 'off');
if has_bg_flow
    bg_fun = @(X) background_flow_velocity(X, params);
    logf(['[FLOW] Background flow ''%s'' active: forces solved from ', ...
          'relative velocity; ambient flow added to plotted field.\n'], ...
        params.bg_flow.type);
else
    bg_fun = [];
end

body_data.X  = X_all(1:N_body, :);
body_data.fx = f_coeffs(1:3:3 * N_body);
body_data.fy = f_coeffs(2:3:3 * N_body);
body_data.fz = f_coeffs(3:3:3 * N_body);
body_data.e2 = params.eps_body^2;

seg_data = struct('active', false);
if M_pili > 0
    N_seg_all = sum(N_segs);
    seg_data.y0 = zeros(N_seg_all, 3);
    seg_data.y1 = zeros(N_seg_all, 3);
    seg_data.f0 = zeros(N_seg_all, 3);
    seg_data.f1 = zeros(N_seg_all, 3);

    f_nodes     = reshape(f_coeffs, 3, []).';
    s_count     = 0;
    pili_offset = N_body;
    for k = 1:M_pili
        N_k = N_segs(k);
        s_idx = s_count + (1:N_k);
        g0 = pili_offset + (1:N_k);
        g1 = g0 + 1;
        seg_data.y0(s_idx, :) = X_all(g0, :);
        seg_data.y1(s_idx, :) = X_all(g1, :);
        seg_data.f0(s_idx, :) = f_nodes(g0, :);
        seg_data.f1(s_idx, :) = f_nodes(g1, :);
        s_count     = s_count + N_k;
        pili_offset = pili_offset + N_k + 1;
    end

    seg_data.v  = seg_data.y0 - seg_data.y1;
    seg_data.L  = sqrt(sum(seg_data.v.^2, 2));
    seg_data.L2 = seg_data.L.^2;
    seg_data.e2 = params.epsilon^2;
    seg_data.has_rotlet = ~isempty(tau);
    if seg_data.has_rotlet
        seg_data.tau_x = tau(1:3:end);
        seg_data.tau_y = tau(2:3:end);
        seg_data.tau_z = tau(3:3:end);
    end
    seg_data.active = true;
end

hydro.geom        = geom;
hydro.body_data   = body_data;
hydro.seg_data    = seg_data;
hydro.bg_fun      = bg_fun;
hydro.X_cm        = X_cm;
hydro.R_body      = R_body;
hydro.L_0         = L_0;
hydro.U_0         = U_0;
hydro.L_to_um     = L_to_um;
hydro.um_to_L     = um_to_L;
hydro.um_s_factor = U_0 * 1e6;
hydro.f_coeffs    = f_coeffs;
hydro.tau         = tau;

end
