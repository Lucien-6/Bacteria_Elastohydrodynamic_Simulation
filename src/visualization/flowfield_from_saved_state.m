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

geom    = reconstruct_geometry(Z, params);
X_cm    = geom.X_cm;
X_all   = geom.X_all;
R_body  = geom.R_body;
N_total = size(X_all, 1);

logf('[FLOW] Computing ODE right-hand side for force extraction...\n');
[dZ, ode_hydro] = ode_coupled_rhs(t_dimless, Z, params);

V    = dZ(1:3);
dqdt = dZ(4:7);
q    = geom.q;

q_conj     = [q(1); -q(2:4)];
omega_quat = 2 * quat_multiply(dqdt, q_conj);
Omega      = omega_quat(2:4);

offset      = 7;
pili_dtheta = cell(M_pili, 1);
pili_dphi   = cell(M_pili, 1);
pili_dpsi   = cell(M_pili, 1);
for k = 1:M_pili
    N_k = N_segs(k);
    pili_dtheta{k} = dZ(offset + 1         : offset + N_k);
    pili_dphi{k}   = dZ(offset + N_k + 1   : offset + 2 * N_k);
    pili_dpsi{k}   = dZ(offset + 2 * N_k + 1 : offset + 3 * N_k);
    offset = offset + 3 * N_k;
end

v_trans = zeros(3 * N_total, 1);
for i = 1:N_body
    r_i = X_all(i, :)' - X_cm;
    v_i = V + cross(Omega, r_i);
    v_trans(3 * (i - 1) + 1 : 3 * i) = v_i;
end

pt_idx = N_body;
for k = 1:M_pili
    N_k   = N_segs(k);
    rot_k = params.rot{k};
    ds    = params.ell(k) / N_k;

    theta_k  = geom.pili_angles{k}.theta;
    phi_k    = geom.pili_angles{k}.phi;
    dtheta_k = pili_dtheta{k};
    dphi_k   = pili_dphi{k};

    v_prev = V + cross(Omega, X_all(pt_idx + 1, :)' - X_cm);
    v_trans(3 * pt_idx + 1 : 3 * pt_idx + 3) = v_prev;

    for j = 1:N_k
        sT = sin(theta_k(j));  cT = cos(theta_k(j));
        sP = sin(phi_k(j));    cP = cos(phi_k(j));

        de3_basis = [cT * cP * dtheta_k(j) - sT * sP * dphi_k(j);
                     cT * sP * dtheta_k(j) + sT * cP * dphi_k(j);
                     -sT * dtheta_k(j)];
        v_next = v_prev + ds * (rot_k' * de3_basis);
        v_trans(3 * (pt_idx + j) + 1 : 3 * (pt_idx + j) + 3) = v_next;
        v_prev = v_next;
    end
    pt_idx = pt_idx + N_k + 1;
end

v_twist = zeros(sum(N_segs), 1);
tw_idx  = 0;
for k = 1:M_pili
    theta_k = geom.pili_angles{k}.theta;
    for j = 1:N_segs(k)
        tw_idx = tw_idx + 1;
        v_twist(tw_idx) = pili_dpsi{k}(j) + cos(theta_k(j)) * pili_dphi{k}(j);
    end
end

has_bg_flow = isfield(params, 'bg_flow') && isfield(params.bg_flow, 'type') ...
    && ~strcmp(params.bg_flow.type, 'off');
if has_bg_flow
    u_bg_all = coupled_background_flow(X_all, geom, params);
    v_trans  = v_trans - u_bg_all(1:3 * N_total);
    v_twist  = v_twist - u_bg_all(3 * N_total + 1:end);
    bg_fun   = @(X) background_flow_velocity(X, params);
    logf(['[FLOW] Background flow ''%s'' active: forces solved from ', ...
          'relative velocity; ambient flow added to plotted field.\n'], ...
        params.bg_flow.type);
else
    bg_fun = [];
end

logf('[FLOW] Solving for force coefficients (%d DOFs)...\n', 3 * N_total);
M_weighted = ode_hydro.M_weighted;
M_twist = ode_hydro.M_twist;
D_twist = ode_hydro.D_twist;
C_twist = ode_hydro.C_twist;
C_twist_reverse = ode_hydro.C_twist_reverse;
if ~isempty(M_twist)
    Minv_v   = M_weighted \ v_trans;
    Minv_C   = M_weighted \ C_twist;
    check_linear_solve(M_weighted, [Minv_v, Minv_C], ...
        [v_trans, C_twist], 'flow-field translational solve');
    S_twist  = M_twist - C_twist_reverse * Minv_C;
    m_twist  = S_twist \ (v_twist - C_twist_reverse * Minv_v);
    check_linear_solve(S_twist, m_twist, ...
        v_twist - C_twist_reverse * Minv_v, ...
        'flow-field twist Schur solve');
    f_coeffs = Minv_v - Minv_C * m_twist;
    tau      = D_twist * m_twist;
    clear Minv_v Minv_C
else
    tau = zeros(0, 1);
    f_coeffs = M_weighted \ v_trans;
    check_linear_solve(M_weighted, f_coeffs, v_trans, ...
        'flow-field translational solve');
end
clear ode_hydro M_weighted M_twist D_twist C_twist C_twist_reverse

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

    s_count     = 0;
    pili_offset = N_body;
    for k = 1:M_pili
        N_k = N_segs(k);
        for s = 1:N_k
            s_count = s_count + 1;
            g0 = pili_offset + s;
            g1 = pili_offset + s + 1;
            seg_data.y0(s_count, :) = X_all(g0, :);
            seg_data.y1(s_count, :) = X_all(g1, :);
            seg_data.f0(s_count, :) = f_coeffs(3 * (g0 - 1) + 1 : 3 * g0)';
            seg_data.f1(s_count, :) = f_coeffs(3 * (g1 - 1) + 1 : 3 * g1)';
        end
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
