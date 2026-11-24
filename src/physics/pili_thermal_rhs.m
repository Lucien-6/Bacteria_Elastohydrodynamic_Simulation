function G = pili_thermal_rhs(F_body_frame, R_body, B_global, params, tau_axial, d3_lab_cell)
%PILI_THERMAL_RHS  Generalised force of the thermal loads carried by the pili
%
%  Description:
%    Pilus thermal-fluctuation scheme (v9.0.0+).  The thermal force on the
%    slow bending modes is realised as a force density on the pilus
%    centrelines together with the optional balancing axial torque density
%    from the constrained probe motion.  The centreline force is stored in
%    the body frame so it co-rotates with the cell during a Brownian
%    interval, while the optional axial torque density is a frame-invariant
%    material-frame scalar applied along each segment's current tangent d3.
%    Both are internal loads of the cell (the body carries the opposite net
%    force and torque; rows 1:6 are zeroed):
%
%        G = B_global(:, force cols) * f_lab
%          + B_global(:, twist cols) * (tau_axial .* d3_lab) ,   G(1:6) = 0 .
%
%    ode_coupled_rhs subtracts G from R_global (an applied generalised force
%    enters R with a minus sign).
%
%    Including the constrained probe's axial torque density cancels the
%    viscous spin drag of a bent filament without exciting twist: it
%    eliminates the ~3 % instantaneous modal response error and ~140 %
%    transient twist-rate leak on thermally bent filaments, reducing twist
%    energy by two orders of magnitude and restoring clean small-dt_B
%    convergence.
%
%  Syntax:
%    G = pili_thermal_rhs(F_body_frame, R_body, B_global, params)
%    G = pili_thermal_rhs(F_body_frame, R_body, B_global, params, tau_axial, d3_lab_cell)
%
%  Inputs:
%    F_body_frame - [3 x N_nodes] Dimensionless force density at the pilus
%                   nodes (all pili, in the order of X_all) in the body frame
%    R_body       - [3 x 3] Body rotation matrix (body frame -> lab frame)
%    B_global     - [n_state x n_force] Operator from build_coupled_lin_sys
%    params       - Parameter structure (N_body, N_segments)
%    tau_axial    - (Optional) [N_seg_all x 1] Dimensionless scalar axial
%                   torque density on each pilus segment (material frame)
%    d3_lab_cell  - (Optional) {M_pili x 1} Cell of [N_k x 3] lab-frame
%                   segment tangents (geom.d3); required when tau_axial is
%                   non-empty
%
%  Outputs:
%    G - [n_state x 1] Generalised force; rows 1:6 are zero
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11
%  Updated: 2026-11-24 — optional balancing axial torque density tau_axial.

n_segs_all = sum(params.N_segments(:));
n_nodes    = n_segs_all + numel(params.N_segments);
if ~isequal(size(F_body_frame), [3, n_nodes])
    error('pili_thermal_rhs:BadSize', ...
        'Thermal pilus force must be 3 x %d (pilus nodes), got %d x %d.', ...
        n_nodes, size(F_body_frame, 1), size(F_body_frame, 2));
end

f_lab = R_body * F_body_frame;
cols  = 3 * params.N_body + (1 : 3 * n_nodes);
G     = B_global(:, cols) * f_lab(:);

if nargin >= 5 && ~isempty(tau_axial)
    tau_vec = tau_axial(:);
    if numel(tau_vec) ~= n_segs_all
        error('pili_thermal_rhs:BadTauSize', ...
            'Thermal pilus axial torque must have %d entries (segments), got %d.', ...
            n_segs_all, numel(tau_vec));
    end
    if nargin < 6 || isempty(d3_lab_cell)
        error('pili_thermal_rhs:MissingTangents', ...
            'd3_lab_cell (geom.d3) is required when tau_axial is non-empty.');
    end
    d3_all = zeros(n_segs_all, 3);
    seg_off = 0;
    for k = 1:numel(params.N_segments)
        N_k = params.N_segments(k);
        d3_all(seg_off + (1:N_k), :) = d3_lab_cell{k};
        seg_off = seg_off + N_k;
    end
    tau_lab = (d3_all .* tau_vec)';                   % [3 x n_segs_all]
    t_cols  = 3 * (params.N_body + n_nodes) + (1 : 3 * n_segs_all);
    G       = G + B_global(:, t_cols) * tau_lab(:);
end

G(1:6) = 0;

end
