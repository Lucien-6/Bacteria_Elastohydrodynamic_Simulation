function x_dot = pili_bending_rates(Z, params, layout, z_dot)
%PILI_BENDING_RATES  Rate of the pilus curvature coordinates for given state rates
%
%  Description:
%    Evaluates d/dt of the curvature coordinates of pili_bending_coords for
%    one or several generalised-rate vectors z_dot of the coupled system,
%
%        x_dot(:, j) = (d x / d angles) * z_dot(7:end, j) ,
%
%    by a central difference of pili_bending_coords along each direction.
%    The curvature depends only on the Euler angles of the segments, so the
%    body rows z_dot(1:6, :) do not enter.  The difference acts on the same
%    discrete curvature operator that the elastic moments use
%    (compute_pili_curvature), which keeps the kinematic map and the elastic
%    force exactly conjugate; the truncation error is O(h^2) with a step of
%    1e-6 rad on the largest angle rate.
%
%  Syntax:
%    x_dot = pili_bending_rates(Z, params, layout, z_dot)
%
%  Inputs:
%    Z      - [(7 + 3*sum(N_segments)) x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%    layout - Layout struct from pili_bending_coords
%    z_dot  - [(6 + 3*sum(N_segments)) x n] Generalised rates
%             [V; Omega; dtheta; dphi; dpsi] (one column per case)
%
%  Outputs:
%    x_dot  - [nx x n] Rates of the curvature coordinates, in the time unit
%             of z_dot
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11
%  Updated: 2026-11-24 — 3-colour sparse curvature Jacobian J_x when n_col > 9
%                        (exact linear map, reduces O(n_col) calls to O(1)).

h_angle = 1e-6;                       % largest angle perturbation [rad]

n_col = size(z_dot, 2);
x_dot = zeros(layout.nx, n_col);
if n_col == 0 || layout.nx == 0
    return;
end
Z = Z(:);

% Fast path for multi-column probes (e.g. compute_effective_coupled_brownian):
%   x_dot = (dx / d_angles) * z_dot(7:end, :)
% Each segment j enters ONLY interior nodes j (if j>=2) and j+1 (if j<=N_k-1),
% so stride-3 coloured central differences build the exact block-banded
% Jacobian J_x in a fixed 18 curvature evaluations per pilus regardless of n_col.
if n_col > 9
    N_segs = params.N_segments(:)';
    M      = params.M_pili;
    [~, ~, pili_angles] = unpack_state(Z, N_segs);

    row_pos = 0;
    col_pos = 0;
    for k = 1:M
        N_k = N_segs(k);
        if N_k < 2
            col_pos = col_pos + 3 * N_k;
            continue;
        end
        n_rows_k = 3 * (N_k - 1);
        n_cols_k = 3 * N_k;
        J_k = zeros(n_rows_k, n_cols_k);

        th0 = pili_angles{k}.theta(:);
        ph0 = pili_angles{k}.phi(:);
        ps0 = pili_angles{k}.psi(:);

        for comp = 1:3
            for color = 0:2
                js = (1 + color) : 3 : N_k;
                if isempty(js)
                    continue;
                end
                th_p = th0;  ph_p = ph0;  ps_p = ps0;
                th_m = th0;  ph_m = ph0;  ps_m = ps0;
                if comp == 1
                    th_p(js) = th_p(js) + h_angle;
                    th_m(js) = th_m(js) - h_angle;
                    c_off = 0;
                elseif comp == 2
                    ph_p(js) = ph_p(js) + h_angle;
                    ph_m(js) = ph_m(js) - h_angle;
                    c_off = N_k;
                else
                    ps_p(js) = ps_p(js) + h_angle;
                    ps_m(js) = ps_m(js) - h_angle;
                    c_off = 2 * N_k;
                end

                [d1p, d2p, d3p] = pili_directors(th_p, ph_p, ps_p);
                kp = compute_pili_curvature(d1p, d2p, d3p, N_k);
                [d1m, d2m, d3m] = pili_directors(th_m, ph_m, ps_m);
                km = compute_pili_curvature(d1m, d2m, d3m, N_k);

                dk = (kp(2:N_k, :) - km(2:N_k, :))' / (2 * h_angle);   % [3 x (N_k-1)]
                for ii = 1:numel(js)
                    jj = js(ii);
                    col_j = c_off + jj;
                    % Interior node indices in 1..(N_k-1) correspond to nodes 2..N_k:
                    %   node jj     -> interior index jj - 1 (when jj >= 2)
                    %   node jj + 1 -> interior index jj     (when jj <= N_k - 1)
                    if jj >= 2
                        r_loc = 3 * (jj - 2) + (1:3);
                        J_k(r_loc, col_j) = dk(:, jj - 1);
                    end
                    if jj <= N_k - 1
                        r_loc = 3 * (jj - 1) + (1:3);
                        J_k(r_loc, col_j) = dk(:, jj);
                    end
                end
            end
        end

        x_dot(row_pos + (1:n_rows_k), :) = ...
            J_k * z_dot(6 + col_pos + (1:n_cols_k), :);
        row_pos = row_pos + n_rows_k;
        col_pos = col_pos + n_cols_k;
    end
    return;
end

for j = 1:n_col
    v = z_dot(7:end, j);
    s = max(abs(v));
    if ~(s > 0)
        continue;                     % no angle rate: curvature is unchanged
    end
    step = h_angle / s;
    Z_p = Z;   Z_p(8:end) = Z_p(8:end) + step * v;
    Z_m = Z;   Z_m(8:end) = Z_m(8:end) - step * v;
    x_dot(:, j) = (pili_bending_coords(Z_p, params, layout) - ...
                   pili_bending_coords(Z_m, params, layout)) / (2 * step);
end

end
