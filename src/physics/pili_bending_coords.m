function [x, layout] = pili_bending_coords(Z, params, layout)
%PILI_BENDING_COORDS  Curvature coordinates of all pili at the interior joints
%
%  Description:
%    Returns the deviation of the material curvature / twist from its
%    intrinsic value at every interior joint (nodes 2..N_k) of every pilus,
%    stacked in the row order of the pilus moment-balance rows of the coupled
%    linear system:
%
%        x = [ kappa_1(2); kappa_2(2); kappa_3(2); kappa_1(3); ... ]  (pilus 1)
%            [ ...                                                 ]  (pilus 2)
%
%    with kappa in units of 1/L_pili (compute_pili_curvature) and the
%    intrinsic curvature params.kb0{k} subtracted.  These are the generalised
%    coordinates in which the elastic energy of the discrete Kirchhoff rod is
%    exactly quadratic,
%
%        E / (F_0 L_0) = (1/2) * sum_r  wK(r) * wE(r) * x(r)^2 ,
%
%    wK = [1, 1, 1/(1+sigma)] (bending, bending, twist) and
%    wE = 1 / (N_k * ell_k).  The moment-balance right-hand side assembled by
%    ode_coupled_rhs is R(rows) = wK .* x, i.e. minus the elastic generalised
%    force divided by wE.
%
%    Node 1 (root) is excluded: its three rows are the root boundary
%    condition, not a joint.  Node N_k+1 (tip) carries no moment.
%
%  Syntax:
%    [x, layout] = pili_bending_coords(Z, params)
%    x           = pili_bending_coords(Z, params, layout)
%
%  Inputs:
%    Z      - [(7 + 3*sum(N_segments)) x 1] Global state vector
%    params - Parameter structure (from init_coupled_params)
%    layout - (Optional) layout struct returned by a previous call; pass it
%             to skip rebuilding the index tables
%
%  Outputs:
%    x      - [nx x 1] Curvature coordinates, nx = 3 * sum(N_segments - 1)
%    layout - Struct with [1 x nx] index tables:
%               .rows   row of lin_sys / R_global of each coordinate
%               .pilus  pilus index
%               .node   node index (2..N_k)
%               .comp   1, 2 (bending about d1, d2) or 3 (twist)
%               .wK     stiffness weight (1, 1, 1/(1+sigma))
%               .wE     energy weight 1 / (N_k * ell_k)
%             and .nx, .n_state.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11

N_segs = params.N_segments(:)';
M      = params.M_pili;

if nargin < 3 || isempty(layout)
    layout = build_layout(params);
end

[~, ~, pili_angles] = unpack_state(Z, N_segs);

x   = zeros(layout.nx, 1);
pos = 0;
for k = 1:M
    N_k = N_segs(k);
    if N_k < 2
        continue;
    end
    [d1, d2, d3] = pili_directors(pili_angles{k}.theta, pili_angles{k}.phi, ...
                                  pili_angles{k}.psi);
    kappa = compute_pili_curvature(d1, d2, d3, N_k) - params.kb0{k};   % [(N_k+1) x 3]
    blk   = kappa(2:N_k, :)';                                          % [3 x (N_k-1)]
    n_blk = 3 * (N_k - 1);
    x(pos + (1:n_blk)) = blk(:);
    pos = pos + n_blk;
end

end


%% ==================== Local Helper Function ====================

function layout = build_layout(params)
%BUILD_LAYOUT  Index tables of the interior-joint curvature coordinates
N_segs = params.N_segments(:)';
M      = params.M_pili;
nx     = 3 * sum(max(N_segs - 1, 0));

rows  = zeros(1, nx);   pilus = zeros(1, nx);   node = zeros(1, nx);
comp  = zeros(1, nx);   wK    = zeros(1, nx);   wE   = zeros(1, nx);

pos        = 0;
row_offset = 6;
for k = 1:M
    N_k = N_segs(k);
    for i = 2:N_k
        idx        = pos + (1:3);
        rows(idx)  = row_offset + 3*(i-1) + (1:3);
        pilus(idx) = k;
        node(idx)  = i;
        comp(idx)  = 1:3;
        wK(idx)    = [1, 1, 1 / (1 + params.sigma)];
        wE(idx)    = 1 / (N_k * params.ell(k));
        pos        = pos + 3;
    end
    row_offset = row_offset + 3 * N_k;
end

layout = struct('rows', rows, 'pilus', pilus, 'node', node, 'comp', comp, ...
                'wK', wK, 'wE', wE, 'nx', nx, 'n_state', 6 + 3 * sum(N_segs));
end
