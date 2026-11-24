function M = chlamy_table_moment_eval(par, t, N_k, d1, d2, d3, d1N, d2N, d3N)
%CHLAMY_TABLE_MOMENT_EVAL  Evaluate the tabulated Chlamydomonas active moment
%
%  Description:
%    Stand-alone evaluator behind chlamy_table_moment_fun.  All drive data
%    travel by value in PAR inside the anonymous handle so saved results
%    remain callable across code updates.
%
%    In-plane component (along body-fixed tilted beat-plane normal n_hat):
%        q_par(t)  = c0 + amp * sum_n [a_n cos(n th) + b_n sin(n th)]
%
%    Optional 3-D out-of-plane component (par.cq_perp, static c0_p == 0):
%        q_perp(t) = amp_perp * (c0_p + sum_n [a_n^p cos(n th) + b_n^p sin(n th)])
%        b_hat(j)  = normalize(d3N(j) x n_hat)
%        M(j)      = q_par(j) * n_hat + q_perp(j) * b_hat(j)
%
%    Because b_hat(j) is strictly perpendicular to the local tangent d3N(j),
%    the out-of-plane drive injects pure out-of-plane bending (toward +n_hat
%    for q_perp > 0) with zero parasitic twist moment (M(3,j) == 0).  In the
%    degenerate limit |d3N x n_hat| < 1e-12, dividing by max(bn, 1e-12)
%    smoothly attenuates b_hat toward zero without NaN.
%
%  Syntax:
%    M = chlamy_table_moment_eval(par, t, N_k, d1, d2, d3, d1N, d2N, d3N)
%
%  Inputs:
%    par - struct with .cq [N_seg x (2 nq + 1)], .w_star, .nq, and optional
%          .n_loc (default [0; 1; 0]), .amp (1), .phase (0),
%          .freq_ratio (1), .cq_perp ([]), .amp_perp (1)
%    t   - Dimensionless time (t_phys / T_0)
%    N_k - Number of balance nodes of the flagellum
%    d1, d2, d3, d1N, d2N, d3N - Root frame and node frames from
%          coupled_active_moment
%
%  Output:
%    M - [3 x N_k] internal active moment in node director frame (d1N,d2N,d3N)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11
%  Updated: 2026-11-24  code review cleanup (direct 3-row projection, docs)

cq     = par.cq;
w_star = par.w_star;
nq     = par.nq;
n_loc  = getopt(par, 'n_loc', [0; 1; 0]);
n_loc  = n_loc(:) / norm(n_loc);
amp    = getopt(par, 'amp', 1);
phase  = getopt(par, 'phase', 0);
fratio = getopt(par, 'freq_ratio', 1);
cq_p   = getopt(par, 'cq_perp', []);
amp_p  = getopt(par, 'amp_perp', 1);

assert(N_k == size(cq, 1), 'chlamy_table_moment_eval:Size', ...
    'Table has %d nodes, flagellum has %d.', size(cq, 1), N_k);

th = w_star * fratio * t + phase;
b = zeros(2*nq + 1, 1);
b(1) = 1;
for n = 1:nq
    b(2*n)   = cos(n * th);
    b(2*n+1) = sin(n * th);
end
q = cq(:, 1) + amp * (cq(:, 2:end) * b(2:end));

has_perp = ~isempty(cq_p) && (amp_p ~= 0);
if has_perp
    assert(isequal(size(cq_p), size(cq)), 'chlamy_table_moment_eval:PerpSize', ...
        'cq_perp size must match cq size.');
    q_perp = amp_p * (cq_p(:, 1) + cq_p(:, 2:end) * b(2:end));
end

use_frame = has_perp || (~isempty(d1) && ~isequal(n_loc, [0; 1; 0]));
if ~use_frame
    M = [zeros(1, N_k); q(:).'; zeros(1, N_k)];
    return
end

assert(~isempty(d1) && ~isempty(d3N), 'chlamy_table_moment_eval:FrameRequired', ...
    'Root and node directors are required for tilted or 3-D active moments.');

% Beat-plane normal in the lab/basis frame, transported by the clamped root frame
n_vec = n_loc(1) * d1(1, :) + n_loc(2) * d2(1, :) + n_loc(3) * d3(1, :);   % [1 x 3]
nv = n_vec.';
comp_n = [(d1N(1:N_k, :) * nv).'; ...
          (d2N(1:N_k, :) * nv).'; ...
          (d3N(1:N_k, :) * nv).'];                                          % [3 x N_k]
M = comp_n .* q(:).';

if has_perp
    % Node-local out-of-plane bending axis b_hat(j) = normalize(d3N(j) x n_vec).
    % Orthogonal to d3N(j) everywhere => zero parasitic twist (row 3 == 0).
    tN = d3N(1:N_k, :);                                                     % [N_k x 3]
    bN = [tN(:, 2) * n_vec(3) - tN(:, 3) * n_vec(2), ...
          tN(:, 3) * n_vec(1) - tN(:, 1) * n_vec(3), ...
          tN(:, 1) * n_vec(2) - tN(:, 2) * n_vec(1)];
    bn = sqrt(sum(bN.^2, 2));
    bN = bN ./ max(bn, 1e-12);
    comp_b = [sum(d1N(1:N_k, :) .* bN, 2).'; ...
              sum(d2N(1:N_k, :) .* bN, 2).'; ...
              zeros(1, N_k)];                                               % [3 x N_k]
    M = M + comp_b .* q_perp(:).';
end
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
