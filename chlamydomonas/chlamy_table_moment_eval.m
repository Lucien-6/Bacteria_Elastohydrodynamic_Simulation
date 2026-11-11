function M = chlamy_table_moment_eval(par, t, N_k, d1, d2, d3, d1N, d2N, d3N)
%CHLAMY_TABLE_MOMENT_EVAL  Evaluate the tabulated Chlamydomonas active moment
%
%  Description:
%    Stand-alone evaluator behind chlamy_table_moment_fun.  All data of the
%    drive travel in the struct PAR, which is stored by value inside the
%    anonymous handle returned by chlamy_table_moment_fun.  A saved
%    simulation_results.mat therefore keeps working after this file is
%    edited, as long as the fields of PAR keep their meaning (missing optional
%    fields fall back to the symmetric drive).
%
%  Syntax:
%    M = chlamy_table_moment_eval(par, t, N_k, d1, d2, d3, d1N, d2N, d3N)
%
%  Inputs:
%    par - struct with .cq [N_seg x (2 nq + 1)] coefficients, .w_star,
%          .nq, and optional .n_loc (default [0; 1; 0]), .amp (1),
%          .phase (0), .freq_ratio (1)
%    t   - Dimensionless time (t_phys / T_0)
%    N_k - Number of balance nodes of the flagellum
%    d1, d2, d3, d1N, d2N, d3N - Root frame and node frames, as passed by
%          coupled_active_moment
%
%  Output:
%    M - [3 x N_k] internal active moment
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11

cq     = par.cq;
w_star = par.w_star;
nq     = par.nq;
n_loc  = getopt(par, 'n_loc', [0; 1; 0]);
n_loc  = n_loc(:) / norm(n_loc);
amp    = getopt(par, 'amp', 1);
phase  = getopt(par, 'phase', 0);
fratio = getopt(par, 'freq_ratio', 1);

assert(N_k == size(cq, 1), 'chlamy_table_moment_fun:Size', ...
    'Table has %d nodes, flagellum has %d.', size(cq, 1), N_k);
th = w_star * fratio * t + phase;
b = zeros(2*nq + 1, 1);
b(1) = 1;
for n = 1:nq
    b(2*n)   = cos(n * th);
    b(2*n+1) = sin(n * th);
end
q = cq(:, 1) + amp * (cq(:, 2:end) * b(2:end));
if ~isempty(d1) && ~isequal(n_loc, [0; 1; 0])
    % beat-plane normal in the lab/basis frame, from the clamped root frame
    n_vec = n_loc(1) * d1(1, :) + n_loc(2) * d2(1, :) + n_loc(3) * d3(1, :);   % [1 x 3]
    comp = [d1N(1:N_k, :) * n_vec.'; d2N(1:N_k, :) * n_vec.'; d3N(1:N_k, :) * n_vec.'];
    comp = reshape(comp, N_k, 3).';                                              % [3 x N_k]
    M = comp .* q(:).';
else
    M = [zeros(1, N_k); q(:).'; zeros(1, N_k)];
end
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
