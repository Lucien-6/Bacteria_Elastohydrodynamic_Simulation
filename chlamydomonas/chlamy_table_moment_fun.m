function fun = chlamy_table_moment_fun(tab, k, opt)
%CHLAMY_TABLE_MOMENT_FUN  Custom active-moment function from a tabulated q(s,t)
%
%  Description:
%    Returns the function handle required by active_moment.type = 'custom'
%    for flagellum K from the table produced by chlamy_inverse_moment (or
%    chlamy_fit_moment):
%
%        q_j(t) = c_0j + amp * sum_n [ a_nj cos(n theta) + b_nj sin(n theta) ],
%        theta  = w * freq_ratio * t + phase,
%
%    with j the balance node (1..N_seg) and w = 2 pi f the beat frequency.
%    The internal moment vector at node j is  q_j * n(t),  where n is the unit
%    normal of the beat plane.  n is body-fixed: it is stored in the frame of
%    the clamped root segment (NORMAL_LOCAL = D_local' * n0, columns d1, d2,
%    d3 of the root frame at t = 0) and followed through the root director
%    frame of the current state, so that the beat plane rotates with the cell.
%    The default NORMAL_LOCAL = [0; 1; 0] is d2 (in-plane bending of the
%    planar gait: M = [0; q; 0]).  The root column is zero (the root internal
%    moment is the clamped reaction).
%
%    Symmetry-breaking handles (all default to the symmetric drive):
%      amp        scales the OSCILLATING part only (the static moment c_0
%                 that holds the mean shape is left unchanged)
%      phase      beat phase shift [rad] of this flagellum
%      freq_ratio beat frequency of this flagellum relative to tab.f_Hz
%                 (a ratio different from 1 makes the drive non-periodic with
%                 respect to the other flagellum)
%
%  Syntax:
%    fun = chlamy_table_moment_fun(tab, k)
%    fun = chlamy_table_moment_fun(tab, k, opt)
%
%  Inputs:
%    tab - struct (.coef, .w_star, .nq, .N_seg) from chlamy_inverse_moment or
%          chlamy_fit_moment
%    k   - Flagellum index (1 or 2)
%    opt - (Optional) struct: .normal_local [3 x 1] (default [0; 1; 0]),
%          .amp (1), .phase (0), .freq_ratio (1)
%
%  Output:
%    fun - Handle fun(t, k, N_k, d1, d2, d3, d1N, d2N, d3N) -> [3 x N_k],
%          t = dimensionless time (t_phys / T_0)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-09
%  Updated: 2026-11-10 (beat-plane normal, amp / phase / frequency handles)
%  Updated: 2026-11-11 (anonymous handle + chlamy_table_moment_eval, so saved
%                       results survive later edits of this file; legacy nested
%                       moment_fun kept on tab/k so pre-2026-11-11 .mat files
%                       remain callable)

if nargin < 3 || isempty(opt), opt = struct(); end
n_loc  = getopt(opt, 'normal_local', [0; 1; 0]);
n_loc  = n_loc(:) / norm(n_loc);

% All drive data are stored by value in PAR.  The handle is an anonymous
% function that calls the stand-alone chlamy_table_moment_eval, NOT a nested
% function: a nested handle saved in simulation_results.mat keeps only a
% snapshot of this workspace and breaks as soon as this file is edited.
par = struct('cq', tab.coef{k}, 'w_star', tab.w_star, 'nq', tab.nq, ...
    'n_loc', n_loc, 'amp', getopt(opt, 'amp', 1), ...
    'phase', getopt(opt, 'phase', 0), ...
    'freq_ratio', getopt(opt, 'freq_ratio', 1));
fun = @(t, ~, N_k, d1, d2, d3, d1N, d2N, d3N) ...
    chlamy_table_moment_eval(par, t, N_k, d1, d2, d3, d1N, d2N, d3N);

    function M = moment_fun(t, ~, N_k, d1, d2, d3, d1N, d2N, d3N)
        % LEGACY PATH - keep the name and the tab/k closure.  Results saved
        % before 2026-11-11 store a handle to this nested function whose
        % workspace is {tab, k}.  Evaluating through tab.coef{k} with the
        % symmetric-drive defaults recovers that era's moment; newer runs
        % never call this path (they use the anonymous handle above).
        par_legacy = struct('cq', tab.coef{k}, 'w_star', tab.w_star, ...
            'nq', tab.nq);
        M = chlamy_table_moment_eval(par_legacy, t, N_k, d1, d2, d3, ...
            d1N, d2N, d3N);
    end
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
