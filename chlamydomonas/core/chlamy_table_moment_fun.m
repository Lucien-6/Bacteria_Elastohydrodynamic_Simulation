function fun = chlamy_table_moment_fun(tab, k, opt)
%CHLAMY_TABLE_MOMENT_FUN  Custom active-moment handle from a tabulated q(s,t)
%
%  Description:
%    Returns the function handle required by active_moment.type = 'custom'
%    for flagellum K from the table produced by chlamy_inverse_moment /
%    chlamy_fit_moment (and optionally extended by chlamy_make_perp_envelope):
%
%        q_par,j(t)  = c_0j + amp * sum_n [ a_nj cos(n th) + b_nj sin(n th) ],
%        q_perp,j(t) = amp_perp * sum_n [ a_nj^p cos(n th) + b_nj^p sin(n th) ],
%        th          = w_star * freq_ratio * t + phase.
%
%    In-plane moment acts along the body-fixed beat-plane normal n(t)
%    (default [0; 1; 0] in the clamped root frame).  When tab.coef_perp{k}
%    is present, the orthogonal bending moment acts along
%    b_hat(j) = normalize(d3N(j) x n(t)) with zero parasitic twist.
%
%  Syntax:
%    fun = chlamy_table_moment_fun(tab, k)
%    fun = chlamy_table_moment_fun(tab, k, opt)
%
%  Inputs:
%    tab - struct (.coef, .w_star, .nq, .N_seg; optional .coef_perp)
%    k   - Flagellum index (1 or 2)
%    opt - (Optional) struct: .normal_local [3 x 1] (default [0; 1; 0]),
%          .amp (1), .amp_perp (1), .phase (0), .freq_ratio (1)
%
%  Output:
%    fun - Handle fun(t, k, N_k, d1, d2, d3, d1N, d2N, d3N) -> [3 x N_k]
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-09
%  Updated: 2026-11-24  input validation, moved to chlamydomonas/core/

arguments
    tab (1,1) struct
    k   (1,1) double {mustBePositive, mustBeInteger}
    opt (1,1) struct = struct()
end

assert(isfield(tab, 'coef') && iscell(tab.coef) && k <= numel(tab.coef), ...
    'chlamy_table_moment_fun:InvalidTable', 'tab.coef{%d} is missing.', k);

n_loc = getopt(opt, 'normal_local', [0; 1; 0]);
n_loc = n_loc(:) / norm(n_loc);

par = struct('cq', tab.coef{k}, 'w_star', tab.w_star, 'nq', tab.nq, ...
    'n_loc', n_loc, 'amp', getopt(opt, 'amp', 1), ...
    'phase', getopt(opt, 'phase', 0), ...
    'freq_ratio', getopt(opt, 'freq_ratio', 1));
if isfield(tab, 'coef_perp') && iscell(tab.coef_perp) && ...
        numel(tab.coef_perp) >= k && ~isempty(tab.coef_perp{k})
    par.cq_perp  = tab.coef_perp{k};
    par.amp_perp = getopt(opt, 'amp_perp', 1);
end
fun = @(t, ~, N_k, d1, d2, d3, d1N, d2N, d3N) ...
    chlamy_table_moment_eval(par, t, N_k, d1, d2, d3, d1N, d2N, d3N);

    function M = moment_fun(t, ~, N_k, d1, d2, d3, d1N, d2N, d3N) %#ok<DEFNU>
        % LEGACY PATH - pre-2026-11-11 saved nested handles (planar only)
        par_legacy = struct('cq', tab.coef{k}, 'w_star', tab.w_star, ...
            'nq', tab.nq);
        M = chlamy_table_moment_eval(par_legacy, t, N_k, d1, d2, d3, ...
            d1N, d2N, d3N);
    end
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
