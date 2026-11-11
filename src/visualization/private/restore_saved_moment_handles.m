function params = restore_saved_moment_handles(params)
%RESTORE_SAVED_MOMENT_HANDLES  Make saved custom active-moment handles callable
%
%  Description:
%    A custom active moment (params.am_cfg{k}.fun) read back from
%    simulation_results.mat must be callable in a fresh session.  Two
%    things can prevent that:
%      1. The handle refers to a function in the chlamydomonas folder, which
%         is not on the path of a post-processing session.  The folder is
%         added to the path here.
%      2. Results written by earlier versions store a handle to a NESTED
%         function (chlamy_table_moment_fun/moment_fun or
%         chlamy_moment_fun/moment_fun).  MATLAB saves only a snapshot of the
%         parent workspace for such a handle and runs the CURRENT source
%         text, so the handle fails ("undefined variable") once the source
%         file has been edited.  Such handles are rebuilt here from the
%         saved workspace values with the stand-alone evaluators
%         chlamy_table_moment_eval / chlamy_moment_eval.  Optional values
%         missing from an old snapshot (n_loc, amp, phase, freq_ratio) take the
%         symmetric-drive defaults, which is what those versions computed.
%
%    Known workspace shapes for chlamy_table_moment_fun/moment_fun:
%      - Earliest saves close over the whole table: fields tab, k
%      - Intermediate saves close over extracted coefficients: cq, w_star, nq
%        (and optionally n_loc, amp, phase, fratio / freq_ratio, opt)
%      - Current saves use an anonymous handle that already stores par by
%        value; those are left unchanged after the path is set.
%
%  Syntax:
%    params = restore_saved_moment_handles(params)
%
%  Inputs:
%    params - Parameter structure loaded from simulation_results.mat
%
%  Outputs:
%    params - Same structure with callable am_cfg{k}.fun handles
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11
%  Updated: 2026-11-11 — rebuild from saved tab+k workspace of Fit-1p00A era

if ~isfield(params, 'am_cfg') || ~iscell(params.am_cfg)
    return
end

proj_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
chlamy_dir = fullfile(proj_root, 'chlamydomonas');
path_ready = false;

for k = 1:numel(params.am_cfg)
    am = params.am_cfg{k};
    if ~(isstruct(am) && isfield(am, 'fun') && isa(am.fun, 'function_handle'))
        continue
    end
    name = func2str(am.fun);
    if ~contains(name, 'chlamy_')
        continue
    end
    if ~path_ready && isfolder(chlamy_dir)
        addpath(chlamy_dir);
        path_ready = true;
    end
    info = functions(am.fun);
    if ~strcmp(info.type, 'nested')
        continue
    end
    ws = struct();
    if isfield(info, 'workspace') && ~isempty(info.workspace)
        ws = info.workspace{1};
    end
    switch name
        case 'chlamy_table_moment_fun/moment_fun'
            par = table_par(ws);
            params.am_cfg{k}.fun = @(t, ~, N_k, d1, d2, d3, d1N, d2N, d3N) ...
                chlamy_table_moment_eval(par, t, N_k, d1, d2, d3, d1N, d2N, d3N);
        case 'chlamy_moment_fun/moment_fun'
            par = wave_par(ws);
            params.am_cfg{k}.fun = @(t, ~, N_k, ~, ~, ~, ~, ~, ~) ...
                chlamy_moment_eval(par, t, N_k);
        otherwise
            continue
    end
    fprintf(['[FLOW] Rebuilt the saved nested active-moment handle of ', ...
        'flagellum %d (%s) as a stand-alone evaluator.\n'], k, name);
end
end


function par = table_par(ws)
%TABLE_PAR  Recover the stand-alone parameter struct from a saved workspace
if isfield(ws, 'par') && isstruct(ws.par) && isfield(ws.par, 'cq')
    par = ws.par;
elseif all(isfield(ws, {'cq', 'w_star', 'nq'}))
    par = struct('cq', ws.cq, 'w_star', ws.w_star, 'nq', ws.nq);
elseif isfield(ws, 'tab') && isfield(ws, 'k')
    tab = ws.tab;
    kk = ws.k;
    if ~(isstruct(tab) && isfield(tab, 'coef') && iscell(tab.coef) && ...
            isnumeric(kk) && isscalar(kk) && kk == round(kk) && ...
            kk >= 1 && kk <= numel(tab.coef))
        error('restore_saved_moment_handles:BadTableHandle', ...
            ['The saved chlamy_table_moment_fun handle has tab/k, but ', ...
             'tab.coef{%g} is missing or invalid.'], double(kk));
    end
    check_fields(tab, {'w_star', 'nq'}, 'chlamy_table_moment_fun.tab');
    par = struct('cq', tab.coef{kk}, 'w_star', tab.w_star, 'nq', tab.nq);
else
    found = '(empty)';
    if ~isempty(fieldnames(ws))
        found = strjoin(fieldnames(ws), ', ');
    end
    error('restore_saved_moment_handles:IncompleteHandle', ...
        ['The saved chlamy_table_moment_fun handle has neither tab+k nor ', ...
         'cq/w_star/nq.  Workspace fields: %s.'], found);
end

% Optional drive knobs may appear as parent locals or inside opt.
if isfield(ws, 'n_loc'), par.n_loc = ws.n_loc; end
if isfield(ws, 'amp'),   par.amp   = ws.amp;   end
if isfield(ws, 'phase'), par.phase = ws.phase; end
if isfield(ws, 'freq_ratio')
    par.freq_ratio = ws.freq_ratio;
elseif isfield(ws, 'fratio')
    par.freq_ratio = ws.fratio;
end
if isfield(ws, 'opt') && isstruct(ws.opt)
    opt = ws.opt;
    if ~isfield(par, 'n_loc') && isfield(opt, 'normal_local')
        par.n_loc = opt.normal_local;
    end
    if ~isfield(par, 'amp') && isfield(opt, 'amp')
        par.amp = opt.amp;
    end
    if ~isfield(par, 'phase') && isfield(opt, 'phase')
        par.phase = opt.phase;
    end
    if ~isfield(par, 'freq_ratio') && isfield(opt, 'freq_ratio')
        par.freq_ratio = opt.freq_ratio;
    end
end
end


function par = wave_par(ws)
if isfield(ws, 'par') && isstruct(ws.par) && isfield(ws.par, 'q_star')
    par = ws.par;
    return
end
need = {'q_star', 'w_star', 'kw', 'phi', 's_taper'};
check_fields(ws, need, 'chlamy_moment_fun');
par = struct('q_star', ws.q_star, 'w_star', ws.w_star, 'kw', ws.kw, ...
    'phi', ws.phi, 's_taper', ws.s_taper, 'standing', false, 'root_env', []);
if isfield(ws, 'standing')
    par.standing = logical(ws.standing);
end
if isfield(ws, 'root_env')
    par.root_env = ws.root_env;
end
end


function check_fields(ws, need, owner)
missing = need(~isfield(ws, need));
if ~isempty(missing)
    found = '(empty)';
    if ~isempty(fieldnames(ws))
        found = strjoin(fieldnames(ws), ', ');
    end
    error('restore_saved_moment_handles:IncompleteHandle', ...
        ['The saved %s handle has no stored value for: %s. ', ...
         'Workspace fields: %s.'], ...
        owner, strjoin(missing, ', '), found);
end
end
