function profile_short_run()
%PROFILE_SHORT_RUN  Profiler ON: short coupled run, then profsave + top-25 print.
%
%  Usage (from project root):
%    profile_short_run
%    matlab -batch "cd('F:\...\Coupled_Body_Pili'); profile_short_run"
%
%  Outputs under output/profile_runs/:
%    profsave_<timestamp>/   HTML profile (if profsave available)
%    profile_info_<timestamp>.mat   raw profile('info') struct
%
%  Author: Lucien (lucien-6@qq.com)
%  Date: 2026-04-21

root = fileparts(mfilename('fullpath'));
if isempty(root)
    root = pwd;
end
cd(root);

cfg = local_build_cfg();

profile clear
profile on -timer cpu

run_simulation(cfg);

profile off
p = profile('info');

outDir = fullfile(root, 'output', 'profile_runs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
stamp = datestr(now, 'yyyymmdd_HHMMSS');
prof_dir = fullfile(outDir, ['profsave_' stamp]);
mat_path = fullfile(outDir, ['profile_info_' stamp '.mat']);

try
    profsave(p, prof_dir);
    fprintf('\n[PROFILE] HTML report: %s\n', prof_dir);
catch ME
    fprintf('\n[PROFILE] profsave failed (%s) — saving .mat only.\n', ME.message);
    prof_dir = '';
end

save(mat_path, 'p', '-v7.3');
fprintf('[PROFILE] Raw struct saved: %s\n', mat_path);

local_print_top_self_time(p, 30);

end

%% ------------------------------------------------------------------------
function cfg = local_build_cfg()
% Minimal cfg aligned with main_coupled_simulation defaults, shortened run.

cfg.case_name = ['ProfileShort-' datestr(now, 'yyyymmdd_HHMMSS')];
cfg.minor_axis = 0.4e-6;
cfg.major_axis = 1.25e-6;
cfg.N_body = 300;
cfg.alpha_0 = 0;

pili_config(1).type = 'linear';
pili_config(1).length_um = 3;
pili_config(1).deflection_z = pi / 2;
pili_config(1).deflection_x = pi / 2;
pili_config(1).site = 3;
pili_config(1).intrinsic_curvature = 'frozen';
pili_config(1).active_moment = struct('type', 'off');

pili_config(2).type = 'linear';
pili_config(2).length_um = 3;
pili_config(2).deflection_z = pi / 2;
pili_config(2).deflection_x = pi / 2;
pili_config(2).site = 5;
pili_config(2).intrinsic_curvature = 'frozen';
pili_config(2).active_moment = struct('type', 'off');

[cfg.Pili_Matrix, cfg.pili_ext] = parse_pili_config(pili_config);

cfg.background_flow.type = 'off';
cfg.seg_per_um = 10;
cfg.E = 3.279e8;
cfg.r_pili = 3e-9;
cfg.sigma_p = 0.2;
cfg.mu = 1e-3;
cfg.rho_bac = 998.232;
cfg.rho_fluid = 998.232;

cfg.T_obs = 0.15;
cfg.AbsTol = 1e-3;
cfg.RelTol = 1e-3;
cfg.delta = pi / 20;
cfg.tlim = 600;
cfg.dt = 0.05;

cfg.enable_brown = false;
cfg.Temper = 293.15;
cfg.dt_B = cfg.dt / 2;
cfg.fixman_midpoint = false;

cfg.preview_initial = false;
cfg.generate_animation = false;
cfg.animation_fps = 30;
cfg.animation_skip = 1;
cfg.auto_confirm = true;
cfg.compute_flowfield = false;

cfg.flowfield_grid_spacing = 0.05;
cfg.flowfield_grid_range = 10.0;
cfg.flowfield_quiver_skip = 5;
cfg.flowfield_chunk_size = 2000;

end

%% ------------------------------------------------------------------------
function local_print_top_self_time(p, nTop)
if ~isstruct(p) || ~isfield(p, 'FunctionTable')
    fprintf('[PROFILE] No FunctionTable in profile info.\n');
    return;
end
T = p.FunctionTable;
n = numel(T);
if n == 0
    fprintf('[PROFILE] Empty FunctionTable.\n');
    return;
end

selfTime = zeros(n, 1);
names = cell(n, 1);
for i = 1:n
    names{i} = T(i).CompleteName;
    if isfield(T, 'SelfTime')
        selfTime(i) = T(i).SelfTime;
    elseif isfield(T, 'TotalTime') && isfield(T, 'ChildrenTime')
        selfTime(i) = T(i).TotalTime - T(i).ChildrenTime;
    else
        selfTime(i) = T(i).TotalTime;
    end
end

[~, ix] = sort(selfTime, 'descend');
nTop = min(nTop, n);

fprintf('\n========== TOP %d BY SELF TIME (CPU) ==========\n', nTop);
fprintf('  # |   self(s) |  total(s) | Function\n');
fprintf('----+-----------+-----------+----------------------------------------\n');
for k = 1:nTop
    j = ix(k);
    fprintf('%3d | %9.4f | %9.4f | %s\n', k, selfTime(j), T(j).TotalTime, names{j});
end
fprintf('================================================\n\n');
end
