function compare_profiles()
%COMPARE_PROFILES  Side-by-side comparison of two profile snapshots

before_file = 'output/profile_runs/profile_info_20260421_144222.mat';
after_file  = 'output/profile_runs/profile_info_20260421_162912.mat';

S_before = load(before_file);
S_after  = load(after_file);

p0 = S_before.p;   % pre-optimisation
p1 = S_after.p;    % post-optimisation

targets = { ...
    'build_coupled_lin_sys', ...
    'build_global_resistance', ...
    'ode_coupled_rhs', ...
    'run_simulation', ...
    'cross', ...
    'blkdiag'};

fprintf('\n=============== Profile comparison ==============\n');
fprintf('%-26s  %10s  %10s  %10s\n', 'Function', 'before(s)', 'after(s)', 'delta(s)');
fprintf('%s\n', repmat('-',1,64));

for i = 1:numel(targets)
    t0 = get_self_time(p0, targets{i});
    t1 = get_self_time(p1, targets{i});
    fprintf('%-26s  %10.3f  %10.3f  %10.3f\n', targets{i}, t0, t1, t1 - t0);
end

% Overall wall-time proxy: TotalTime of the top-level function
fprintf('%s\n', repmat('-',1,64));
t0 = p0.FunctionTable(1).TotalTime;
t1 = p1.FunctionTable(1).TotalTime;
fprintf('%-26s  %10.3f  %10.3f  %10.3f  (%.1f%% reduction)\n', ...
    'TOP total', t0, t1, t1 - t0, 100*(t0 - t1)/t0);
fprintf('=================================================\n\n');

end

function t = get_self_time(info, name)
ft = info.FunctionTable;
t  = 0;
for i = 1:numel(ft)
    fn = ft(i).FunctionName;
    if contains(fn, [name, '>']) || strcmp(fn, name)
        t = t + ft(i).TotalTime;  % self time = TotalTime for leaf-like funcs
    end
end
end
