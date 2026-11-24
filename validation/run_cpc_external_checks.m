function summary = run_cpc_external_checks(varargin)
%RUN_CPC_EXTERNAL_CHECKS  Drivers for the two external CPC validation cases
%
%  Description:
%    Runs cpc_helix_propulsive_matrix and cpc_cantilever_relaxation and
%    writes both tables under output/cpc_checks/.  A failing pass line is
%    reported in the summary; the functions themselves do not throw.
%
%  Syntax:
%    summary = run_cpc_external_checks()
%    summary = run_cpc_external_checks('IncludeBenchmark', true)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-15
%  Updated: 2026-11-24 — optional IncludeBenchmark for v10.0.0 performance suite.

ip = inputParser;
addParameter(ip, 'OutDir', '');
addParameter(ip, 'IncludeBenchmark', false, @(x) islogical(x) || isnumeric(x));
parse(ip, varargin{:});
if isempty(ip.Results.OutDir)
    args = {};
else
    args = {'OutDir', ip.Results.OutDir};
end

summary.helix = cpc_helix_propulsive_matrix(args{:});
summary.cantilever = cpc_cantilever_relaxation(args{:});
summary.passed = summary.helix.passed && summary.cantilever.passed;
if ip.Results.IncludeBenchmark
    summary.benchmark = cpc_performance_benchmark(args{:});
    summary.passed = summary.passed && summary.benchmark.passed;
end
fprintf('External CPC checks %s.\n', tern(summary.passed, 'PASSED', 'FAILED'));
end


function word = tern(flag, yes, no)
if flag
    word = yes;
else
    word = no;
end
end
