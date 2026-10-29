function spec = flowfield_state_data(results, params, index, kind)
%FLOWFIELD_STATE_DATA  Restore one saved output as a flow-field input
%
%  Description:
%    Checks the version-1 flow_state history and returns the complete state,
%    a parameter copy, and the one recorded Brownian force for that output.
%    Every output, including the last, uses flow_state.brown_force(:,index).
%
%  Syntax:
%    spec = flowfield_state_data(results, params, index, kind)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28

if ~isfield(results, 'flow_state') || ~isstruct(results.flow_state) || ...
        ~isfield(results.flow_state, 'version') || results.flow_state.version ~= 1
    error('flowfield_state_data:MissingFlowState', ...
        ['This results file has no version-1 flow_state history. ', ...
         'Historical flow fields require a simulation saved by the current program.']);
end

fs = results.flow_state;
required = {'Z', 'rot', 'brown_force', ...
    'brownian_endpoint_convention', 'state_dimension', 'n_times'};
missing = required(~isfield(fs, required));
if ~isempty(missing)
    error('flowfield_state_data:MissingFields', ...
        'flow_state is missing field(s): %s.', strjoin(missing, ', '));
end

if ~isfield(results, 't') || ~isfield(results, 't_dimless')
    error('flowfield_state_data:MissingTime', ...
        'results.t and results.t_dimless are required.');
end

n_times = fs.n_times;
if ~isscalar(index) || ~isfinite(index) || index ~= round(index) || ...
        index < 1 || index > n_times
    error('flowfield_state_data:BadIndex', ...
        'Output index must be an integer from 1 to %d.', n_times);
end
index = double(index);

if numel(results.t) ~= n_times || numel(results.t_dimless) ~= n_times || ...
        size(fs.Z, 2) ~= n_times || size(fs.brown_force, 2) ~= n_times || ...
        size(fs.Z, 1) ~= fs.state_dimension || size(fs.brown_force, 1) ~= 6
    error('flowfield_state_data:InconsistentHistory', ...
        'flow_state arrays do not match the saved time vector.');
end

expected_nZ = 7 + 3 * sum(params.N_segments);
if fs.state_dimension ~= expected_nZ
    error('flowfield_state_data:StateDimension', ...
        'flow_state.state_dimension = %d, expected %d from params.', ...
        fs.state_dimension, expected_nZ);
end

if ~iscell(fs.rot) || numel(fs.rot) ~= params.M_pili
    error('flowfield_state_data:RotCount', ...
        'flow_state.rot must contain one basis per pilus.');
end
for k = 1:params.M_pili
    Rk = fs.rot{k};
    if size(Rk, 1) ~= 3 || size(Rk, 2) ~= 3 || size(Rk, 3) ~= n_times
        error('flowfield_state_data:RotShape', ...
            'flow_state.rot{%d} must be 3 x 3 x %d.', k, n_times);
    end
    R_now = Rk(:, :, index);
    if any(~isfinite(R_now(:)))
        error('flowfield_state_data:NonFiniteBasis', ...
            'Pilus %d basis at output %d is not finite.', k, index);
    end
end

if any(~isfinite(fs.Z(:, index))) || any(~isfinite(fs.brown_force(:, index)))
    error('flowfield_state_data:NonFiniteState', ...
        'Saved state or Brownian force at output %d is not finite.', index);
end

kind = lower(char(kind));
if strcmp(kind, 'final') && index ~= n_times
    error('flowfield_state_data:FinalIndex', ...
        'The final output is index %d.', n_times);
elseif ~strcmp(kind, 'history') && ~strcmp(kind, 'final')
    error('flowfield_state_data:BadKind', ...
        'kind must be ''history'' or ''final''.');
end
force = fs.brown_force(:, index);
convention = fs.brownian_endpoint_convention;

params_out = params;
if isfield(params_out, 'brown_enabled') && ~params_out.brown_enabled
    force = zeros(6, 1);
end
params_out.brown_force_dimless = force(:);
for k = 1:params.M_pili
    params_out.rot{k} = fs.rot{k}(:, :, index);
end
if ~isfield(params_out, 'logf') || isempty(params_out.logf)
    params_out.logf = @(varargin) [];
end
if ~isfield(params_out, 'logf_file') || isempty(params_out.logf_file)
    params_out.logf_file = @(varargin) [];
end

spec.Z           = fs.Z(:, index);
spec.params      = params_out;
spec.t_dimless   = results.t_dimless(index);
spec.t_s         = results.t(index);
spec.index       = index;
spec.convention  = convention;

end
