function [spacing_um, length_um] = flowfield_stream_overrides(flow_opts)
%FLOWFIELD_STREAM_OVERRIDES  Positive streamline spacing and length, or unset
%
%  Description:
%    Reads flow_opts.stream_spacing_um and flow_opts.stream_length_um.
%    A missing or empty field is unset and returns [].  A supplied value
%    must be a positive finite scalar.  When both are supplied, the length
%    must be at least half the spacing, otherwise every curve would be
%    discarded as a stub.
%
%  Syntax:
%    [spacing_um, length_um] = flowfield_stream_overrides(flow_opts)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-13

if nargin < 1 || isempty(flow_opts)
    flow_opts = struct();
end
spacing_um = read_positive(flow_opts, 'stream_spacing_um');
length_um = read_positive(flow_opts, 'stream_length_um');
if ~isempty(spacing_um) && ~isempty(length_um) && length_um < 0.5 * spacing_um
    error('flowfield_stream_overrides:BadLength', ...
        ['stream_length_um (%.4g) must be at least half of ', ...
         'stream_spacing_um (%.4g).'], length_um, spacing_um);
end
end


function value = read_positive(s, name)
value = [];
if ~isstruct(s) || ~isfield(s, name) || isempty(s.(name))
    return
end
value = s.(name);
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
    error('flowfield_stream_overrides:BadValue', ...
        '%s must be a positive finite scalar.', name);
end
value = double(value);
end
