function s = normalise_json_struct_array(v)
%NORMALISE_JSON_STRUCT_ARRAY  Convert jsondecode output to a struct array
%
%  Description:
%    MATLAB's jsondecode returns a struct array only if every element of
%    a JSON array has exactly the same fields in the same order.  When
%    optional fields (e.g. 'label') are present on some entries and not
%    others, the decoder returns a heterogeneous cell array of structs.
%    This helper unifies the two cases into a single struct array by
%    taking the union of fields across all elements and filling missing
%    ones with [].
%
%  Syntax:
%    s = normalise_json_struct_array(v)
%
%  Inputs:
%    v - Any value; typically the result of jsondecode on a JSON array.
%        Accepted shapes:
%          - [] (empty):           returns empty struct array
%          - struct array:         returned as-is
%          - cell array of struct: merged into union struct array
%          - single struct:        wrapped to a 1-element struct array
%
%  Output:
%    s - [1 x N] struct array whose field set is the union of all input
%        elements' fields.  Missing fields are filled with [].
%
%  Example:
%    % JSON: [{"a":1,"b":2}, {"a":3}] -> cell array from jsondecode
%    v = {struct('a', 1, 'b', 2); struct('a', 3)};
%    s = normalise_json_struct_array(v);
%    % s(1).a = 1, s(1).b = 2, s(2).a = 3, s(2).b = []
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-05-04

if isempty(v)
    s = struct([]);
    return;
end

if isstruct(v)
    s = v(:)';           % ensure row struct array
    return;
end

if ~iscell(v)
    error('normalise_json_struct_array:UnsupportedType', ...
        'Input must be a struct array, cell array of structs, or empty; got %s.', ...
        class(v));
end

% Validate every element is a struct
for ii = 1:numel(v)
    if ~isstruct(v{ii})
        error('normalise_json_struct_array:NonStructCellElement', ...
            'Cell element %d is %s, expected struct.', ii, class(v{ii}));
    end
end

% Collect union of field names, preserving first-seen order
all_fields = {};
for ii = 1:numel(v)
    fn = fieldnames(v{ii});
    for kk = 1:numel(fn)
        if ~any(strcmp(fn{kk}, all_fields))
            all_fields{end+1} = fn{kk}; %#ok<AGROW>
        end
    end
end

% Build uniform struct array with missing fields filled by []
N = numel(v);
s(N) = struct();
for ii = 1:N
    entry = v{ii};
    for kk = 1:numel(all_fields)
        fname = all_fields{kk};
        if isfield(entry, fname)
            s(ii).(fname) = entry.(fname);
        else
            s(ii).(fname) = [];
        end
    end
end
s = s(:)';

end
