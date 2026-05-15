function pts = read_geometry_csv_mat(file_path, varargin)
%READ_GEOMETRY_CSV_MAT  Read 3D point cloud data from CSV or MAT files
%
%  Description:
%    Reads point cloud data from a user-supplied geometry file in either
%    CSV or MAT format and returns an N-by-3 matrix of 3D coordinates.
%
%    CSV format (.csv / .txt):
%      - Three columns: X, Y, Z (header row optional).
%      - Comma, space, or tab delimited.
%      - Comment lines starting with '#' or '%' are skipped.
%      - Blank/empty lines are skipped.
%      - Header auto-detection: if the first non-comment row contains any
%        non-numeric token, it is treated as a header and skipped.
%
%    MAT format (.mat):
%      - Prefers a variable named 'pts', 'points', 'xyz', 'coords', or
%        'coordinates' (case-insensitive).
%      - If no preferred name exists, uses the first numeric variable and
%        issues a warning listing all available variable names.
%      - Accepts N-by-3 or 3-by-N layouts (transposed automatically when
%        one dimension is clearly 3 and the other is larger).
%      - 3-by-3 matrices are REJECTED (ambiguous); user must transpose.
%      - N-by-2 (planar) data is REJECTED by default; pass
%        'AllowPlanar2D', true to permit Z=0 padding (v7.1 stricter).
%
%  Syntax:
%    pts = read_geometry_csv_mat(file_path)
%    pts = read_geometry_csv_mat(file_path, 'Name', Value, ...)
%
%  Inputs:
%    file_path - Path to a CSV or MAT file containing point cloud data
%
%  Name-Value Pairs:
%    'AllowPlanar2D' - Logical (default false): accept N-by-2 input and
%                      auto-append Z=0 column (with a warning).
%
%  Output:
%    pts - [N x 3] matrix of 3D point coordinates
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-28
%  Updated: 2026-05-04 — Stricter shape validation, robust CSV reader via
%                        readmatrix, preferred MAT variable names.

%% Parse options
p = inputParser;
addParameter(p, 'AllowPlanar2D', false, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

%% Validate file existence
if ~isfile(file_path)
    error('read_geometry_csv_mat:FileNotFound', ...
          'Geometry file not found: %s', file_path);
end

%% Dispatch by file extension
[~, ~, ext] = fileparts(file_path);
ext = lower(ext);

switch ext
    case {'.csv', '.txt'}
        pts = read_csv_points(file_path);
    case '.mat'
        pts = read_mat_points(file_path, opts);
    otherwise
        error('read_geometry_csv_mat:UnsupportedFormat', ...
              ['Unsupported file format ''%s''. ', ...
               'Expected .csv, .txt, or .mat.'], ext);
end

%% Shape normalisation (N-by-3 / N-by-2 / 3-by-N)
pts = normalise_shape(pts, file_path, opts);

%% Final validation: must be N x 3 finite numeric
if ~isnumeric(pts) || ~ismatrix(pts)
    error('read_geometry_csv_mat:NonNumeric', ...
          'Point data in %s is not a numeric 2-D matrix.', file_path);
end
if size(pts, 2) ~= 3
    error('read_geometry_csv_mat:InvalidDimensions', ...
          'Point data in %s must have 3 columns, got %d.', ...
          file_path, size(pts, 2));
end
if any(~isfinite(pts(:)))
    error('read_geometry_csv_mat:NonFiniteValues', ...
          'Point data in %s contains NaN or Inf values.', file_path);
end
if size(pts, 1) < 3
    error('read_geometry_csv_mat:InsufficientPoints', ...
          'Point data in %s must contain at least 3 points, got %d.', ...
          file_path, size(pts, 1));
end

end


%% ==================== Local Helper Functions ====================

function pts = read_csv_points(file_path)
%READ_CSV_POINTS  Read points from CSV/TXT file robustly.
%
%  Strategy:
%    1. Try readmatrix with auto delimiter / header detection.
%    2. On failure, fall back to a line-by-line reader that strips
%       comments (#, %), unifies delimiters, and uses sscanf.

    % ---- Primary path: readmatrix (handles auto delimiter + header) ----
    try
        opts = detectImportOptions(file_path, ...
            'FileType',        'text', ...
            'CommentStyle',    {'#', '%'}, ...
            'ExpectedNumVariables', 3);
        % Force all columns numeric; readmatrix returns NaN on non-numeric
        opts.VariableTypes(:) = {'double'};
        M = readmatrix(file_path, opts);

        % Drop rows where every column is NaN (typically header line)
        all_nan_rows = all(isnan(M), 2);
        M(all_nan_rows, :) = [];

        % Drop any remaining row containing NaN (mixed row — malformed)
        M(any(isnan(M), 2), :) = [];

        if ~isempty(M)
            pts = M;
            return;
        end
    catch
        % fall through to manual parser
    end

    % ---- Fallback path: manual parser ----
    fid = fopen(file_path, 'r', 'n', 'UTF-8');
    if fid == -1
        error('read_geometry_csv_mat:CannotOpenCSV', ...
              'Cannot open CSV file: %s', file_path);
    end
    cleanup_fid = onCleanup(@() fclose(fid));

    raw_lines = cell(0, 1);
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), continue; end
        line = strtrim(line);
        if isempty(line),                            continue; end
        if startsWith(line, '#') || startsWith(line, '%'), continue; end
        raw_lines{end+1, 1} = line; %#ok<AGROW>
    end

    if isempty(raw_lines)
        error('read_geometry_csv_mat:EmptyCSV', ...
              'No data lines found in CSV file: %s', file_path);
    end

    % Normalise separators: commas and tabs → spaces
    raw_lines = regexprep(raw_lines, '[,\t]', ' ');

    % Detect header: first row not parseable as >=3 numerics is treated as header
    first_vals = sscanf(raw_lines{1}, '%f');
    start_row  = 1;
    if numel(first_vals) < 3
        start_row = 2;
    end

    % Pre-allocate exactly once (vectorised)
    n_rows = numel(raw_lines) - start_row + 1;
    pts    = nan(n_rows, 3);
    row_i  = 0;
    for i = start_row:numel(raw_lines)
        vals = sscanf(raw_lines{i}, '%f');
        if numel(vals) >= 3
            row_i = row_i + 1;
            pts(row_i, :) = vals(1:3)';
        end
    end
    pts = pts(1:row_i, :);

end


function pts = read_mat_points(file_path, ~)
%READ_MAT_POINTS  Read points from MAT file, preferring known variable names.

    data   = load(file_path);
    fields = fieldnames(data);

    if isempty(fields)
        error('read_geometry_csv_mat:EmptyMAT', ...
              'MAT file contains no variables: %s', file_path);
    end

    % Preferred variable names (case-insensitive)
    preferred = {'pts', 'points', 'xyz', 'coords', 'coordinates'};
    chosen = '';
    for kk = 1:numel(fields)
        if any(strcmpi(fields{kk}, preferred))
            chosen = fields{kk};
            break;
        end
    end

    if isempty(chosen)
        % Fall back to first numeric variable
        for kk = 1:numel(fields)
            if isnumeric(data.(fields{kk}))
                chosen = fields{kk};
                break;
            end
        end
        if numel(fields) > 1
            warning('read_geometry_csv_mat:MultipleMATVars', ...
                ['MAT file %s contains %d variables [%s]; ', ...
                 'using ''%s''. Rename your variable to ''pts'' / ', ...
                 '''points'' / ''XYZ'' / ''coords'' to suppress this warning.'], ...
                 file_path, numel(fields), strjoin(fields, ', '), chosen);
        end
    end

    if isempty(chosen)
        error('read_geometry_csv_mat:NoNumericMATVar', ...
              'MAT file %s contains no numeric variable.', file_path);
    end

    M = data.(chosen);
    if ~isnumeric(M)
        error('read_geometry_csv_mat:NonNumericMAT', ...
              'Variable ''%s'' in MAT file %s is not numeric.', ...
              chosen, file_path);
    end
    pts = M;

end


function pts = normalise_shape(M, file_path, opts)
%NORMALISE_SHAPE  Reduce N-by-3 / 3-by-N / N-by-2 inputs to N-by-3.

    if ~isnumeric(M) || ~ismatrix(M)
        error('read_geometry_csv_mat:NonNumericMatrix', ...
              'Input data from %s is not a 2-D numeric matrix.', file_path);
    end

    [nr, nc] = size(M);

    if nr == 3 && nc == 3
        error('read_geometry_csv_mat:Ambiguous3x3', ...
            ['Data from %s is 3-by-3, which is ambiguous between ', ...
             '3 points (N-by-3) and 3 coordinates (3-by-N). ', ...
             'Resolve by supplying at least 4 points or transposing.'], ...
             file_path);
    end

    if nc == 3
        pts = M;
        return;
    end

    if nr == 3 && nc > 3
        pts = M.';   % 3-by-N -> N-by-3 (unambiguous)
        return;
    end

    if nc == 2
        if opts.AllowPlanar2D
            warning('read_geometry_csv_mat:PlanarInput', ...
                ['Data from %s has 2 columns; appending Z=0 because ', ...
                 'AllowPlanar2D was enabled. This is typically a mistake ', ...
                 'for 3-D geometry.'], file_path);
            pts = [M, zeros(nr, 1)];
            return;
        else
            error('read_geometry_csv_mat:PlanarInputRejected', ...
                ['Data from %s has 2 columns (planar). 3-D geometry is ', ...
                 'required. If this is intentional, call ', ...
                 'read_geometry_csv_mat(..., ''AllowPlanar2D'', true).'], ...
                 file_path);
        end
    end

    error('read_geometry_csv_mat:InvalidShape', ...
        ['Data from %s has size [%d x %d]. ', ...
         'Expected N-by-3, 3-by-N, or N-by-2 (with AllowPlanar2D).'], ...
         file_path, nr, nc);

end
