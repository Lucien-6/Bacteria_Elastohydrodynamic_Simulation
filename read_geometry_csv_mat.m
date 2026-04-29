function pts = read_geometry_csv_mat(file_path)
%READ_GEOMETRY_CSV_MAT  Read 3D point cloud data from CSV or MAT files
%
%  Description:
%    Reads point cloud data from a user-supplied geometry file in either
%    CSV or MAT format and returns an N-by-3 matrix of 3D coordinates.
%
%    CSV format:
%      - Three columns: X, Y, Z (header row is optional).
%      - Comment lines starting with '#' or '%' are skipped.
%      - Blank/empty lines are skipped.
%      - If the first row is non-numeric (header), it is automatically
%        detected and skipped.
%
%    MAT format:
%      - Must contain exactly one numeric variable (the first one found).
%      - The variable must have 2 or 3 columns; if 2 columns, a zero
%        third column is appended (Z = 0).
%      - A 3-by-N layout (3 rows, N > 3 columns) is also accepted
%        (transposed automatically).
%      - A 3-by-3 matrix is ambiguous; a warning is issued and the data
%        is interpreted as N-by-3 (3 points).
%
%  Syntax:
%    pts = read_geometry_csv_mat(file_path)
%
%  Input:
%    file_path - Path to a CSV or MAT file containing point cloud data
%
%  Output:
%    pts - [N x 3] matrix of 3D point coordinates
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-28

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
        pts = read_mat_points(file_path);
    otherwise
        error('read_geometry_csv_mat:UnsupportedFormat', ...
              ['Unsupported file format ''%s''. ', ...
               'Expected .csv, .txt, or .mat.'], ext);
end

%% Final validation: must be N x 3 finite numeric
if size(pts, 2) ~= 3
    error('read_geometry_csv_mat:InvalidDimensions', ...
          'Point data must have 3 columns, got %d.', size(pts, 2));
end
if any(~isfinite(pts(:)))
    error('read_geometry_csv_mat:NonFiniteValues', ...
          'Point data contains NaN or Inf values.');
end
if size(pts, 1) < 3
    error('read_geometry_csv_mat:InsufficientPoints', ...
          'Point data must contain at least 3 points, got %d.', size(pts, 1));
end

end


%% ==================== Local Helper Functions ====================

function pts = read_csv_points(file_path)
%READ_CSV_POINTS  Read points from CSV/TXT file with auto-header detection

    % Read raw lines (text mode, UTF-8)
    fid = fopen(file_path, 'r', 'n', 'UTF-8');
    if fid == -1
        error('read_geometry_csv_mat:CannotOpenCSV', ...
              'Cannot open CSV file: %s', file_path);
    end

    raw_lines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), continue; end
        line = strtrim(line);
        % Skip blank lines
        if isempty(line), continue; end
        % Skip comment lines
        if startsWith(line, '#') || startsWith(line, '%'), continue; end
        raw_lines{end+1} = line; %#ok<AGROW>
    end
    fclose(fid);

    if isempty(raw_lines)
        error('read_geometry_csv_mat:EmptyCSV', ...
              'No data lines found in CSV file: %s', file_path);
    end

    % Detect header: if the first line cannot be parsed as 3 numbers,
    % treat it as a header and skip it.
    start_row = 1;
    first_vals = sscanf(raw_lines{1}, '%f');
    if numel(first_vals) < 3
        start_row = 2;
    end

    % Parse numeric data
    all_data = [];
    for i = start_row:numel(raw_lines)
        vals = sscanf(raw_lines{i}, '%f');
        if numel(vals) >= 3
            all_data = [all_data; vals(1:3)']; %#ok<AGROW>
        end
    end

    pts = all_data;

end


function pts = read_mat_points(file_path)
%READ_MAT_POINTS  Read points from MAT file

    data = load(file_path);
    fields = fieldnames(data);

    if isempty(fields)
        error('read_geometry_csv_mat:EmptyMAT', ...
              'MAT file contains no variables: %s', file_path);
    end

    % Use the first variable found
    M = data.(fields{1});

    if ~isnumeric(M)
        error('read_geometry_csv_mat:NonNumericMAT', ...
              'Variable ''%s'' in MAT file is not numeric.', fields{1});
    end

    % Accept either N-by-3 or 3-by-N layouts
    if size(M, 1) == 3 && size(M, 2) > 3
        pts = M';    % Transpose 3-by-N to N-by-3 (unambiguous)
    elseif size(M, 1) == 3 && size(M, 2) == 3
        % Ambiguous: both 3 points (3x3) and 3 coordinates (3x3) are possible.
        % Assume N-by-3 (3 points in 3D) as the more common use case.
        warning('read_geometry_csv_mat:AmbiguousShape', ...
            ['Variable ''%s'' is 3-by-3, which is ambiguous. ', ...
             'Interpreting as 3 points (N-by-3). If you intended 3-by-N, ', ...
             'transpose the data before saving.'], fields{1});
        pts = M;
    elseif size(M, 2) == 3
        pts = M;     % Already N-by-3
    elseif size(M, 2) == 2
        pts = [M, zeros(size(M, 1), 1)];   % Append Z = 0
    else
        error('read_geometry_csv_mat:InvalidMATShape', ...
              ['Variable ''%s'' has size [%d x %d]. ', ...
               'Expected N-by-3, 3-by-N, or N-by-2.'], ...
               fields{1}, size(M, 1), size(M, 2));
    end

end
