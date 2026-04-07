function postprocess_cases()
% POSTPROCESS_CASES  Batch post-processing of coupled body-pili simulation results.
%
%   Usage:
%       postprocess_cases()
%
%   Workflow:
%       1. User selects a parent directory containing one or more Case subfolders.
%       2. Valid subfolders (those containing timeseries.csv) are shown in a
%          multi-select list — user chooses which cases to process.
%       3. For each selected case the script:
%            (a) Loads  timeseries.csv  and  config.json
%            (b) Computes four kinematic metrics:
%                  - Sinking velocity statistics
%                  - Maximum horizontal displacement
%                  - 3-D trajectory tortuosity
%                  - Translational diffusion coefficient  (MSD fit)
%       4. All results are exported to a professionally formatted Excel workbook.
%
%   Data source per case folder  (non-recursive — subfolders are ignored):
%       timeseries.csv  ->  t_s, X_um, Y_um, Z_um,
%                           Vx_um_s, Vy_um_s, Vz_um_s,
%                           Vsink_um_s, tilt_rad, tilt_deg
%       config.json     ->  nested JSON  (same structure as config_example.json)
%
%   Requirements: MATLAB R2018b+, Microsoft Excel installed (for COM formatting).
%
%   Author : Lucien  <lucien-6@qq.com>
%   Date   : 2026-04-07

%% ---- Banner ------------------------------------------------------------
sep = repmat('=', 1, 64);
fprintf('\n%s\n', sep);
fprintf('  Coupled Body-Pili Simulation  |  Batch Post-Processor\n');
fprintf('  Author : Lucien  <lucien-6@qq.com>\n');
fprintf('  Run at : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('%s\n\n', sep);

%% ---- Step 1: Case folder selection ------------------------------------
fprintf('[1/4] Selecting case folders ...\n');
case_folders = select_case_folders();
if isempty(case_folders)
    fprintf('      No cases selected — exiting.\n\n');
    return;
end
N = numel(case_folders);
fprintf('      %d case(s) selected.\n\n', N);

%% ---- Step 2: Output file path ----------------------------------------
fprintf('[2/4] Choosing output file path ...\n');
default_name = sprintf('postprocess_summary_%s.xlsx', datestr(now, 'yyyymmdd'));
[out_file, out_dir] = uiputfile( ...
    {'*.xlsx', 'Excel Workbook (*.xlsx)'}, ...
    'Save Summary Workbook As', default_name);
if isequal(out_file, 0)
    fprintf('      Save cancelled — exiting.\n\n');
    return;
end
out_path = fullfile(out_dir, out_file);
fprintf('      Output : %s\n\n', out_path);

%% ---- Step 3: Process loop with waitbar --------------------------------
fprintf('[3/4] Processing cases ...\n');
fprintf('      %s\n', repmat('-', 1, 56));

hWB = waitbar(0, 'Initialising ...', ...
    'Name', 'Batch Post-Processor  —  Coupled Body-Pili Simulation', ...
    'CreateCancelBtn', 'setappdata(gcbf,''cancelled'',true)');
setappdata(hWB, 'cancelled', false);
cleanupWB = onCleanup(@() close_waitbar_safe(hWB));

all_rows  = cell(N, 1);
n_ok      = 0;
n_err     = 0;

for i = 1:N
    % Cancellation check
    if getappdata(hWB, 'cancelled')
        fprintf('\n      Processing cancelled by user at case %d / %d.\n', i, N);
        break;
    end

    folder_path = case_folders{i};
    [~, folder_name] = fileparts(folder_path);

    % Progress bar — start of case
    pct_start = (i - 1) / N;
    waitbar(pct_start, hWB, ...
        sprintf('[%d / %d]  Loading  »  %s', i, N, folder_name));
    fprintf('      [%d/%d]  %s\n', i, N, folder_name);

    % Load, compute, build row
    try
        data   = load_case_data(folder_path);
        waitbar(pct_start + 0.4/N, hWB, ...
            sprintf('[%d / %d]  Computing  »  %s', i, N, folder_name));
        mets   = compute_metrics(data);
        row    = build_result_row(folder_name, folder_path, data, mets);
        row.Status = 'OK';
        n_ok   = n_ok + 1;
        fprintf('             Status : OK   |  V_sink_steady = %+.4e µm/s  |  Tortuosity = %.4f\n', ...
            row.Vsink_Steady_um_s, row.Tortuosity);
    catch ME
        row        = empty_result_row(folder_name, folder_path);
        row.Status = ['ERROR: ' ME.message];
        n_err      = n_err + 1;
        fprintf('             Status : ERROR  >>  %s\n', ME.message);
    end

    all_rows{i} = row;

    % Progress bar — end of case
    waitbar(i / N, hWB, ...
        sprintf('[%d / %d]  Done  »  %s', i, N, folder_name));
    fprintf('\n');
end

fprintf('      %s\n', repmat('-', 1, 56));
fprintf('      Completed : %d OK,  %d ERROR(S)\n\n', n_ok, n_err);

%% ---- Step 4: Export to Excel -----------------------------------------
fprintf('[4/4] Exporting to Excel ...\n');
waitbar(1, hWB, 'Writing Excel workbook  —  please wait ...');

valid_rows = all_rows(~cellfun(@isempty, all_rows));
if isempty(valid_rows)
    fprintf('      No results to export.\n\n');
    return;
end

export_to_excel(valid_rows, out_path);

%% ---- Done --------------------------------------------------------------
fprintf('\n%s\n', sep);
fprintf('  Post-processing COMPLETE\n');
fprintf('  Processed : %d / %d case(s)   (%d error(s))\n', n_ok, N, n_err);
fprintf('  Output    : %s\n', out_path);
fprintf('%s\n\n', sep);

% Ask user before opening (prevents race condition with COM cleanup)
answer = questdlg( ...
    sprintf(['Post-processing complete!\n\n' ...
    'Processed:  %d / %d case(s)   (%d error(s))\n\n' ...
    'Output file:\n%s\n\n' ...
    'Open the Excel workbook now?'], ...
    n_ok, N, n_err, out_path), ...
    'Post-Processing Complete', ...
    'Open Excel', 'Close', 'Open Excel');

if strcmp(answer, 'Open Excel')
    % Allow the COM Excel process to fully exit before opening a new instance
    pause(2);
    try
        winopen(out_path);
    catch
        fprintf('[INFO] Could not auto-open the file.\n       Path: %s\n\n', out_path);
    end
end
end


%% =========================================================================
%  LOCAL FUNCTION: select_case_folders
%  Step 1 — browse for parent directory, then multi-select valid Case folders.
%% =========================================================================
function folders = select_case_folders()
folders = {};

% Browse for parent directory
start_dir = pwd;
parent = uigetdir(start_dir, ...
    'Select the parent directory containing Case result subfolders');
if isequal(parent, 0)
    return;
end

% Scan first-level subdirectories only
d = dir(parent);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.', '..'}));

if isempty(d)
    errordlg( ...
        sprintf('No subfolders found in:\n%s', parent), ...
        'Post-Processor — No Subfolders');
    return;
end

% Keep only folders that contain timeseries.csv
valid_names = {};
valid_paths = {};
for k = 1:numel(d)
    csv_test = fullfile(parent, d(k).name, 'timeseries.csv');
    if isfile(csv_test)
        valid_names{end+1} = d(k).name;              %#ok<AGROW>
        valid_paths{end+1} = fullfile(parent, d(k).name); %#ok<AGROW>
    end
end

if isempty(valid_names)
    errordlg( ...
        sprintf(['No valid Case folders found in:\n%s\n\n' ...
        '(Each Case folder must contain timeseries.csv)'], parent), ...
        'Post-Processor — No Valid Cases');
    return;
end

% Custom multi-select dialog (fully English, live selection count)
sel_idx = show_case_selector_dlg(valid_names, parent);
if isempty(sel_idx)
    return;
end

folders = valid_paths(sel_idx);
end


%% =========================================================================
%  LOCAL FUNCTION: load_case_data
%  Reads timeseries.csv and config.json from a single Case folder.
%  Returns struct with fields: ts (table), cfg (struct), folder (char).
%% =========================================================================
function data = load_case_data(folder_path)
data.folder = folder_path;

% timeseries.csv
csv_path = fullfile(folder_path, 'timeseries.csv');
if ~isfile(csv_path)
    error('postprocess_cases:missingCSV', ...
        'timeseries.csv not found in:\n  %s', folder_path);
end
data.ts = readtable(csv_path, 'TextType', 'string');

required_cols = {'t_s', 'X_um', 'Y_um', 'Z_um', 'Vsink_um_s', 'tilt_rad'};
for rc = required_cols
    if ~ismember(rc{1}, data.ts.Properties.VariableNames)
        error('postprocess_cases:missingColumn', ...
            'Column "%s" missing in timeseries.csv', rc{1});
    end
end

% config.json (optional — use empty struct if absent)
json_path = fullfile(folder_path, 'config.json');
if isfile(json_path)
    fid = fopen(json_path, 'r', 'n', 'UTF-8');
    raw = fread(fid, '*char')';
    fclose(fid);
    data.cfg = jsondecode(raw);
else
    data.cfg = struct();
    warning('postprocess_cases:noConfig', ...
        'config.json not found in:\n  %s', folder_path);
end
end


%% =========================================================================
%  LOCAL FUNCTION: compute_metrics
%  Computes all four kinematic indicators from loaded Case data.
%% =========================================================================
function m = compute_metrics(data)
ts = data.ts;
t  = ts.t_s(:);
X  = ts.X_um(:);
Y  = ts.Y_um(:);
Z  = ts.Z_um(:);
Vs = ts.Vsink_um_s(:);
NT = numel(t);

% (1) Sinking velocity statistics
m.Vsink_Mean_um_s   = mean(Vs);
m.Vsink_Std_um_s    = std(Vs);
steady_idx          = max(1, round(0.9 * NT)) : NT;
m.Vsink_Steady_um_s = mean(Vs(steady_idx));

% (2) Maximum horizontal displacement from initial position
horiz_disp          = sqrt((X - X(1)).^2 + (Y - Y(1)).^2);
m.Max_Horiz_Disp_um = max(horiz_disp);

% (3) 3-D trajectory tortuosity  =  path arc-length / end-to-end distance
seg_len           = sqrt(diff(X).^2 + diff(Y).^2 + diff(Z).^2);
m.Path_Length_um  = sum(seg_len);
m.Direct_Dist_um  = sqrt((X(end)-X(1))^2 + (Y(end)-Y(1))^2 + (Z(end)-Z(1))^2);
if m.Direct_Dist_um > 0
    m.Tortuosity = m.Path_Length_um / m.Direct_Dist_um;
else
    m.Tortuosity = NaN;
end

% (4) Translational diffusion coefficient via time-averaged MSD fitting
%
%   TAMSD(tau) = < |r(t+tau) - r(t)|^2 >_t        (3-D, all components)
%   Model      : MSD(tau) = 6*D*tau + V^2*tau^2
%
%   Lag range  : tau = dt  ..  T_obs/4   (statistical reliability)
%   Fit method : unconstrained least-squares  A \ b
dt      = mean(diff(t));
max_lag = max(1, floor(NT / 4));
lags    = (1 : max_lag)';
tau_vec = lags * dt;
msd_vec = zeros(max_lag, 1);

for li = 1:max_lag
    lag = lags(li);
    dr2 = (X(lag+1:end) - X(1:end-lag)).^2 + ...
        (Y(lag+1:end) - Y(1:end-lag)).^2 + ...
        (Z(lag+1:end) - Z(1:end-lag)).^2;
    msd_vec(li) = mean(dr2);
end

A     = [tau_vec, tau_vec.^2];    % design matrix  (no intercept)
coeff = A \ msd_vec;

m.D_trans_um2_s = coeff(1) / 6;
m.V_msd_um_s    = sqrt(max(0, coeff(2)));

msd_pred = A * coeff;
ss_res   = sum((msd_vec - msd_pred).^2);
ss_tot   = sum((msd_vec - mean(msd_vec)).^2);
m.MSD_R2 = ternary(ss_tot > 0, 1 - ss_res/ss_tot, NaN);

% Store MSD curve for potential downstream use
m.tau_vec_s   = tau_vec;
m.msd_vec_um2 = msd_vec;
end


%% =========================================================================
%  LOCAL FUNCTION: build_result_row
%  Combines config metadata + computed metrics into one result struct.
%% =========================================================================
function row = build_result_row(folder_name, folder_path, data, mets)
cfg = data.cfg;
ts  = data.ts;

row.Case_Name     = folder_name;
row.Source_Folder = folder_path;

% Simulation info
row.T_obs_s = ts.t_s(end);
row.N_steps = height(ts);

% Body geometry
row.Major_Axis_um = safe_get(cfg, {'body','major_axis'}, NaN) * 1e6;
row.Minor_Axis_um = safe_get(cfg, {'body','minor_axis'}, NaN) * 1e6;
row.Alpha0_deg    = rad2deg(safe_get(cfg, {'body','alpha_0'}, NaN));

% Pili
pili = safe_get(cfg, {'pili'}, []);
if isstruct(pili) && numel(pili) > 0 && isfield(pili(1), 'length_um')
    row.N_Pili         = numel(pili);
    lengths            = arrayfun(@(p) safe_get(p, {'length_um'}, NaN), pili);
    row.Pili_Length_um = mat2str(lengths, 4);
else
    row.N_Pili         = 0;
    row.Pili_Length_um = 'N/A';
end
row.Pili_E_Pa = safe_get(cfg, {'material','E'}, NaN);

% Fluid & Brownian
row.Mu_Pa_s          = safe_get(cfg, {'fluid','mu'}, NaN);
row.Rho_Bac          = safe_get(cfg, {'fluid','rho_bac'}, NaN);
row.Brownian_Enabled = logical_to_str(safe_get(cfg, {'brownian','enable'}, false));

% Kinematic metrics
row.Vsink_Mean_um_s   = mets.Vsink_Mean_um_s;
row.Vsink_Steady_um_s = mets.Vsink_Steady_um_s;
row.Vsink_Std_um_s    = mets.Vsink_Std_um_s;
row.Max_Horiz_Disp_um = mets.Max_Horiz_Disp_um;
row.Path_Length_um    = mets.Path_Length_um;
row.Direct_Dist_um    = mets.Direct_Dist_um;
row.Tortuosity        = mets.Tortuosity;
row.D_trans_um2_s     = mets.D_trans_um2_s;
row.V_msd_um_s        = mets.V_msd_um_s;
row.MSD_R2            = mets.MSD_R2;

row.Status = 'OK';
end


%% =========================================================================
%  LOCAL FUNCTION: empty_result_row
%  All-NaN row for cases that could not be processed.
%% =========================================================================
function row = empty_result_row(folder_name, folder_path)
row.Case_Name          = folder_name;
row.Source_Folder      = folder_path;
row.T_obs_s            = NaN;
row.N_steps            = NaN;
row.Major_Axis_um      = NaN;
row.Minor_Axis_um      = NaN;
row.Alpha0_deg         = NaN;
row.N_Pili             = NaN;
row.Pili_Length_um     = 'N/A';
row.Pili_E_Pa          = NaN;
row.Mu_Pa_s            = NaN;
row.Rho_Bac            = NaN;
row.Brownian_Enabled   = 'N/A';
row.Vsink_Mean_um_s    = NaN;
row.Vsink_Steady_um_s  = NaN;
row.Vsink_Std_um_s     = NaN;
row.Max_Horiz_Disp_um  = NaN;
row.Path_Length_um     = NaN;
row.Direct_Dist_um     = NaN;
row.Tortuosity         = NaN;
row.D_trans_um2_s      = NaN;
row.V_msd_um_s         = NaN;
row.MSD_R2             = NaN;
row.Status             = 'ERROR';
end


%% =========================================================================
%  LOCAL FUNCTION: export_to_excel
%  Writes raw data with writecell, then applies COM formatting in one pass.
%% =========================================================================
function export_to_excel(rows, out_path)
N = numel(rows);

% ---- Column definitions  {field, label, unit, num_format} -------------
cols = { ...
    'Case_Name',         'Case Name',                '—',     ''; ...
    'Source_Folder',     'Source Folder',            '—',     ''; ...
    'T_obs_s',           'T_obs',                    's',     '0.00'; ...
    'N_steps',           'N Steps',                  '—',     '0'; ...
    'Major_Axis_um',     'Major Axis',               'µm',    '0.000'; ...
    'Minor_Axis_um',     'Minor Axis',               'µm',    '0.000'; ...
    'Alpha0_deg',        'Initial Tilt α₀',          '°',     '0.00'; ...
    'N_Pili',            'N Pili',                   '—',     '0'; ...
    'Pili_Length_um',    'Pili Length',              'µm',    ''; ...
    'Pili_E_Pa',         'Pili Young''s Modulus E',  'Pa',    '0.000E+00'; ...
    'Mu_Pa_s',           'Fluid Viscosity µ',        'Pa·s',  '0.000E+00'; ...
    'Rho_Bac',           'Bacterium Density',        'kg/m³', '0.0'; ...
    'Brownian_Enabled',  'Brownian Motion',          '—',     ''; ...
    'Vsink_Mean_um_s',   'V_sink Mean',              'µm/s',  '0.0000'; ...
    'Vsink_Steady_um_s', 'V_sink Steady-State',      'µm/s',  '0.0000'; ...
    'Vsink_Std_um_s',    'V_sink Std Dev',           'µm/s',  '0.0000'; ...
    'Max_Horiz_Disp_um', 'Max Horiz. Displacement',  'µm',    '0.0000'; ...
    'Path_Length_um',    '3D Path Length',           'µm',    '0.0000'; ...
    'Direct_Dist_um',    'Direct Distance',          'µm',    '0.0000'; ...
    'Tortuosity',        'Tortuosity',               '—',     '0.0000'; ...
    'D_trans_um2_s',     'D_trans  (MSD fit)',       'µm²/s', '0.000E+00'; ...
    'V_msd_um_s',        'V_drift  (MSD fit)',       'µm/s',  '0.0000'; ...
    'MSD_R2',            'MSD Fit R²',               '—',     '0.0000'; ...
    'Status',            'Status',                   '—',     ''; ...
    };
n_cols = size(cols, 1);

header1 = cols(:,2)';   % label row
header2 = cols(:,3)';   % unit row

% Build data cell array
data_cell = cell(N, n_cols);
for ri = 1:N
    for ci = 1:n_cols
        val = rows{ri}.(cols{ci,1});
        if isnumeric(val) && isscalar(val) && isnan(val)
            data_cell{ri,ci} = '';
        else
            data_cell{ri,ci} = val;
        end
    end
end

% Check if the file is writable (might be locked by a previous run)
if isfile(out_path)
    fid_test = fopen(out_path, 'a');
    if fid_test == -1
        error('postprocess_cases:fileLocked', ...
            ['Cannot write to:\n  %s\n\n' ...
            'The file may be open in Excel.\n' ...
            'Please close it and run the script again.'], out_path);
    end
    fclose(fid_test);
    delete(out_path);
end

% Write raw data
writecell([header1; header2; data_cell], out_path, 'Sheet', 'Summary');
fprintf('      Raw data written  (%d row(s),  %d column(s)).\n', N, n_cols);

% Apply COM formatting
try
    format_excel_com(out_path, rows, cols, N, n_cols);
    fprintf('      Excel formatting applied.\n');
catch ME_com
    fprintf('      [WARN] COM formatting skipped : %s\n', ME_com.message);
    fprintf('             Data is saved but without visual formatting.\n');
end
end


%% =========================================================================
%  LOCAL FUNCTION: format_excel_com
%  Opens the .xlsx via Windows COM and applies all visual formatting.
%  Uses try-finally to guarantee the COM process is always closed.
%% =========================================================================
function format_excel_com(out_path, rows, cols, N, n_cols)

% Ensure absolute Windows path for Excel COM
if ~java.io.File(out_path).isAbsolute()
    out_path = fullfile(pwd, out_path);
end
abs_path = strrep(out_path, '/', '\');

% Colour palette  (Excel BGR integer — convert via rgb2xl)
C_HEADER_BG  = rgb2xl([ 44,  62,  80]);   % dark navy
C_HEADER_FG  = rgb2xl([255, 255, 255]);    % white
C_UNIT_BG    = rgb2xl([236, 240, 241]);    % light grey
C_UNIT_FG    = rgb2xl([127, 140, 141]);    % grey text
C_OK_ODD     = rgb2xl([212, 239, 223]);    % pale green (odd rows)
C_OK_EVEN    = rgb2xl([232, 248, 238]);    % lighter green (even rows)
C_ERROR_BG   = rgb2xl([250, 219, 216]);    % pale red
C_WARN_BG    = rgb2xl([253, 243, 207]);    % pale amber
C_WHITE      = rgb2xl([255, 255, 255]);
C_CASE_BG    = rgb2xl([245, 245, 245]);    % very light grey for Case Name col

last_col = col_letter(n_cols);

xl = [];

try
    % ---- Launch Excel (headless) ----------------------------------------
    xl = actxserver('Excel.Application');
    xl.Visible        = false;
    xl.DisplayAlerts  = false;   % suppresses ALL dialogs incl. "file in use"
    xl.ScreenUpdating = false;

    wb = xl.Workbooks.Open(abs_path);
    ws = wb.Worksheets.Item('Summary');
    ws.Activate();

    % ---- Freeze top 2 rows ----------------------------------------------
    try
        ws.Cells.Item(3, 1).Select();
        xl.ActiveWindow.FreezePanes = true;
    catch; end   % non-fatal in headless mode

    % ---- Row 1: header labels -------------------------------------------
    r1 = ws.Range(sprintf('A1:%s1', last_col));
    r1.Interior.Color      = C_HEADER_BG;
    r1.Font.Color          = C_HEADER_FG;
    r1.Font.Bold           = true;
    r1.Font.Size           = 11;
    r1.Font.Name           = 'Arial';
    r1.HorizontalAlignment = -4108;   % xlCenter
    r1.RowHeight           = 22;

    % ---- Row 2: units ---------------------------------------------------
    r2 = ws.Range(sprintf('A2:%s2', last_col));
    r2.Interior.Color      = C_UNIT_BG;
    r2.Font.Color          = C_UNIT_FG;
    r2.Font.Italic         = true;
    r2.Font.Size           = 9;
    r2.Font.Name           = 'Arial';
    r2.HorizontalAlignment = -4108;
    r2.RowHeight           = 15;

    % ---- Data rows: colour by status + alternating ----------------------
    for ri = 1:N
        er  = ri + 2;   % Excel row index (data starts at row 3)
        rng = ws.Range(sprintf('A%d:%s%d', er, last_col, er));
        rng.Font.Name  = 'Arial';
        rng.Font.Size  = 10;
        rng.RowHeight  = 16;

        sv = rows{ri}.Status;
        if strcmp(sv, 'OK')
            if mod(ri, 2) == 1
                rng.Interior.Color = C_OK_ODD;
            else
                rng.Interior.Color = C_OK_EVEN;
            end
        elseif strncmp(sv, 'ERROR', 5)
            rng.Interior.Color = C_ERROR_BG;
        else
            rng.Interior.Color = C_WARN_BG;
        end
    end

    % ---- Number formats (column by column) ------------------------------
    for ci = 1:n_cols
        fmt = cols{ci, 4};
        if ~isempty(fmt) && N > 0
            col_rng = ws.Range( ...
                sprintf('%s3:%s%d', col_letter(ci), col_letter(ci), N+2));
            col_rng.NumberFormat = fmt;
        end
    end

    % ---- Bold Case Name column (column A, data rows) --------------------
    if N > 0
        ws.Range(sprintf('A3:A%d', N+2)).Font.Bold = true;
    end

    % ---- Source Folder column: smaller font, wrap text ------------------
    if N > 0
        sfr = ws.Range(sprintf('B3:B%d', N+2));
        sfr.Font.Size = 8;
        sfr.WrapText  = false;
    end

    % ---- Borders (applied per-side to avoid collection-level issues) ----
    all_rng = ws.Range(sprintf('A1:%s%d', last_col, N+2));
    border_ids = [7, 8, 9, 10, 11, 12];   % Left Top Bottom Right InsideV InsideH
    for bid = border_ids
        try
            brd = all_rng.Borders.Item(bid);
            brd.LineStyle = 1;   % xlContinuous
            brd.Weight    = 2;   % xlThin
        catch; end
    end

    % ---- Auto-fit column widths -----------------------------------------
    ws.Columns.AutoFit();

    % ---- Cap Source Folder column width to 50 pts -----------------------
    try ws.Columns.Item(2).ColumnWidth = 50; catch; end

    % ---- Save and close -------------------------------------------------
    wb.Save();
    wb.Close(false);
    wb = [];  

    xl.Quit();
    xl.delete();
    xl = [];  
    drawnow;

catch ME
    % Guaranteed cleanup on any error
    try
        if ~isempty(xl)
            xl.DisplayAlerts = false;
            xl.Quit();
            xl.delete();
        end
    catch; end
    rethrow(ME);
end
end


%% =========================================================================
%  HELPER: col_letter  —  1-based column index  ->  Excel column letter(s)
%% =========================================================================
function letter = col_letter(idx)
letter = '';
while idx > 0
    r      = mod(idx - 1, 26);
    letter = [char(r + 65) letter]; %#ok<AGROW>
    idx    = floor((idx - 1) / 26);
end
end


%% =========================================================================
%  HELPER: rgb2xl  —  [R G B] (0-255)  ->  Excel BGR integer
%% =========================================================================
function c = rgb2xl(rgb)
c = rgb(1) + rgb(2)*256 + rgb(3)*65536;
end


%% =========================================================================
%  HELPER: safe_get  —  safely traverse nested struct fields
%  Example:  safe_get(cfg, {'body','major_axis'}, NaN)
%% =========================================================================
function val = safe_get(s, path_cell, default_val)
val = s;
for k = 1:numel(path_cell)
    fld = path_cell{k};
    if isstruct(val) && isfield(val, fld)
        val = val.(fld);
    else
        val = default_val;
        return;
    end
end
end


%% =========================================================================
%  HELPER: logical_to_str
%% =========================================================================
function s = logical_to_str(v)
if (islogical(v) || isnumeric(v)) && isscalar(v)
    s = ternary(logical(v), 'Yes', 'No');
else
    s = 'N/A';
end
end


%% =========================================================================
%  HELPER: ternary  —  inline if-else  (avoids nested ternary anti-pattern)
%% =========================================================================
function out = ternary(cond, val_true, val_false)
if cond
    out = val_true;
else
    out = val_false;
end
end


%% =========================================================================
%  HELPER: close_waitbar_safe  —  used by onCleanup to safely close waitbar
%% =========================================================================
function close_waitbar_safe(h)
try
    if ishandle(h)
        delete(h);
    end
catch
end
end


%% =========================================================================
%  DIALOG: show_case_selector_dlg
%  Minimal, clean selector with only the essential controls.
%% =========================================================================
function sel_idx = show_case_selector_dlg(valid_names, parent_folder)

sel_idx = [];
n_cases = numel(valid_names);

FIG_W = 660;
FIG_H = 470;
PAD   = 18;
HDR_H = 66;
BTN_H = 34;

C_FIG_BG      = [246 248 251] ./ 255;
C_PANEL_BG    = [1 1 1];
C_HEADER_BG   = [35 52 69] ./ 255;
C_HEADER_FG   = [1 1 1];
C_SUBTITLE_FG = [190 203 214] ./ 255;
C_TEXT_MAIN   = [46 62 80] ./ 255;
C_TEXT_MUTED  = [108 122 137] ./ 255;
C_BORDER      = [220 226 232] ./ 255;
C_OK_BG       = [41 128 185] ./ 255;
C_OK_FG       = [1 1 1];
C_LIST_BG     = [244 247 250] ./ 255;

hFig = figure( ...
    'Name',            'Post-Processor - Select Cases', ...
    'MenuBar',         'none', ...
    'ToolBar',         'none', ...
    'NumberTitle',     'off', ...
    'Resize',          'off', ...
    'WindowStyle',     'modal', ...
    'Color',           C_FIG_BG, ...
    'Position',        [100 100 FIG_W FIG_H], ...
    'KeyPressFcn',     @(~, evt) dlg_keypress(hFig, evt));
movegui(hFig, 'center');

setappdata(hFig, 'n_cases',   n_cases);
setappdata(hFig, 'confirmed', false);

uipanel(hFig, ...
    'Units',           'pixels', ...
    'Position',        [0 FIG_H-HDR_H FIG_W HDR_H], ...
    'BackgroundColor', C_HEADER_BG, ...
    'BorderType',      'none');

uicontrol(hFig, 'Style', 'text', ...
    'String',              'Select Cases', ...
    'BackgroundColor',     C_HEADER_BG, ...
    'ForegroundColor',     C_HEADER_FG, ...
    'FontName',            'Arial', ...
    'FontSize',            15, ...
    'FontWeight',          'bold', ...
    'HorizontalAlignment', 'left', ...
    'Position',            [PAD FIG_H-40 240 22]);

uicontrol(hFig, 'Style', 'text', ...
    'String',              'Choose one or more valid case folders for post-processing.', ...
    'BackgroundColor',     C_HEADER_BG, ...
    'ForegroundColor',     C_SUBTITLE_FG, ...
    'FontName',            'Arial', ...
    'FontSize',            9.5, ...
    'HorizontalAlignment', 'left', ...
    'Position',            [PAD FIG_H-60 420 18]);

uipanel(hFig, ...
    'Units',           'pixels', ...
    'Position',        [PAD 74 FIG_W - 2*PAD FIG_H - HDR_H - 90], ...
    'BackgroundColor', C_PANEL_BG, ...
    'BorderType',      'line', ...
    'BorderColor',     C_BORDER);

uicontrol(hFig, 'Style', 'text', ...
    'String',              'Parent Folder', ...
    'BackgroundColor',     C_PANEL_BG, ...
    'ForegroundColor',     C_TEXT_MUTED, ...
    'FontName',            'Arial', ...
    'FontSize',            9, ...
    'FontWeight',          'bold', ...
    'HorizontalAlignment', 'left', ...
    'Position',            [PAD+16 FIG_H-HDR_H-34 FIG_W-2*PAD-32 16]);

uicontrol(hFig, 'Style', 'text', ...
    'String',              truncate_middle_path(parent_folder, 88), ...
    'BackgroundColor',     C_PANEL_BG, ...
    'ForegroundColor',     C_TEXT_MAIN, ...
    'FontName',            'Arial', ...
    'FontSize',            9.5, ...
    'HorizontalAlignment', 'left', ...
    'Position',            [PAD+16 FIG_H-HDR_H-56 FIG_W-2*PAD-32 18]);

hList = uicontrol(hFig, 'Style', 'listbox', ...
    'String',              valid_names, ...
    'Max',                 n_cases, ...
    'Min',                 0, ...
    'Value',               [], ...
    'FontName',            'Arial', ...
    'FontSize',            10.5, ...
    'BackgroundColor',     C_LIST_BG, ...
    'ForegroundColor',     C_TEXT_MAIN, ...
    'Position',            [PAD+16 138 FIG_W-2*PAD-32 FIG_H-HDR_H-218], ...
    'Callback',            @(~,~) dlg_list_changed(hFig));
setappdata(hFig, 'hList', hList);

hStatus = uicontrol(hFig, 'Style', 'text', ...
    'String',              sprintf('0 of %d case(s) selected', n_cases), ...
    'BackgroundColor',     C_PANEL_BG, ...
    'ForegroundColor',     C_TEXT_MUTED, ...
    'FontName',            'Arial', ...
    'FontSize',            9, ...
    'HorizontalAlignment', 'left', ...
    'Position',            [PAD+16 102 220 18]);
setappdata(hFig, 'hStatus', hStatus);

uicontrol(hFig, 'Style', 'text', ...
    'String',              'Use Ctrl / Shift for multi-selection.', ...
    'BackgroundColor',     C_PANEL_BG, ...
    'ForegroundColor',     C_TEXT_MUTED, ...
    'FontName',            'Arial', ...
    'FontSize',            8.5, ...
    'HorizontalAlignment', 'right', ...
    'Position',            [FIG_W-2*PAD-220 102 220 18]);

btn_gap   = 10;
btn_w_sel = 88;
btn_w_clr = 88;
btn_w_can = 76;
btn_w_run = 142;
btn_x_run = FIG_W - PAD - btn_w_run;
btn_x_can = btn_x_run - btn_gap - btn_w_can;
btn_x_clr = btn_x_can - btn_gap - btn_w_clr;
btn_x_sel = btn_x_clr - btn_gap - btn_w_sel;

uicontrol(hFig, 'Style', 'pushbutton', ...
    'String',          'Select All', ...
    'FontName',        'Arial', ...
    'FontSize',        9.5, ...
    'Position',        [btn_x_sel 24 btn_w_sel BTN_H], ...
    'Callback',        @(~,~) dlg_btn_select_all(hFig));

uicontrol(hFig, 'Style', 'pushbutton', ...
    'String',          'Clear All', ...
    'FontName',        'Arial', ...
    'FontSize',        9.5, ...
    'Position',        [btn_x_clr 24 btn_w_clr BTN_H], ...
    'Callback',        @(~,~) dlg_btn_clear_all(hFig));

uicontrol(hFig, 'Style', 'pushbutton', ...
    'String',          'Cancel', ...
    'FontName',        'Arial', ...
    'FontSize',        10, ...
    'Position',        [btn_x_can 24 btn_w_can BTN_H], ...
    'Callback',        @(~,~) dlg_btn_cancel(hFig));

uicontrol(hFig, 'Style', 'pushbutton', ...
    'String',          'Process Selected', ...
    'FontName',        'Arial', ...
    'FontSize',        10, ...
    'FontWeight',      'bold', ...
    'ForegroundColor', C_OK_FG, ...
    'BackgroundColor', C_OK_BG, ...
    'Position',        [btn_x_run 24 btn_w_run BTN_H], ...
    'Callback',        @(~,~) dlg_btn_ok(hFig));

set(hFig, 'CloseRequestFcn', @(~,~) dlg_btn_cancel(hFig));
dlg_list_changed(hFig);

uiwait(hFig);

if ishandle(hFig) && getappdata(hFig, 'confirmed')
    hList = getappdata(hFig, 'hList');
    sel_idx = sort(get(hList, 'Value'));
end
if ishandle(hFig)
    delete(hFig);
end
end


%% =========================================================================
%  DIALOG CALLBACKS
%% =========================================================================

function dlg_list_changed(hFig)
hList = getappdata(hFig, 'hList');
hStatus = getappdata(hFig, 'hStatus');
n_cases = getappdata(hFig, 'n_cases');
n_sel = numel(get(hList, 'Value'));
set(hStatus, 'String', sprintf('%d of %d case(s) selected', n_sel, n_cases));
end

function dlg_btn_select_all(hFig)
hList = getappdata(hFig, 'hList');
n_cases = getappdata(hFig, 'n_cases');
set(hList, 'Value', 1:n_cases);
dlg_list_changed(hFig);
end

function dlg_btn_clear_all(hFig)
hList = getappdata(hFig, 'hList');
set(hList, 'Value', []);
dlg_list_changed(hFig);
end

function dlg_btn_ok(hFig)
hList = getappdata(hFig, 'hList');
if isempty(get(hList, 'Value'))
    uiwait(msgbox('Please select at least one case before proceeding.', ...
        'No Case Selected', 'warn'));
    return;
end
setappdata(hFig, 'confirmed', true);
uiresume(hFig);
end

function dlg_btn_cancel(hFig)
setappdata(hFig, 'confirmed', false);
if strcmp(get(hFig, 'BeingDeleted'), 'off')
    uiresume(hFig);
end
end

function dlg_keypress(hFig, evt)
switch evt.Key
    case {'return', 'space'}
        dlg_btn_ok(hFig);
    case 'escape'
        dlg_btn_cancel(hFig);
end
end


%% =========================================================================
%  DIALOG HELPERS
%% =========================================================================

function out = truncate_middle_path(in, max_len)
if numel(in) <= max_len
    out = in;
    return;
end

head_len = floor((max_len - 3) / 2);
tail_len = max_len - 3 - head_len;
out = [in(1:head_len) '...' in(end-tail_len+1:end)];
end
