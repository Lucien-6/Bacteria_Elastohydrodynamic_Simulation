function Pili = custom_pili_building(custom_pili_cfg, attachment_sites, ...
        seg_per_um, L_0, dis, E, r_pili) %#ok<INUSD>
%CUSTOM_PILI_BUILDING  Initialize custom pili from user-provided centerline files
%
%  Description:
%    Constructs pili structures from user-provided centerline point files,
%    applying the same output format as pili_building.m for downstream
%    compatibility.  For each pilus:
%      1. Read and deduplicate raw centerline points.
%      2. Enforce the "uniform segment length" constraint: the user input
%         MUST discretize the centerline such that every segment has the
%         same length (relative tolerance 1e-3).  Non-uniform input is
%         rejected with an error.
%      3. Resample to the target segment density seg_per_um via PCHIP.
%      4. Align the root tangent to the attachment site direction (with
%         logging of the alignment angle; warn if > 45 deg).
%      5. Apply optional local-basis deflection rotations around the
%         tangent (deflection_z) and bending-plane normal (deflection_x).
%      6. Translate to (site.position + dis * site.direction).
%
%  Syntax:
%    Pili = custom_pili_building(custom_pili_cfg, attachment_sites, ...
%             seg_per_um, L_0, dis, E, r_pili)
%
%  Inputs:
%    custom_pili_cfg   - [M x 1] struct array (or cell array), one per pilus:
%        .points_file      - Path to CSV/MAT file with centerline [m]
%        .site_index       - Index into attachment_sites struct array
%        .deflection_z     - Rotation about the aligned tangent [rad]
%                             (default: 0). Applied in the pilus local basis.
%        .deflection_x     - Rotation about the bending-plane normal [rad]
%                             (default: 0). Applied in the pilus local basis.
%        .bending_stiffness - Custom A [N*m^2] (optional, default: E*I)
%    attachment_sites  - Struct array from custom_body_building:
%        .position         - [3 x 1] Root position [m]
%        .direction        - [3 x 1] Unit outward direction
%        .label            - String label
%    seg_per_um        - Segments per micron for discretisation [seg/um]
%    L_0               - Characteristic length scale [m] (unused; retained
%                         for API parity with pili_building)
%    dis               - Gap between body surface and pili root [m]
%    E                 - Young's modulus [Pa]
%    r_pili            - Pili radius [m]
%
%  Outputs:
%    Pili - [1 x M] Structure array (same fields as pili_building output):
%        .N, .L, .ell, .x, .v, .kb, .kb0, .A, .r, .f_bend,
%        .attach_idx, .root_local, .root_dir, .type, .is_clamped
%
%    Semantic corrections (v7.1):
%      - .root_dir now stores the ACTUAL first-segment tangent after all
%        transformations (alignment + deflection), matching reconstruction
%        logic in init_coupled_params.  Previously set to site.direction,
%        which diverged from the true tangent whenever deflection != 0.
%
%  See also: pili_building, read_geometry_csv_mat, rotation_between_vectors,
%            normalise_json_struct_array
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-28
%  Updated: 2026-05-04 — Local-basis deflection, correct root_dir, uniform
%                         segment-length check, centerline deduplication,
%                         alignment-angle logging.

%% Input validation
if nargin < 7
    error('custom_pili_building:InsufficientInputs', ...
          'custom_pili_building requires 7 input arguments.');
end

%% Default bending stiffness
I_pili    = pi * r_pili^4 / 4;    % Second moment of area [m^4]
A_default = E * I_pili;            % Bending stiffness [N*m^2]

%% Normalise config input (handle JSON cell-array case)
custom_pili_cfg = normalise_json_struct_array(custom_pili_cfg);
NPili = numel(custom_pili_cfg);

%% Pre-allocate output structure array (same fields as pili_building)
template = struct('N', [], 'L', [], 'ell', [], ...
    'x', [], 'v', [], 'kb', [], 'kb0', [], ...
    'A', [], 'r', [], 'f_bend', [], ...
    'attach_idx', [], 'root_local', [], 'root_dir', [], ...
    'type', [], 'is_clamped', []);
if NPili == 0
    Pili = repmat(template, 1, 0);
    return;
end
Pili = repmat(template, 1, NPili);

%% Uniform-segment-length tolerance (relative)
UNIFORM_SEG_TOL = 1e-3;

%% Construct each pilus
for j = 1:NPili

    cfg_j    = custom_pili_cfg(j);
    site_idx = cfg_j.site_index;

    % Deflection angles (optional)
    deflection_z = 0;
    if isfield(cfg_j, 'deflection_z') && ~isempty(cfg_j.deflection_z)
        deflection_z = cfg_j.deflection_z;
    end
    deflection_x = 0;
    if isfield(cfg_j, 'deflection_x') && ~isempty(cfg_j.deflection_x)
        deflection_x = cfg_j.deflection_x;
    end

    % Bending stiffness (optional)
    if isfield(cfg_j, 'bending_stiffness') && ~isempty(cfg_j.bending_stiffness)
        A = cfg_j.bending_stiffness;
    else
        A = A_default;
    end

    % --- Read raw centerline from file ---
    raw_pts = read_geometry_csv_mat(cfg_j.points_file);   % [N_raw x 3] [m]

    % --- Remove consecutive duplicates (M5) ---
    raw_pts = remove_consecutive_duplicates(raw_pts, cfg_j.points_file, j);

    % --- Enforce uniform segment length (new requirement) ---
    seg_vecs = diff(raw_pts, 1, 1);
    seg_lens = sqrt(sum(seg_vecs.^2, 2));
    if any(seg_lens <= 0)
        error('custom_pili_building:ZeroSegment', ...
            'Pilus %d (%s): at least one segment has zero length even after deduplication.', ...
            j, cfg_j.points_file);
    end
    seg_len_mean = mean(seg_lens);
    seg_len_spread = (max(seg_lens) - min(seg_lens)) / seg_len_mean;
    if seg_len_spread > UNIFORM_SEG_TOL
        error('custom_pili_building:NonUniformSegments', ...
            ['Pilus %d (%s): segment lengths are not uniform ', ...
             '(min = %.4e m, max = %.4e m, mean = %.4e m, ', ...
             'relative spread = %.3e > tolerance %.1e). ', ...
             'User-supplied pili centerlines MUST be pre-discretized ', ...
             'with equal-length segments.'], ...
             j, cfg_j.points_file, ...
             min(seg_lens), max(seg_lens), seg_len_mean, ...
             seg_len_spread, UNIFORM_SEG_TOL);
    end

    L_actual = sum(seg_lens);

    % --- Determine target segment count (consistent with pili_building) ---
    L_um  = L_actual * 1e6;
    N     = max(ceil(L_um * seg_per_um) + 1, 3);
    ell   = L_actual / (N - 1);

    % --- Resample to uniform arc length (PCHIP) ---
    s_raw = [0; cumsum(seg_lens)];
    s_new = linspace(0, L_actual, N)';

    x_resampled = zeros(N, 3);
    for dim = 1:3
        x_resampled(:, dim) = interp1(s_raw, raw_pts(:, dim), s_new, 'pchip');
    end

    % --- Compute root tangent direction from resampled geometry ---
    root_tangent = x_resampled(2, :) - x_resampled(1, :);
    rt_norm      = norm(root_tangent);
    if rt_norm < eps
        error('custom_pili_building:DegenerateFirstSegment', ...
            'Pilus %d: first resampled segment has zero length.', j);
    end
    root_tangent = root_tangent / rt_norm;

    % --- Align root tangent to attachment site direction ---
    site     = attachment_sites(site_idx);
    site_dir = site.direction(:);

    % Rotation from current root tangent to site direction
    R_align = rotation_between_vectors(root_tangent.', site_dir);

    % Log alignment angle (M6)
    cos_theta_align = max(-1, min(1, dot(root_tangent.', site_dir)));
    theta_align_deg = acos(cos_theta_align) * 180 / pi;

    % Apply alignment around the origin (shape preserved; root moves)
    x_aligned = (R_align * x_resampled.').';

    % --- Build local basis for deflection (M2) ---
    %   e3 = site_dir                (aligned tangent)
    %   e1 = unit perpendicular via Gram-Schmidt (bending-plane normal)
    %   e2 = e3 x e1                 (binormal in local plane)
    e3 = site_dir;
    helper = [1; 0; 0];
    if abs(dot(helper, e3)) > 0.9
        helper = [0; 1; 0];
    end
    e1 = helper - (helper' * e3) * e3;
    e1 = e1 / norm(e1);
    e2 = cross(e3, e1);
    B  = [e1, e2, e3];    % lab <- local basis matrix (columns = e1, e2, e3)

    % Deflection rotations expressed in the LOCAL basis (e1, e2, e3):
    %   deflection_z : rotation about e3 (twist around aligned tangent)
    %   deflection_x : rotation about e1 (tilt in the bending plane)
    %
    % NOTE ON SEMANTICS:
    %   The axis-bindings above are NOT the same as in parametric
    %   pili_type_building.m.  There the pilus is built along +X and the
    %   matrices rotate about lab-frame Z (RZ, named deflection_z) and
    %   lab-frame X (RX, named deflection_x).  Because X IS the initial
    %   tangent in the parametric construction, parametric deflection_x
    %   is a TWIST around the tangent while parametric deflection_z is a
    %   TILT.  In the custom local-basis convention used here, the axes
    %   are swapped (deflection_z = twist, deflection_x = tilt).  The
    %   parameter NAMES are kept for JSON-compatibility only; users
    %   migrating a deflection setting from parametric to custom must
    %   swap deflection_x <-> deflection_z.  RX_loc retains the parametric
    %   left-handed sign for deflection_x.
    RZ_loc = [ cos(deflection_z), -sin(deflection_z), 0;
               sin(deflection_z),  cos(deflection_z), 0;
               0,                   0,                  1];
    RX_loc = [1,  0,                  0;
              0,  cos(deflection_x),  sin(deflection_x);
              0, -sin(deflection_x),  cos(deflection_x)];
    R_deflect_local = RX_loc * RZ_loc;
    R_deflect       = B * R_deflect_local * B.';     % local -> lab conjugation

    % Apply deflection about the first-node location (node 1 stays put)
    pivot     = x_aligned(1, :);
    x_deflect = (R_deflect * (x_aligned - pivot).').' + pivot;

    % --- Translate to (site.position + dis * site.direction) ---
    root_pos = site.position(:) + dis * site.direction(:);
    offset   = root_pos.' - x_deflect(1, :);
    x_final  = x_deflect + repmat(offset, N, 1);

    % --- Compute true first-segment tangent AFTER all transformations (M3) ---
    root_tangent_final = x_final(2, :) - x_final(1, :);
    rtf_norm           = norm(root_tangent_final);
    if rtf_norm < eps
        error('custom_pili_building:DegeneratePostTransform', ...
            'Pilus %d: first segment collapsed after transformations.', j);
    end
    root_tangent_final = (root_tangent_final / rtf_norm).';   % [3 x 1]

    % --- Pack as [3 x N] column-major format (same as pili_building) ---
    x      = x_final.';
    v      = zeros(3, N);
    f_bend = zeros(3, N);

    % --- Store in structure ---
    Pili(j).N          = N;
    Pili(j).L          = L_actual;
    Pili(j).ell        = ell;
    Pili(j).x          = x;
    Pili(j).v          = v;
    Pili(j).kb         = zeros(3, 0);   % Recomputed by init_coupled_params
    Pili(j).kb0        = zeros(3, 0);   % Recomputed by init_coupled_params
    Pili(j).A          = A;
    Pili(j).r          = r_pili;
    Pili(j).f_bend     = f_bend;
    Pili(j).attach_idx = site_idx;
    Pili(j).root_local = root_pos;                    % [3 x 1]
    Pili(j).root_dir   = root_tangent_final;          % [3 x 1], ACTUAL tangent
    Pili(j).type       = 0;                           % 0 = custom
    Pili(j).is_clamped = true;

    % --- Diagnostic logging ---
    fprintf(['  [Custom Pili %d] L = %.3f um, N = %d, site = %d (''%s''), ', ...
             'defl = [%.3f, %.3f] rad, align_angle = %.2f deg\n'], ...
        j, L_um, N, site_idx, site.label, deflection_z, deflection_x, ...
        theta_align_deg);

    if theta_align_deg > 45.0
        warning('custom_pili_building:LargeAlignmentAngle', ...
            ['Pilus %d: root tangent was rotated by %.2f deg (>45 deg) ', ...
             'to match site direction. Verify the source centerline is ', ...
             'pre-oriented consistently with attachment_sites.direction.'], ...
            j, theta_align_deg);
    end

end

end


%% ==================== Local Helper Functions ====================

function pts = remove_consecutive_duplicates(pts, fname, idx)
%REMOVE_CONSECUTIVE_DUPLICATES  Drop rows equal to their immediate predecessor.
%
%  Tolerance is scaled by the bounding-box extent (same strategy as
%  custom_body_building).  Throws an error if fewer than 2 distinct
%  points remain.

    if size(pts, 1) < 2
        error('custom_pili_building:TooFewPoints', ...
            'Pilus %d (%s): need at least 2 points, got %d.', ...
            idx, fname, size(pts, 1));
    end

    tol = 1e-6 * max(max(pts) - min(pts));
    tol = max(tol, 1e-15);

    diffs    = sqrt(sum(diff(pts, 1, 1).^2, 2));
    keep     = [true; diffs > tol];
    n_before = size(pts, 1);
    pts      = pts(keep, :);
    n_after  = size(pts, 1);

    if n_after < n_before
        warning('custom_pili_building:DuplicateCenterlinePoints', ...
            ['Pilus %d (%s): %d consecutive duplicate centerline point(s) ', ...
             'removed (tol = %.3e m).'], ...
            idx, fname, n_before - n_after, tol);
    end

    if n_after < 2
        error('custom_pili_building:InsufficientPointsAfterDedup', ...
            'Pilus %d (%s): fewer than 2 distinct points after deduplication.', ...
            idx, fname);
    end

end
