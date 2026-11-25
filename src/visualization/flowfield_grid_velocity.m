function planes = flowfield_grid_velocity(hydro, params, flow_opts, logf, silent)
%FLOWFIELD_GRID_VELOCITY  Evaluate XZ, YZ, and XY Stokes velocity planes
%
%  Description:
%    Builds XZ, YZ, and XY grids (through the body centre) that cover the
%    configured window and the whole cell (body and flagella), and evaluates
%    the regularised Stokeslet, Stokeslet-segment, and rotlet kernels
%    outside the body.  The body interior is left blank (NaN).  The ambient
%    flow is added when hydro.bg_fun is present.
%
%    All three planes are evaluated in a single batched pass over exterior
%    grid points so body and pilus source arrays are uploaded to the GPU
%    (or prepared on CPU) once per frame.  FP64 precision is preserved
%    throughout.
%
%  Syntax:
%    planes = flowfield_grid_velocity(hydro, params, flow_opts, logf, silent)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28
%  Updated: 2026-10-30
%  Updated: 2026-11-11 — added the XY plane (Z = Z_cm)
%  Updated: 2026-11-13 — optional streamline spacing and length overrides
%  Updated: 2026-11-24 — batched 3-plane FP64 GPU/BLAS3 evaluation,
%           analytic capsule section outlines, 2x2 surface corner grids

if nargin < 4 || isempty(logf), logf = @(varargin) []; end
if nargin < 5 || isempty(silent), silent = false; end
grid_spacing_um = get_flow_opt(flow_opts, 'grid_spacing_um', 0.05);
grid_range_um   = get_flow_opt(flow_opts, 'grid_range_um',   10.0);
quiver_skip     = get_flow_opt(flow_opts, 'quiver_skip',     5);
chunk_size      = get_flow_opt(flow_opts, 'chunk_size',      2000);
[stream_spacing_um, stream_length_um] = flowfield_stream_overrides(flow_opts);
if ~isempty(stream_spacing_um) || ~isempty(stream_length_um)
    if isempty(stream_spacing_um)
        sp_txt = 'default';
    else
        sp_txt = sprintf('%.4g um', stream_spacing_um);
    end
    if isempty(stream_length_um)
        len_txt = 'default';
    else
        len_txt = sprintf('%.4g um', stream_length_um);
    end
    logf('[FLOW] Streamlines: spacing %s, length %s\n', sp_txt, len_txt);
end

L_to_um     = hydro.L_to_um;
um_to_L     = hydro.um_to_L;
um_s_factor = hydro.um_s_factor;
X_cm        = hydro.X_cm;
R_body      = hydro.R_body;
X_cm_um     = X_cm * L_to_um;
if isfield(flow_opts, 'body_vis') && isstruct(flow_opts.body_vis) && ~isempty(flow_opts.body_vis)
    body_vis = flow_opts.body_vis;
else
    body_vis = build_body_visualization(params, hydro.L_0);
end
gpu_active = resolve_gpu_flag(flow_opts);

if isfield(flow_opts, 'fixed_below_um') && ~isempty(flow_opts.fixed_below_um)
    % Same span on every frame; the window translates with the centre of mass.
    % The step count is taken from this span, not from (hi - lo): adding and
    % subtracting the centre of mass makes (hi - lo) flicker by an ulp, and
    % ceil then adds a whole grid step on some frames.
    below = flow_opts.fixed_below_um(:);
    above = flow_opts.fixed_above_um(:);
    n_fixed = max(1, ceil((above + below) / grid_spacing_um - 1e-9));
    lo_um = (X_cm_um(:) - below).';
    hi_um = lo_um + (n_fixed.' * grid_spacing_um);
else
    n_fixed = [];
    [lo_um, hi_um] = flowfield_frame_window(hydro.geom, params, flow_opts);
end

logf(['[FLOW] Grid X [%.2f, %.2f], Y [%.2f, %.2f], Z [%.2f, %.2f] um, ', ...
    'spacing = %.2f um, quiver skip = %d\n'], ...
    lo_um(1), hi_um(1), lo_um(2), hi_um(2), lo_um(3), hi_um(3), ...
    grid_spacing_um, quiver_skip);

x_vec = grid_samples(lo_um(1), hi_um(1), grid_spacing_um, n_fixed, 1);
y_vec = grid_samples(lo_um(2), hi_um(2), grid_spacing_um, n_fixed, 2);
z_vec = grid_samples(lo_um(3), hi_um(3), grid_spacing_um, n_fixed, 3);

nx = numel(x_vec);
ny = numel(y_vec);
nz = numel(z_vec);

% ---- Build dimensionless plane coordinates and body-interior masks ----
pts_xz_dl = [repelem(x_vec(:), nz) * um_to_L, ...
             repmat(X_cm(2), nx * nz, 1), ...
             repmat(z_vec(:), nx, 1) * um_to_L];
[mask_xz_vec, ~] = body_inside_mask(pts_xz_dl, X_cm, R_body, body_vis);
mask_xz = reshape(mask_xz_vec, [nz, nx]);
out_xz  = find(~mask_xz_vec);

pts_yz_dl = [repmat(X_cm(1), ny * nz, 1), ...
             repelem(y_vec(:), nz) * um_to_L, ...
             repmat(z_vec(:), ny, 1) * um_to_L];
[mask_yz_vec, ~] = body_inside_mask(pts_yz_dl, X_cm, R_body, body_vis);
mask_yz = reshape(mask_yz_vec, [nz, ny]);
out_yz  = find(~mask_yz_vec);

pts_xy_dl = [repelem(x_vec(:), ny) * um_to_L, ...
             repmat(y_vec(:), nx, 1) * um_to_L, ...
             repmat(X_cm(3), nx * ny, 1)];
[mask_xy_vec, ~] = body_inside_mask(pts_xy_dl, X_cm, R_body, body_vis);
mask_xy = reshape(mask_xy_vec, [ny, nx]);
out_xy  = find(~mask_xy_vec);

% ---- Single batched evaluation across all 3 planes ----
n_xz = numel(out_xz);
n_yz = numel(out_yz);
X_eval_all = [pts_xz_dl(out_xz, :); pts_yz_dl(out_yz, :); pts_xy_dl(out_xy, :)];

prog_tag = '[FLOW] ';
if silent, prog_tag = ''; end
logf('[FLOW] Computing XZ, YZ, and XY flow fields (%d exterior points)...\n', ...
    size(X_eval_all, 1));
u_all = evaluate_points_velocity(X_eval_all, X_cm, ...
    hydro.body_data, hydro.seg_data, um_s_factor, chunk_size, ...
    prog_tag, hydro.bg_fun, gpu_active);

% ---- Scatter batched velocities back to the three planes ----
i1 = n_xz;
i2 = n_xz + n_yz;
[Ux_xz, Uy_xz, Uz_xz, Umag_xz] = scatter_plane_velocity( ...
    u_all(1:i1, :), out_xz, [nz, nx]);
[Ux_yz, Uy_yz, Uz_yz, Umag_yz] = scatter_plane_velocity( ...
    u_all(i1 + 1:i2, :), out_yz, [nz, ny]);
[Ux_xy, Uy_xy, Uz_xy, Umag_xy] = scatter_plane_velocity( ...
    u_all(i2 + 1:end, :), out_xy, [ny, nx]);

logf('[FLOW] XZ velocity range: [%.4e, %.4e] um/s\n', ...
    min(Umag_xz(:), [], 'omitnan'), max(Umag_xz(:), [], 'omitnan'));
logf('[FLOW] YZ velocity range: [%.4e, %.4e] um/s\n', ...
    min(Umag_yz(:), [], 'omitnan'), max(Umag_yz(:), [], 'omitnan'));
logf('[FLOW] XY velocity range: [%.4e, %.4e] um/s\n', ...
    min(Umag_xy(:), [], 'omitnan'), max(Umag_xy(:), [], 'omitnan'));

[cap_xz_h, cap_xz_v] = body_section_outline(body_vis, R_body, X_cm_um, 'xz');
[cap_yz_h, cap_yz_v] = body_section_outline(body_vis, R_body, X_cm_um, 'yz');
[cap_xy_h, cap_xy_v] = body_section_outline(body_vis, R_body, X_cm_um, 'xy');

% 2x2 corner coordinates suffice for texture-mapped 3D surfaces
HG_xz = [x_vec(1), x_vec(end); x_vec(1), x_vec(end)];
VG_xz = [X_cm_um(2), X_cm_um(2); X_cm_um(2), X_cm_um(2)];
ZG_xz = [z_vec(1), z_vec(1); z_vec(end), z_vec(end)];

HG_yz = [X_cm_um(1), X_cm_um(1); X_cm_um(1), X_cm_um(1)];
VG_yz = [y_vec(1), y_vec(end); y_vec(1), y_vec(end)];
ZG_yz = [z_vec(1), z_vec(1); z_vec(end), z_vec(end)];

HG_xy = [x_vec(1), x_vec(end); x_vec(1), x_vec(end)];
VG_xy = [y_vec(1), y_vec(1); y_vec(end), y_vec(end)];
ZG_xy = [X_cm_um(3), X_cm_um(3); X_cm_um(3), X_cm_um(3)];

planes.body_vis = body_vis;
planes.axis_lim_um = [x_vec(1), x_vec(end); ...
                      y_vec(1), y_vec(end); ...
                      z_vec(1), z_vec(end)];
planes.half_range = max(abs(planes.axis_lim_um - X_cm_um), [], 'all');
planes.grid_spacing_um = grid_spacing_um;
planes.grid_range_um = grid_range_um;
planes.quiver_skip = quiver_skip;
planes.stream_spacing_um = stream_spacing_um;
planes.stream_length_um = stream_length_um;
planes.X_cm_um = X_cm_um;
planes.xz = pack_plane(x_vec, z_vec, HG_xz, VG_xz, ZG_xz, ...
    Ux_xz, Uy_xz, Uz_xz, Umag_xz, mask_xz, cap_xz_h, cap_xz_v);
planes.yz = pack_plane(y_vec, z_vec, HG_yz, VG_yz, ZG_yz, ...
    Ux_yz, Uy_yz, Uz_yz, Umag_yz, mask_yz, cap_yz_h, cap_yz_v);
planes.xy = pack_plane(x_vec, y_vec, HG_xy, VG_xy, ZG_xy, ...
    Ux_xy, Uy_xy, Uz_xy, Umag_xy, mask_xy, cap_xy_h, cap_xy_v);
samples = sqrt(sum(u_all.^2, 2));
planes.samples = samples(isfinite(samples));

end


function plane = pack_plane(h_vec, v_vec, HG, VG, ZG, Ux, Uy, Uz, Umag, mask, cap_h, cap_v)
plane.h_vec = h_vec;
plane.v_vec = v_vec;
plane.HG = HG;
plane.VG = VG;
plane.ZG = ZG;
plane.Ux = Ux;
plane.Uy = Uy;
plane.Uz = Uz;
plane.Umag = Umag;
plane.mask = mask;
plane.cap_h = cap_h;
plane.cap_v = cap_v;
end


function [Ux, Uy, Uz, Umag] = scatter_plane_velocity(u_ev, outside_idx, grid_size)
Ux = NaN(grid_size);
Uy = NaN(grid_size);
Uz = NaN(grid_size);
Ux(outside_idx) = u_ev(:, 1);
Uy(outside_idx) = u_ev(:, 2);
Uz(outside_idx) = u_ev(:, 3);
Umag = sqrt(Ux.^2 + Uy.^2 + Uz.^2);
end


%% =====================================================================
%  LOCAL FUNCTION: resolve_gpu_flag
%  =====================================================================
function gpu_active = resolve_gpu_flag(flow_opts)
%RESOLVE_GPU_FLAG  Resolve 'auto' / true / false GPU mode with cached probe
    req = 'auto';
    if isfield(flow_opts, 'use_gpu') && ~isempty(flow_opts.use_gpu)
        req = flow_opts.use_gpu;
    end
    if islogical(req) || isnumeric(req)
        want_gpu = logical(req(1));
    elseif ischar(req) || isstring(req)
        want_gpu = ~strcmpi(char(req), 'false') && ~strcmpi(char(req), 'off');
    else
        want_gpu = false;
    end
    if ~want_gpu
        gpu_active = false;
        return;
    end
    gpu_active = has_usable_gpu();
end


function tf = has_usable_gpu()
    persistent cached_tf
    if ~isempty(cached_tf)
        tf = cached_tf;
        return;
    end
    tf = false;
    try
        if exist('gpuDeviceCount', 'file') == 2 && gpuDeviceCount > 0
            tf = true;
        end
    catch
        tf = false;
    end
    cached_tf = tf;
end


%% =====================================================================
%  LOCAL FUNCTION: evaluate_points_velocity
%  =====================================================================
function u_um = evaluate_points_velocity(X_eval, X_cm, ...
    body_data, seg_data, um_s_factor, chunk_size, progress_tag, bg_fun, gpu_active)
%EVALUATE_POINTS_VELOCITY  Batched FP64 RSM + RSS + Rotlet velocity kernel
%
%  Evaluates the disturbance velocity at all N_eval exterior points in
%  body-CM-shifted coordinates (eliminating large-translation cancellation),
%  using BLAS3 matrix-matrix contractions and optional GPU execution.
%  Adds the ambient flow bg_fun(X_eval) in unshifted lab coordinates.

    N_eval = size(X_eval, 1);
    if N_eval == 0
        u_um = zeros(0, 3);
        return;
    end

    if gpu_active
        try
            u_dl = evaluate_points_kernel(X_eval, X_cm, body_data, seg_data, ...
                chunk_size, progress_tag, true);
        catch
            u_dl = evaluate_points_kernel(X_eval, X_cm, body_data, seg_data, ...
                chunk_size, progress_tag, false);
        end
    else
        u_dl = evaluate_points_kernel(X_eval, X_cm, body_data, seg_data, ...
            chunk_size, progress_tag, false);
    end

    if ~isempty(bg_fun)
        u_dl = u_dl + bg_fun(X_eval);
    end
    u_um = u_dl * um_s_factor;
end


function u_ev = evaluate_points_kernel(X_eval, X_cm, body_data, seg_data, ...
    chunk_size, progress_tag, gpu_active)
    N_eval = size(X_eval, 1);
    show_progress = ~isempty(progress_tag);

    % Shift coordinates relative to body CM so O(1) values keep BLAS3
    % contractions free of catastrophic cancellation after long drift.
    cm_row = X_cm(:).';
    X_ev_rel = X_eval - cm_row;
    X_b_rel  = body_data.X - cm_row;
    bf_mat   = [body_data.fx(:), body_data.fy(:), body_data.fz(:)];
    e2_b     = body_data.e2;
    n_body   = size(X_b_rel, 1);

    has_pili = seg_data.active;
    n_seg    = 0;
    if has_pili
        n_seg  = numel(seg_data.L);
        e2_p   = seg_data.e2;
        L_r    = seg_data.L.';
        L2_r   = seg_data.L2.';
        inv_L  = 1 ./ L_r;
        inv_L2 = 1 ./ L2_r;
        inv3L2 = inv_L2 / 3;
        v_mat  = seg_data.v;
        v1_r   = v_mat(:, 1).';
        v2_r   = v_mat(:, 2).';
        v3_r   = v_mat(:, 3).';
        y0_rel = seg_data.y0 - cm_row;
        y0_1   = y0_rel(:, 1).';
        y0_2   = y0_rel(:, 2).';
        y0_3   = y0_rel(:, 3).';
        f0_mat = seg_data.f0;
        f1_mat = seg_data.f1;
        f0x = f0_mat(:, 1).';  f0y = f0_mat(:, 2).';  f0z = f0_mat(:, 3).';
        f1x = f1_mat(:, 1).';  f1y = f1_mat(:, 2).';  f1z = f1_mat(:, 3).';
        v_dot_f0 = v1_r .* f0x + v2_r .* f0y + v3_r .* f0z;
        v_dot_f1 = v1_r .* f1x + v2_r .* f1y + v3_r .* f1z;

        has_rot = seg_data.has_rotlet;
        if has_rot
            tx = seg_data.tau_x(:);
            ty = seg_data.tau_y(:);
            tz = seg_data.tau_z(:);
            tau_mat = [tx, ty, tz];
            % tau x y0_rel and tau x v_mat for rotlet BLAS3 contraction
            t_cross_y0 = [y0_rel(:, 3) .* ty - y0_rel(:, 2) .* tz, ...
                          y0_rel(:, 1) .* tz - y0_rel(:, 3) .* tx, ...
                          y0_rel(:, 2) .* tx - y0_rel(:, 1) .* ty];
            t_cross_v  = [v_mat(:, 3) .* ty - v_mat(:, 2) .* tz, ...
                          v_mat(:, 1) .* tz - v_mat(:, 3) .* tx, ...
                          v_mat(:, 2) .* tx - v_mat(:, 1) .* ty];
        end
    end

    % Use larger chunks on GPU (bounded to stay safely inside VRAM)
    if gpu_active
        max_src = max([n_body, n_seg, 1]);
        eff_chunk = min(max(chunk_size, 16384), max(2048, floor(8e6 / max_src)));
    else
        eff_chunk = chunk_size;
    end

    % Upload source arrays to GPU once for all three planes
    X_bx = X_b_rel(:, 1).';
    X_by = X_b_rel(:, 2).';
    X_bz = X_b_rel(:, 3).';
    bfx_r = bf_mat(:, 1).';
    bfy_r = bf_mat(:, 2).';
    bfz_r = bf_mat(:, 3).';

    if gpu_active
        X_bx = gpuArray(X_bx);     X_by = gpuArray(X_by);     X_bz = gpuArray(X_bz);
        X_b_rel = gpuArray(X_b_rel);
        bf_mat  = gpuArray(bf_mat);
        bfx_r = gpuArray(bfx_r);   bfy_r = gpuArray(bfy_r);   bfz_r = gpuArray(bfz_r);
        if has_pili
            L_r    = gpuArray(L_r);     L2_r   = gpuArray(L2_r);
            inv_L  = gpuArray(inv_L);   inv_L2 = gpuArray(inv_L2);  inv3L2 = gpuArray(inv3L2);
            v_mat  = gpuArray(v_mat);
            v1_r   = gpuArray(v1_r);    v2_r   = gpuArray(v2_r);    v3_r   = gpuArray(v3_r);
            y0_rel = gpuArray(y0_rel);
            y0_1   = gpuArray(y0_1);    y0_2   = gpuArray(y0_2);    y0_3   = gpuArray(y0_3);
            f0_mat = gpuArray(f0_mat);  f1_mat = gpuArray(f1_mat);
            f0x    = gpuArray(f0x);     f0y    = gpuArray(f0y);     f0z    = gpuArray(f0z);
            f1x    = gpuArray(f1x);     f1y    = gpuArray(f1y);     f1z    = gpuArray(f1z);
            v_dot_f0 = gpuArray(v_dot_f0);
            v_dot_f1 = gpuArray(v_dot_f1);
            if has_rot
                tau_mat    = gpuArray(tau_mat);
                t_cross_y0 = gpuArray(t_cross_y0);
                t_cross_v  = gpuArray(t_cross_v);
            end
        end
    end

    u_ev = zeros(N_eval, 3);
    n_chunks = ceil(N_eval / eff_chunk);
    if show_progress
        t_bar = tic;
    end

    chunk_idx = 0;
    for ci = 1:eff_chunk:N_eval
        chunk_idx = chunk_idx + 1;
        ce  = min(ci + eff_chunk - 1, N_eval);
        idx = ci:ce;

        if gpu_active
            X_c = gpuArray(X_ev_rel(idx, :));
        else
            X_c = X_ev_rel(idx, :);
        end
        x1 = X_c(:, 1);
        x2 = X_c(:, 2);
        x3 = X_c(:, 3);

        % ---- RSM body contribution (BLAS3 contraction) ----
        rx = x1 - X_bx;
        ry = x2 - X_by;
        rz = x3 - X_bz;
        r2 = rx.^2 + ry.^2 + rz.^2;
        inv_R = 1 ./ sqrt(r2 + e2_b);
        H2    = inv_R .* inv_R .* inv_R;
        H1    = (r2 + 2 * e2_b) .* H2;
        H2_rd = H2 .* (rx .* bfx_r + ry .* bfy_r + rz .* bfz_r);
        u_c   = H1 * bf_mat + X_c .* sum(H2_rd, 2) - H2_rd * X_b_rel;

        % ---- RSS pili + Rotlet contribution ----
        if has_pili
            x01 = x1 - y0_1;
            x02 = x2 - y0_2;
            x03 = x3 - y0_3;
            r0_sq = x01.^2 + x02.^2 + x03.^2;

            x11 = x01 + v1_r;  x12 = x02 + v2_r;  x13 = x03 + v3_r;
            r1_sq = x11.^2 + x12.^2 + x13.^2;

            R0 = sqrt(r0_sq + e2_p);
            R1 = sqrt(r1_sq + e2_p);
            inv_R0 = 1 ./ R0;
            inv_R1 = 1 ./ R1;
            d0 = x01 .* v1_r + x02 .* v2_r + x03 .* v3_r;

            H0_v = max(L_r .* R0 + d0,        1e-300);
            H1_v = max(L_r .* R1 + d0 + L2_r, 1e-300);
            d0_L2 = d0 .* inv_L2;

            T0m1 = log(H1_v ./ H0_v) .* inv_L;
            T0m3 = inv_R0 ./ H0_v - inv_R1 ./ H1_v;
            T1m1 = (R1 - R0) .* inv_L2 - d0_L2 .* T0m1;
            T1m3 = (inv_R0 - inv_R1) .* inv_L2 - d0_L2 .* T0m3;
            neg_inv_R1_L2 = -inv_R1 .* inv_L2;
            T2m3 = neg_inv_R1_L2 + T0m1 .* inv_L2 - d0_L2 .* T1m3;
            T3m3 = neg_inv_R1_L2 + 2 * T1m1 .* inv_L2 - d0_L2 .* T2m3;

            ci2 = L_r .* (T1m1 + e2_p * T1m3);
            c12 = L_r .* T1m3;
            c22 = L_r .* T2m3;
            c32 = L_r .* T3m3;

            ci1 = L_r .* (T0m1 + e2_p * T0m3) - ci2;
            c11 = L_r .* T0m3 - c12;
            c21 = c12 - c22;
            c31 = c22 - c32;

            x0_f0 = x01 .* f0x + x02 .* f0y + x03 .* f0z;
            x0_f1 = x01 .* f1x + x02 .* f1y + x03 .* f1z;

            a_x0 = (c11 .* x0_f0 + c21 .* v_dot_f0) + (c12 .* x0_f1 + c22 .* v_dot_f1);
            a_v  = (c21 .* x0_f0 + c31 .* v_dot_f0) + (c22 .* x0_f1 + c32 .* v_dot_f1);

            u_c = u_c + ci1 * f0_mat + ci2 * f1_mat ...
                + X_c .* sum(a_x0, 2) - a_x0 * y0_rel + a_v * v_mat;

            if has_rot
                inv_R0_3   = inv_R0 .* inv_R0 .* inv_R0;
                inv_R1_3   = inv_R1 .* inv_R1 .* inv_R1;
                d_sq       = max(L2_r .* (r0_sq + e2_p) - d0.^2, e2_p .* L2_r);
                inv_d_sq_3 = (1 / 3) ./ d_sq;
                d1_v       = d0 + L2_r;

                T0m5 = inv_d_sq_3 .* (d1_v .* inv_R1_3 - d0 .* inv_R0_3 + 2 * L2_r .* T0m3);
                T1m5 = (inv_R0_3 - inv_R1_3) .* inv3L2 - d0_L2 .* T0m5;

                Lca = L_r .* (2 * T0m3 + 3 * e2_p * T0m5);
                Lcb = L_r .* (2 * T1m3 + 3 * e2_p * T1m5);
                T_ca = Lca * tau_mat;
                u_rot_x = [x3 .* T_ca(:, 2) - x2 .* T_ca(:, 3), ...
                           x1 .* T_ca(:, 3) - x3 .* T_ca(:, 1), ...
                           x2 .* T_ca(:, 1) - x1 .* T_ca(:, 2)];
                u_c = u_c + 0.5 * (u_rot_x - Lca * t_cross_y0 + Lcb * t_cross_v);
            end
        end

        if gpu_active
            u_ev(idx, :) = gather(u_c);
        else
            u_ev(idx, :) = u_c;
        end

        if show_progress
            render_progress_bar(progress_tag, chunk_idx, n_chunks, t_bar);
        end
    end
end


%% =====================================================================
%  LOCAL FUNCTION: build_body_visualization
%  =====================================================================
function body_vis = build_body_visualization(params, L_0)
%BUILD_BODY_VISUALIZATION  Mode-aware body descriptor for flow-field plot.
%
%  Returns a struct with the data that evaluate_points_velocity,
%  build_body_mesh, and body_section_outline need to render the body and
%  mask its interior.  For parametric mode, stores capsule SDF params.
%  For custom mode, stores dimensionless vertices/faces and an alphaShape
%  reconstructed in the body-local frame for inShape() masking.

    body_vis = struct('type', '');

    body_type = 'capsule';
    if isfield(params, 'body_type') && ~isempty(params.body_type)
        body_type = params.body_type;
    end
    body_vis.type = body_type;

    switch body_type
        case 'capsule'
            body_vis.b_dl = params.dim.minor_axis / L_0;
            body_vis.c_dl = params.dim.major_axis / L_0 - body_vis.b_dl;
            body_vis.a_um = params.dim.major_axis * 1e6;
            body_vis.b_um = params.dim.minor_axis * 1e6;

        case 'custom'
            if ~isfield(params, 'body_verts_local') || isempty(params.body_verts_local)
                error('compute_plot_flowfield:MissingCustomVerts', ...
                    ['Custom geometry mode expects params.body_verts_local ', ...
                     '(set in init_coupled_params).']);
            end
            body_vis.verts_local_dl = params.body_verts_local;     % [L_0]
            body_vis.faces          = params.body_faces;           % index triangulation
            body_vis.verts_m        = params.body_verts;           % [m]

            % A critical alpha that merely connects the surface points is a
            % thin crust: inShape is then false in the cell interior and the
            % flow field is drawn inside the body.  Grow alpha until the
            % shape is a solid volume that contains the centroid.
            body_vis.alpha_shape_dl = [];
            try
                body_vis.alpha_shape_dl = solid_alpha_shape(body_vis.verts_local_dl);
            catch
                body_vis.alpha_shape_dl = [];
            end

            % Cache a bounding sphere for fallback masking
            centroid = mean(body_vis.verts_local_dl, 1);
            body_vis.bound_center_dl = centroid(:);
            body_vis.bound_radius_dl = max(sqrt(sum( ...
                (body_vis.verts_local_dl - centroid).^2, 2)));

        otherwise
            error('compute_plot_flowfield:UnknownBodyType', ...
                'Unknown params.body_type = ''%s''.', body_type);
    end
end


%% =====================================================================
%  LOCAL FUNCTION: body_inside_mask
%  =====================================================================
function [inside_mask_vec, dist_field] = body_inside_mask( ...
    grid_pts_dl, X_cm, R_body, body_vis)
%BODY_INSIDE_MASK  Return interior mask + SDF-like field for outline plotting.

    % Transform grid to body-local dimensionless frame
    X_local = (grid_pts_dl - X_cm(:).') * R_body;

    switch body_vis.type
        case 'capsule'
            x_cl = max(-body_vis.c_dl, min(body_vis.c_dl, X_local(:, 1)));
            d_cap = sqrt((X_local(:, 1) - x_cl).^2 + ...
                          X_local(:, 2).^2 + X_local(:, 3).^2);
            inside_mask_vec = d_cap <= body_vis.b_dl;
            dist_field      = d_cap - body_vis.b_dl;     % SDF: < 0 interior

        case 'custom'
            if ~isempty(body_vis.alpha_shape_dl)
                try
                    inside_mask_vec = inShape(body_vis.alpha_shape_dl, ...
                        X_local(:, 1), X_local(:, 2), X_local(:, 3));
                catch
                    inside_mask_vec = bound_sphere_inside(X_local, body_vis);
                end
            else
                inside_mask_vec = bound_sphere_inside(X_local, body_vis);
            end
            dist_field = -(double(inside_mask_vec) - 0.5);

        otherwise
            error('compute_plot_flowfield:UnknownBodyType', ...
                'Unknown body_vis.type = ''%s''.', body_vis.type);
    end

    inside_mask_vec = logical(inside_mask_vec);
end


function mask_vec = bound_sphere_inside(X_local, body_vis)
    d = sqrt(sum((X_local - body_vis.bound_center_dl.').^2, 2));
    mask_vec = d <= body_vis.bound_radius_dl;
end


%% =====================================================================
%  LOCAL FUNCTION: body_section_outline
%  =====================================================================
function [h, v] = body_section_outline(body_vis, R_body, X_cm_um, plane_name)
%BODY_SECTION_OUTLINE  Smooth body silhouette on the XZ, YZ, or XY plane
%
%  For parametric capsules, evaluates the exact analytical plane cut in
%  polar coordinates about X_cm_um (C1-continuous cylinder + spherical-cap
%  intersection).  For custom triangulated meshes, slices only triangles
%  that straddle the plane and resamples the resulting closed loop with a
%  periodic spline.

    if isempty(body_vis.type)
        error('flowfield_grid_velocity:BadBodyVis', 'body_vis.type is empty.');
    end

    if strcmp(plane_name, 'xz')
        fixed_axis = 2;  h_axis = 1;  v_axis = 3;
    elseif strcmp(plane_name, 'xy')
        fixed_axis = 3;  h_axis = 1;  v_axis = 2;
    else
        fixed_axis = 1;  h_axis = 2;  v_axis = 3;
    end

    if strcmp(body_vis.type, 'capsule')
        [h, v] = analytic_capsule_section(body_vis.a_um, body_vis.b_um, ...
            R_body, X_cm_um, h_axis, v_axis);
        return;
    end

    V = (R_body * (body_vis.verts_m * 1e6).' + X_cm_um(:)).';
    F = body_vis.faces;
    [h, v] = mesh_plane_loops(V, F, fixed_axis, X_cm_um(fixed_axis), h_axis, v_axis);
end


function [h, v] = analytic_capsule_section(a_um, b_um, R_body, X_cm_um, h_axis, v_axis)
%ANALYTIC_CAPSULE_SECTION  Exact central-plane cross-section of a capsule
    c_um = max(0, a_um - b_um);
    theta = linspace(0, 2 * pi, 241);
    theta = theta(1:end - 1);
    ct = cos(theta);
    st = sin(theta);

    eh_loc = R_body(h_axis, :).';
    ev_loc = R_body(v_axis, :).';
    d_loc  = eh_loc * ct + ev_loc * st;

    ax = abs(d_loc(1, :));
    s  = hypot(d_loc(2, :), d_loc(3, :));

    on_cyl = (b_um * ax <= c_um * s) & (s > 0);
    r = zeros(size(theta));
    r(on_cyl) = b_um ./ s(on_cyl);
    s_cap = s(~on_cyl);
    disc = max(b_um^2 - (c_um * s_cap).^2, 0);
    r(~on_cyl) = c_um * ax(~on_cyl) + sqrt(disc);

    h = [X_cm_um(h_axis) + r .* ct, nan];
    v = [X_cm_um(v_axis) + r .* st, nan];
end


function [h, v] = mesh_plane_loops(V, F, fixed_axis, fixed_value, h_axis, v_axis)
%MESH_PLANE_LOOPS  Smooth closed cuts of a triangle mesh by one plane
    s_all = V(:, fixed_axis) - fixed_value;
    sF = s_all(F);
    cross_tri = any(sF < 0, 2) & any(sF >= 0, 2);
    F_hit = F(cross_tri, :);
    nF = size(F_hit, 1);
    if nF == 0
        h = nan;
        v = nan;
        return
    end

    P = zeros(nF, 2);
    Q = zeros(nF, 2);
    nseg = 0;
    edges = [1, 2; 2, 3; 3, 1];
    for i = 1:nF
        tri = V(F_hit(i, :), :);
        s = tri(:, fixed_axis) - fixed_value;
        hits = zeros(2, 2);
        nh = 0;
        for e = 1:3
            a = edges(e, 1);
            b = edges(e, 2);
            if (s(a) < 0 && s(b) >= 0) || (s(b) < 0 && s(a) >= 0)
                denom = s(a) - s(b);
                if denom == 0
                    t = 0;
                else
                    t = s(a) / denom;
                end
                pt = tri(a, :) + t * (tri(b, :) - tri(a, :));
                nh = nh + 1;
                if nh <= 2
                    hits(nh, :) = pt([h_axis, v_axis]);
                end
            end
        end
        if nh >= 2 && norm(hits(1, :) - hits(2, :)) > 1e-9
            nseg = nseg + 1;
            P(nseg, :) = hits(1, :);
            Q(nseg, :) = hits(2, :);
        end
    end
    if nseg == 0
        h = nan;
        v = nan;
        return
    end
    P = P(1:nseg, :);
    Q = Q(1:nseg, :);
    scale = max(max(abs([P; Q]), [], 'all'), 1);
    tol = 1e-4 * scale;
    used = false(nseg, 1);
    h = [];
    v = [];
    for i = 1:nseg
        if used(i)
            continue
        end
        used(i) = true;
        loop = [P(i, :); Q(i, :)];
        for step = 1:nseg
            tail = loop(end, :);
            best = 0;
            best_d = tol;
            take_q = false;
            for j = 1:nseg
                if used(j)
                    continue
                end
                d1 = hypot(P(j, 1) - tail(1), P(j, 2) - tail(2));
                d2 = hypot(Q(j, 1) - tail(1), Q(j, 2) - tail(2));
                if d1 < best_d
                    best = j;
                    best_d = d1;
                    take_q = true;
                end
                if d2 < best_d
                    best = j;
                    best_d = d2;
                    take_q = false;
                end
            end
            if best == 0
                break
            end
            used(best) = true;
            if take_q
                loop = [loop; Q(best, :)]; %#ok<AGROW>
            else
                loop = [loop; P(best, :)]; %#ok<AGROW>
            end
        end
        if size(loop, 1) >= 4 && norm(loop(1, :) - loop(end, :)) <= tol
            loop = smooth_closed_loop(loop(1:end - 1, :));
            h = [h, loop(:, 1).', nan]; %#ok<AGROW>
            v = [v, loop(:, 2).', nan]; %#ok<AGROW>
        end
    end
    if isempty(h)
        h = nan;
        v = nan;
    end
end


function loop = smooth_closed_loop(loop)
%SMOOTH_CLOSED_LOOP  Periodic spline resample of one closed section
    n = size(loop, 1);
    if n < 4
        return
    end
    closed = [loop; loop(1, :)];
    seg = sqrt(sum(diff(closed, 1, 1).^2, 2));
    keep = seg > 1e-9;
    if ~all(keep)
        closed = closed([true; keep], :);
        if size(closed, 1) < 5
            return
        end
        if norm(closed(1, :) - closed(end, :)) > 1e-9
            closed = [closed; closed(1, :)];
        end
        seg = sqrt(sum(diff(closed, 1, 1).^2, 2));
    end
    s = [0; cumsum(seg)];
    if s(end) <= 0
        return
    end
    s = s / s(end);
    [su, ia] = unique(s, 'stable');
    closed = closed(ia, :);
    if numel(su) < 4
        return
    end
    n_out = max(160, 4 * n);
    sq = linspace(0, 1, n_out + 1).';
    sq = sq(1:end - 1);
    s3 = [su(1:end - 1) - 1; su; su(2:end) + 1];
    p3 = [closed(1:end - 1, :); closed; closed(2:end, :)];
    loop = [interp1(s3, p3(:, 1), sq, 'spline'), ...
            interp1(s3, p3(:, 2), sq, 'spline')];
end


function v = get_flow_opt(s, fname, default)
%GET_FLOW_OPT  Read field from options struct with fallback default
    if isfield(s, fname), v = s.(fname); else, v = default; end
end


function v = grid_samples(lo, hi, spacing, n_fixed, axis_index)
%GRID_SAMPLES  Uniform samples on one axis
%  n_fixed, when set, is the step count of the co-moving window.  It does
%  not depend on the centre of mass, so every frame has the same span.
    if ~(isfinite(lo) && isfinite(spacing) && spacing > 0)
        error('flowfield_grid_velocity:BadGrid', ...
            'Grid limits and spacing must be finite, and spacing must be positive.');
    end
    if ~isempty(n_fixed)
        n = n_fixed(axis_index);
    else
        if ~(isfinite(hi))
            error('flowfield_grid_velocity:BadGrid', ...
                'Grid limits and spacing must be finite, and spacing must be positive.');
        end
        if hi < lo
            tmp = lo;
            lo = hi;
            hi = tmp;
        end
        n = max(1, ceil((hi - lo) / spacing));
    end
    v = lo + (0:n) * spacing;
end


function shp = solid_alpha_shape(verts)
%SOLID_ALPHA_SHAPE  Alpha shape whose interior is the cell volume
%
%  The smallest one-region alpha is often only a crust on a surface point
%  cloud.  Increase it until the centroid is inside a volume that is a
%  substantial fraction of the bounding box.  The convex hull is the
%  fallback and is itself a solid set.

    shp = alphaShape(verts);
    alpha_crit = criticalAlpha(shp, 'one-region');
    if ~(isfinite(alpha_crit) && alpha_crit > 0)
        alpha_crit = criticalAlpha(shp, 'all-points');
    end
    centroid = mean(verts, 1);
    bb = max(verts, [], 1) - min(verts, [], 1);
    bb_vol = max(prod(bb), eps);
    chosen = Inf;
    if isfinite(alpha_crit) && alpha_crit > 0
        multipliers = logspace(0, log10(40), 48);
        for ai = 1:numel(multipliers)
            shp.Alpha = alpha_crit * multipliers(ai);
            if numRegions(shp) ~= 1
                continue
            end
            if ~inShape(shp, centroid(1), centroid(2), centroid(3))
                continue
            end
            vol = volume(shp);
            if isfinite(vol) && vol > 0.15 * bb_vol
                chosen = shp.Alpha;
                break
            end
        end
    end
    shp.Alpha = chosen;
end
