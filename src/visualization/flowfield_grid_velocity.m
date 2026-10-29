function planes = flowfield_grid_velocity(hydro, params, flow_opts, logf, silent)
%FLOWFIELD_GRID_VELOCITY  Evaluate XZ and YZ Stokes velocity planes
%
%  Description:
%    Builds the body-centred XZ and YZ grids and evaluates the regularised
%    Stokeslet, Stokeslet-segment, and rotlet kernels outside the body.
%    Interior samples are NaN.  The ambient flow is added when hydro.bg_fun
%    is present.
%
%  Syntax:
%    planes = flowfield_grid_velocity(hydro, params, flow_opts, logf, silent)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28

if nargin < 4 || isempty(logf), logf = @(varargin) []; end
if nargin < 5 || isempty(silent), silent = false; end
grid_spacing_um = get_flow_opt(flow_opts, 'grid_spacing_um', 0.05);
grid_range_um   = get_flow_opt(flow_opts, 'grid_range_um',   10.0);
quiver_skip     = get_flow_opt(flow_opts, 'quiver_skip',     5);
chunk_size      = get_flow_opt(flow_opts, 'chunk_size',      2000);

L_to_um     = hydro.L_to_um;
um_to_L     = hydro.um_to_L;
um_s_factor = hydro.um_s_factor;
X_cm        = hydro.X_cm;
R_body      = hydro.R_body;
X_cm_um     = X_cm * L_to_um;
half_range  = grid_range_um / 2;
body_vis    = build_body_visualization(params, hydro.L_0);

logf('[FLOW] Grid: %.1f x %.1f um, spacing = %.2f um, quiver skip = %d\n', ...
    grid_range_um, grid_range_um, grid_spacing_um, quiver_skip);

logf('[FLOW] Computing XZ-plane flow field...\n');
x_vec_xz = (X_cm_um(1) - half_range) : grid_spacing_um : (X_cm_um(1) + half_range);
z_vec_xz = (X_cm_um(3) - half_range) : grid_spacing_um : (X_cm_um(3) + half_range);
[XG_xz, ZG_xz] = meshgrid(x_vec_xz, z_vec_xz);
YG_xz = X_cm_um(2) * ones(size(XG_xz));
xz_tag = '[FLOW-XZ] ';
if silent, xz_tag = ''; end
[Ux_xz, Uy_xz, Uz_xz, Umag_xz, mask_xz, dist_xz] = evaluate_plane_velocity( ...
    [XG_xz(:)*um_to_L, YG_xz(:)*um_to_L, ZG_xz(:)*um_to_L], ...
    [numel(z_vec_xz), numel(x_vec_xz)], X_cm, R_body, body_vis, ...
    hydro.body_data, hydro.seg_data, um_s_factor, chunk_size, xz_tag, hydro.bg_fun);
logf('[FLOW] XZ velocity range: [%.4e, %.4e] um/s\n', ...
    min(Umag_xz(:), [], 'omitnan'), max(Umag_xz(:), [], 'omitnan'));

logf('[FLOW] Computing YZ-plane flow field...\n');
y_vec_yz = (X_cm_um(2) - half_range) : grid_spacing_um : (X_cm_um(2) + half_range);
z_vec_yz = (X_cm_um(3) - half_range) : grid_spacing_um : (X_cm_um(3) + half_range);
[YG_yz, ZG_yz] = meshgrid(y_vec_yz, z_vec_yz);
XG_yz = X_cm_um(1) * ones(size(YG_yz));
yz_tag = '[FLOW-YZ] ';
if silent, yz_tag = ''; end
[Ux_yz, Uy_yz, Uz_yz, Umag_yz, mask_yz, dist_yz] = evaluate_plane_velocity( ...
    [XG_yz(:)*um_to_L, YG_yz(:)*um_to_L, ZG_yz(:)*um_to_L], ...
    [numel(z_vec_yz), numel(y_vec_yz)], X_cm, R_body, body_vis, ...
    hydro.body_data, hydro.seg_data, um_s_factor, chunk_size, yz_tag, hydro.bg_fun);
logf('[FLOW] YZ velocity range: [%.4e, %.4e] um/s\n', ...
    min(Umag_yz(:), [], 'omitnan'), max(Umag_yz(:), [], 'omitnan'));

[cap_xz_h, cap_xz_v] = extract_body_outline( ...
    dist_xz, [numel(z_vec_xz), numel(x_vec_xz)], x_vec_xz, z_vec_xz, body_vis);
[cap_yz_h, cap_yz_v] = extract_body_outline( ...
    dist_yz, [numel(z_vec_yz), numel(y_vec_yz)], y_vec_yz, z_vec_yz, body_vis);

planes.body_vis = body_vis;
planes.half_range = half_range;
planes.grid_spacing_um = grid_spacing_um;
planes.grid_range_um = grid_range_um;
planes.quiver_skip = quiver_skip;
planes.X_cm_um = X_cm_um;
planes.xz = pack_plane(x_vec_xz, z_vec_xz, XG_xz, YG_xz, ZG_xz, ...
    Ux_xz, Uy_xz, Uz_xz, Umag_xz, mask_xz, cap_xz_h, cap_xz_v);
planes.yz = pack_plane(y_vec_yz, z_vec_yz, XG_yz, YG_yz, ZG_yz, ...
    Ux_yz, Uy_yz, Uz_yz, Umag_yz, mask_yz, cap_yz_h, cap_yz_v);
samples = [Umag_xz(isfinite(Umag_xz)); Umag_yz(isfinite(Umag_yz))];
planes.samples = samples(:);

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


%% =====================================================================
%  LOCAL FUNCTION: evaluate_plane_velocity
%  =====================================================================
function [Ux, Uy, Uz, Umag, inside_mask, dist_field] = evaluate_plane_velocity( ...
    grid_pts_dl, grid_size, X_cm, R_body, body_vis, ...
    body_data, seg_data, um_s_factor, chunk_size, progress_tag, bg_fun)
%EVALUATE_PLANE_VELOCITY  Stokes velocity on a planar evaluation grid
%
%  Performs body-interior masking and evaluates the RSM (body) plus
%  RSS+Rotlet (pili) velocity kernels at all exterior grid points.  When a
%  background-flow evaluator bg_fun is supplied, the ambient velocity is
%  added so the returned field is the TOTAL fluid velocity.
%
%  Inputs:
%    grid_pts_dl   [N x 3]  Grid points in dimensionless units
%    grid_size     [1 x 2]  [Nrows, Ncols] of the grid
%    X_cm          [3 x 1]  Body CM (dimensionless)
%    R_body        [3 x 3]  Body rotation matrix
%    body_vis      struct   Body visualisation descriptor (see
%                             build_body_visualization): .type,
%                             .b_dl, .c_dl  (parametric) OR
%                             .alpha_shape_dl, .verts_local_dl, .faces (custom)
%    body_data     struct   .X, .fx, .fy, .fz, .e2
%    seg_data      struct   Pili segment data (.active, .y0, .y1, ...)
%    um_s_factor   scalar   Velocity unit conversion factor
%    chunk_size    scalar   Memory-management chunk size
%    progress_tag  char     (optional) Progress bar label, e.g. '[FLOW-XZ] '
%    bg_fun        handle   (optional) X [N x 3] -> [N x 3] dimensionless
%                             background velocity; [] when no ambient flow
%
%  Outputs:
%    Ux, Uy, Uz  [Nrows x Ncols]  Velocity components [um/s] (NaN inside)
%    Umag         [Nrows x Ncols]  Velocity magnitude [um/s]
%    inside_mask  [Nrows x Ncols]  Logical mask (true = interior)
%    dist_field   [N x 1]          Signed-distance-like field used for
%                                   contour extraction.  Level-set convention
%                                   depends on body_vis.type (see
%                                   extract_body_outline).

if nargin < 10 || isempty(progress_tag), progress_tag = ''; end
if nargin < 11, bg_fun = []; end
show_progress = ~isempty(progress_tag);

Nrows = grid_size(1);
Ncols = grid_size(2);

% ---- Body interior masking (mode-aware) ----
[inside_mask_vec, dist_field] = body_inside_mask(grid_pts_dl, X_cm, R_body, body_vis);
inside_mask = reshape(inside_mask_vec, [Nrows, Ncols]);
outside_idx = find(~inside_mask(:));
N_eval      = length(outside_idx);
X_eval      = grid_pts_dl(outside_idx, :);

% ---- Initialise output ----
ux_ev = zeros(N_eval, 1);
uy_ev = zeros(N_eval, 1);
uz_ev = zeros(N_eval, 1);

% ---- Progress bar setup ----
if show_progress
    n_chunks_body = ceil(N_eval / chunk_size);
    n_chunks_pili = 0;
    if seg_data.active
        n_chunks_pili = ceil(N_eval / chunk_size);
    end
    n_chunks_total = n_chunks_body + n_chunks_pili;
    flow_tick      = 0;
    t_flow_bar     = tic;
end

% ---- RSM body contribution ----
X_body = body_data.X;
bfx = body_data.fx;  bfy = body_data.fy;  bfz = body_data.fz;
e2_b = body_data.e2;

for ci = 1:chunk_size:N_eval
    ce  = min(ci + chunk_size - 1, N_eval);
    idx = ci:ce;

    rx = bsxfun(@minus, X_eval(idx,1), X_body(:,1)');
    ry = bsxfun(@minus, X_eval(idx,2), X_body(:,2)');
    rz = bsxfun(@minus, X_eval(idx,3), X_body(:,3)');
    r2 = rx.^2 + ry.^2 + rz.^2;

    denom = (r2 + e2_b).^1.5;
    H1    = (r2 + 2*e2_b) ./ denom;
    H2    = 1 ./ denom;
    rdotf = rx .* bfx' + ry .* bfy' + rz .* bfz';

    ux_ev(idx) = ux_ev(idx) + H1*bfx + sum(H2 .* rx .* rdotf, 2);
    uy_ev(idx) = uy_ev(idx) + H1*bfy + sum(H2 .* ry .* rdotf, 2);
    uz_ev(idx) = uz_ev(idx) + H1*bfz + sum(H2 .* rz .* rdotf, 2);

    if show_progress
        flow_tick = flow_tick + 1;
        render_progress_bar(progress_tag, flow_tick, n_chunks_total, t_flow_bar);
    end
end

% ---- RSS pili + Rotlet contribution ----
if seg_data.active
    e2_p = seg_data.e2;
    L_r  = seg_data.L';   L2_r = seg_data.L2';
    v1_r = seg_data.v(:,1)';
    v2_r = seg_data.v(:,2)';
    v3_r = seg_data.v(:,3)';

    f0x = seg_data.f0(:,1)';  f0y = seg_data.f0(:,2)';  f0z = seg_data.f0(:,3)';
    f1x = seg_data.f1(:,1)';  f1y = seg_data.f1(:,2)';  f1z = seg_data.f1(:,3)';

    for ci = 1:chunk_size:N_eval
        ce  = min(ci + chunk_size - 1, N_eval);
        idx = ci:ce;

        x01 = bsxfun(@minus, X_eval(idx,1), seg_data.y0(:,1)');
        x02 = bsxfun(@minus, X_eval(idx,2), seg_data.y0(:,2)');
        x03 = bsxfun(@minus, X_eval(idx,3), seg_data.y0(:,3)');
        r0_sq = x01.^2 + x02.^2 + x03.^2;

        x11 = x01 + v1_r;  x12 = x02 + v2_r;  x13 = x03 + v3_r;
        r1_sq = x11.^2 + x12.^2 + x13.^2;

        R0 = sqrt(r0_sq + e2_p);
        R1 = sqrt(r1_sq + e2_p);
        d0 = x01.*v1_r + x02.*v2_r + x03.*v3_r;

        H0_v = max(L_r.*R0 + d0,        1e-300);
        H1_v = max(L_r.*R1 + d0 + L2_r,  1e-300);

        inv_L  = 1 ./ L_r;
        inv_L2 = 1 ./ L2_r;
        d0_L2  = d0 .* inv_L2;

        T0m1 = (log(H1_v) - log(H0_v)) .* inv_L;
        T0m3 = 1./(R0.*H0_v) - 1./(R1.*H1_v);
        T1m1 = (R1 - R0) .* inv_L2 - d0_L2 .* T0m1;
        T1m3 = (1./R0 - 1./R1) .* inv_L2 - d0_L2 .* T0m3;
        T2m3 = -1./(R1.*L2_r) + T0m1.*inv_L2 - d0_L2 .* T1m3;
        T3m3 = -1./(R1.*L2_r) + 2*T1m1.*inv_L2 - d0_L2 .* T2m3;

        ci2 = L_r .* (T1m1 + e2_p * T1m3);
        c12 = L_r .* T1m3;
        c22 = L_r .* T2m3;
        c32 = L_r .* T3m3;

        ci1 = L_r .* (T0m1 + e2_p * T0m3) - ci2;
        c11 = L_r .* T0m3 - c12;
        c21 = L_r .* T1m3 - c22;
        c31 = L_r .* T2m3 - c32;

        % M1 * f0 (start-node kernel)
        ux_M1 = sum( ...
            (ci1+c11.*x01.*x01+c21.*2.*x01.*v1_r+c31.*v1_r.*v1_r).*f0x + ...
            (c11.*x01.*x02+c21.*(x01.*v2_r+x02.*v1_r)+c31.*v1_r.*v2_r).*f0y + ...
            (c11.*x01.*x03+c21.*(x01.*v3_r+x03.*v1_r)+c31.*v1_r.*v3_r).*f0z, 2);
        uy_M1 = sum( ...
            (c11.*x02.*x01+c21.*(x02.*v1_r+x01.*v2_r)+c31.*v2_r.*v1_r).*f0x + ...
            (ci1+c11.*x02.*x02+c21.*2.*x02.*v2_r+c31.*v2_r.*v2_r).*f0y + ...
            (c11.*x02.*x03+c21.*(x02.*v3_r+x03.*v2_r)+c31.*v2_r.*v3_r).*f0z, 2);
        uz_M1 = sum( ...
            (c11.*x03.*x01+c21.*(x03.*v1_r+x01.*v3_r)+c31.*v3_r.*v1_r).*f0x + ...
            (c11.*x03.*x02+c21.*(x03.*v2_r+x02.*v3_r)+c31.*v3_r.*v2_r).*f0y + ...
            (ci1+c11.*x03.*x03+c21.*2.*x03.*v3_r+c31.*v3_r.*v3_r).*f0z, 2);

        % M2 * f1 (end-node kernel)
        ux_M2 = sum( ...
            (ci2+c12.*x01.*x01+c22.*2.*x01.*v1_r+c32.*v1_r.*v1_r).*f1x + ...
            (c12.*x01.*x02+c22.*(x01.*v2_r+x02.*v1_r)+c32.*v1_r.*v2_r).*f1y + ...
            (c12.*x01.*x03+c22.*(x01.*v3_r+x03.*v1_r)+c32.*v1_r.*v3_r).*f1z, 2);
        uy_M2 = sum( ...
            (c12.*x02.*x01+c22.*(x02.*v1_r+x01.*v2_r)+c32.*v2_r.*v1_r).*f1x + ...
            (ci2+c12.*x02.*x02+c22.*2.*x02.*v2_r+c32.*v2_r.*v2_r).*f1y + ...
            (c12.*x02.*x03+c22.*(x02.*v3_r+x03.*v2_r)+c32.*v2_r.*v3_r).*f1z, 2);
        uz_M2 = sum( ...
            (c12.*x03.*x01+c22.*(x03.*v1_r+x01.*v3_r)+c32.*v3_r.*v1_r).*f1x + ...
            (c12.*x03.*x02+c22.*(x03.*v2_r+x02.*v3_r)+c32.*v3_r.*v2_r).*f1y + ...
            (ci2+c12.*x03.*x03+c22.*2.*x03.*v3_r+c32.*v3_r.*v3_r).*f1z, 2);

        ux_ev(idx) = ux_ev(idx) + ux_M1 + ux_M2;
        uy_ev(idx) = uy_ev(idx) + uy_M1 + uy_M2;
        uz_ev(idx) = uz_ev(idx) + uz_M1 + uz_M2;

        % Rotlet contribution.  The stored formulas are Cortez's full-curl
        % rotlet; physical torque-induced velocity uses half of this kernel.
        if seg_data.has_rotlet
            R0_3 = R0.^3;  R1_3 = R1.^3;

            d_sq     = max(L2_r.*(r0_sq + e2_p) - d0.^2, e2_p.*L2_r);
            inv_d_sq = 1 ./ d_sq;
            d1_v     = d0 + L2_r;
            inv3L2   = inv_L2 / 3;

            T0m5 = inv_d_sq/3 .* (d1_v./R1_3 - d0./R0_3) ...
                 + 2*L2_r.*inv_d_sq/3 .* T0m3;
            T1m5 = (1./R0_3 - 1./R1_3) .* inv3L2 - d0_L2 .* T0m5;

            ca  = 2*T0m3 + 3*e2_p*T0m5;
            cb  = 2*T1m3 + 3*e2_p*T1m5;
            Rxy =  L_r .* (ca.*x03 + cb.*v3_r);
            Rxz = -L_r .* (ca.*x02 + cb.*v2_r);
            Ryz =  L_r .* (ca.*x01 + cb.*v1_r);

            tx = seg_data.tau_x';
            ty = seg_data.tau_y';
            tz = seg_data.tau_z';

            ux_ev(idx) = ux_ev(idx) + 0.5 * sum( Rxy.*ty + Rxz.*tz, 2);
            uy_ev(idx) = uy_ev(idx) + 0.5 * sum(-Rxy.*tx + Ryz.*tz, 2);
            uz_ev(idx) = uz_ev(idx) + 0.5 * sum(-Rxz.*tx - Ryz.*ty, 2);
        end

        if show_progress
            flow_tick = flow_tick + 1;
            render_progress_bar(progress_tag, flow_tick, n_chunks_total, t_flow_bar);
        end
    end
end

% ---- Add prescribed background flow (total = ambient + disturbance) ----
if ~isempty(bg_fun) && N_eval > 0
    u_bg_ev = bg_fun(X_eval);              % [N_eval x 3], dimensionless
    ux_ev = ux_ev + u_bg_ev(:, 1);
    uy_ev = uy_ev + u_bg_ev(:, 2);
    uz_ev = uz_ev + u_bg_ev(:, 3);
end

% ---- Map to grid and convert units ----
Ux = NaN(Nrows, Ncols);
Uy = NaN(Nrows, Ncols);
Uz = NaN(Nrows, Ncols);
Ux(outside_idx) = ux_ev * um_s_factor;
Uy(outside_idx) = uy_ev * um_s_factor;
Uz(outside_idx) = uz_ev * um_s_factor;
Umag = sqrt(Ux.^2 + Uy.^2 + Uz.^2);

end


%% =====================================================================
%  LOCAL FUNCTION: build_body_visualization
%  =====================================================================
function body_vis = build_body_visualization(params, L_0)
%BUILD_BODY_VISUALIZATION  Mode-aware body descriptor for flow-field plot.
%
%  Returns a struct with the data that evaluate_plane_velocity,
%  build_body_mesh, and extract_body_outline need to render the body and
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

            % Build alphaShape in dimensionless body-local frame (used for
            % inShape() interior masking).  If construction fails, we fall
            % back to a bounding-sphere approximation at evaluation time.
            body_vis.alpha_shape_dl = [];
            try
                shp_local = alphaShape(body_vis.verts_local_dl);
                % Prefer the critical alpha that yields a single watertight
                % region — matches the construction in custom_body_building.
                alpha_crit = criticalAlpha(shp_local, 'all-points');
                if isfinite(alpha_crit) && alpha_crit > 0
                    shp_local.Alpha = alpha_crit;
                end
                body_vis.alpha_shape_dl = shp_local;
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
%
%  The SDF sign convention is: values <= 0 indicate interior, values > 0
%  indicate exterior.  For parametric capsule, dist_field is the classic
%  signed distance (radius - b_dl).  For custom mode, dist_field is the
%  double() of the inside mask minus 0.5, so contourc at level 0 gives
%  the body silhouette on the plane.

    % Transform grid to body-local dimensionless frame
    X_local = (R_body.' * (grid_pts_dl - X_cm(:).').').';

    switch body_vis.type
        case 'capsule'
            x_cl = max(-body_vis.c_dl, min(body_vis.c_dl, X_local(:,1)));
            d_cap = sqrt((X_local(:,1) - x_cl).^2 + ...
                          X_local(:,2).^2 + X_local(:,3).^2);
            inside_mask_vec = d_cap <= body_vis.b_dl;
            dist_field      = d_cap - body_vis.b_dl;     % SDF: < 0 interior

        case 'custom'
            if ~isempty(body_vis.alpha_shape_dl)
                try
                    inside_mask_vec = inShape(body_vis.alpha_shape_dl, ...
                        X_local(:,1), X_local(:,2), X_local(:,3));
                catch
                    inside_mask_vec = bound_sphere_inside(X_local, body_vis);
                end
            else
                inside_mask_vec = bound_sphere_inside(X_local, body_vis);
            end
            % Binary SDF-like field centred on 0.5 -> 0 at boundary
            dist_field = double(inside_mask_vec) - 0.5;   % < 0 exterior
            dist_field = -dist_field;                     % < 0 interior (match capsule convention)

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
%  LOCAL FUNCTION: extract_body_outline
%  =====================================================================
function [cap_h, cap_v] = extract_body_outline(dist_field, grid_size, h_vec, v_vec, body_vis)
%EXTRACT_BODY_OUTLINE  Largest body-silhouette contour on a 2D plane.
%
%  For parametric capsule, dist_field is the classic SDF (< 0 interior),
%  so we use level 0.  For custom mode, dist_field is binary-derived;
%  using level 0 still produces a usable outline because dist_field was
%  shifted to (0.5 - mask).

    dist_grid = reshape(dist_field, grid_size);
    C = contourc(h_vec, v_vec, dist_grid, [0 0]);

    max_pts = 0;
    cap_h   = [];
    cap_v   = [];
    kk = 1;
    while kk < size(C, 2)
        n_pts = C(2, kk);
        if n_pts > max_pts
            max_pts = n_pts;
            cap_h = C(1, kk+1 : kk+n_pts);
            cap_v = C(2, kk+1 : kk+n_pts);
        end
        kk = kk + n_pts + 1;
    end

    % Suppress unused-parameter warning (body_vis reserved for future
    % custom-mode outline refinements, e.g. projected convex hull).
    if isempty(body_vis.type)
        error('flowfield_grid_velocity:BadBodyVis', 'body_vis.type is empty.');
    end
end


function v = get_flow_opt(s, fname, default)
%GET_FLOW_OPT  Read field from options struct with fallback default
    if isfield(s, fname), v = s.(fname); else, v = default; end
end
