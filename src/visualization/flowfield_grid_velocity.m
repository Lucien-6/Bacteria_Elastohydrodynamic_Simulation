function planes = flowfield_grid_velocity(hydro, params, flow_opts, logf, silent)
%FLOWFIELD_GRID_VELOCITY  Evaluate XZ and YZ Stokes velocity planes
%
%  Description:
%    Builds XZ and YZ grids that cover the configured window and the whole
%    cell (body and flagella), and evaluates the regularised Stokeslet,
%    Stokeslet-segment, and rotlet kernels outside the body.  The body
%    interior is left blank (NaN).  The ambient flow is added when
%    hydro.bg_fun is present.
%
%  Syntax:
%    planes = flowfield_grid_velocity(hydro, params, flow_opts, logf, silent)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28
%  Updated: 2026-10-30

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
[lo_um, hi_um] = grid_window_um(hydro, body_vis, X_cm_um, half_range, ...
    grid_spacing_um, L_to_um);

logf(['[FLOW] Grid X [%.2f, %.2f], Y [%.2f, %.2f], Z [%.2f, %.2f] um, ', ...
    'spacing = %.2f um, quiver skip = %d\n'], ...
    lo_um(1), hi_um(1), lo_um(2), hi_um(2), lo_um(3), hi_um(3), ...
    grid_spacing_um, quiver_skip);

logf('[FLOW] Computing XZ-plane flow field...\n');
x_vec_xz = axis_samples(lo_um(1), hi_um(1), grid_spacing_um);
z_vec_xz = axis_samples(lo_um(3), hi_um(3), grid_spacing_um);
[XG_xz, ZG_xz] = meshgrid(x_vec_xz, z_vec_xz);
YG_xz = X_cm_um(2) * ones(size(XG_xz));
xz_tag = '[FLOW-XZ] ';
if silent, xz_tag = ''; end
[Ux_xz, Uy_xz, Uz_xz, Umag_xz, mask_xz, ~] = evaluate_plane_velocity( ...
    [XG_xz(:)*um_to_L, YG_xz(:)*um_to_L, ZG_xz(:)*um_to_L], ...
    [numel(z_vec_xz), numel(x_vec_xz)], X_cm, R_body, body_vis, ...
    hydro.body_data, hydro.seg_data, um_s_factor, chunk_size, xz_tag, hydro.bg_fun);
logf('[FLOW] XZ velocity range: [%.4e, %.4e] um/s\n', ...
    min(Umag_xz(:), [], 'omitnan'), max(Umag_xz(:), [], 'omitnan'));

logf('[FLOW] Computing YZ-plane flow field...\n');
y_vec_yz = axis_samples(lo_um(2), hi_um(2), grid_spacing_um);
z_vec_yz = axis_samples(lo_um(3), hi_um(3), grid_spacing_um);
[YG_yz, ZG_yz] = meshgrid(y_vec_yz, z_vec_yz);
XG_yz = X_cm_um(1) * ones(size(YG_yz));
yz_tag = '[FLOW-YZ] ';
if silent, yz_tag = ''; end
[Ux_yz, Uy_yz, Uz_yz, Umag_yz, mask_yz, ~] = evaluate_plane_velocity( ...
    [XG_yz(:)*um_to_L, YG_yz(:)*um_to_L, ZG_yz(:)*um_to_L], ...
    [numel(z_vec_yz), numel(y_vec_yz)], X_cm, R_body, body_vis, ...
    hydro.body_data, hydro.seg_data, um_s_factor, chunk_size, yz_tag, hydro.bg_fun);
logf('[FLOW] YZ velocity range: [%.4e, %.4e] um/s\n', ...
    min(Umag_yz(:), [], 'omitnan'), max(Umag_yz(:), [], 'omitnan'));

[cap_xz_h, cap_xz_v] = body_section_outline(body_vis, R_body, X_cm_um, 'xz');
[cap_yz_h, cap_yz_v] = body_section_outline(body_vis, R_body, X_cm_um, 'yz');

planes.body_vis = body_vis;
planes.axis_lim_um = [x_vec_xz(1), x_vec_xz(end); ...
                      y_vec_yz(1), y_vec_yz(end); ...
                      z_vec_xz(1), z_vec_xz(end)];
planes.half_range = max(abs(planes.axis_lim_um - X_cm_um), [], 'all');
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
function [h, v] = body_section_outline(body_vis, R_body, X_cm_um, plane_name)
%BODY_SECTION_OUTLINE  Smooth body silhouette on the XZ or YZ plane
%
%  The curve is the plane cut of the body surface, resampled with a
%  periodic spline.  A contour of the velocity-grid mask is stair-stepped
%  at the flow spacing and is not used.

    [V, F] = body_triangles_um(body_vis, R_body, X_cm_um);
    if strcmp(plane_name, 'xz')
        fixed_axis = 2;
        fixed_value = X_cm_um(2);
        h_axis = 1;
        v_axis = 3;
    else
        fixed_axis = 1;
        fixed_value = X_cm_um(1);
        h_axis = 2;
        v_axis = 3;
    end
    [h, v] = mesh_plane_loops(V, F, fixed_axis, fixed_value, h_axis, v_axis);
    if isempty(body_vis.type)
        error('flowfield_grid_velocity:BadBodyVis', 'body_vis.type is empty.');
    end
end


function [V, F] = body_triangles_um(body_vis, R_body, X_cm_um)
%BODY_TRIANGLES_UM  Closed body surface in lab-frame micrometres
    if strcmp(body_vis.type, 'custom') && ~isempty(body_vis.faces)
        V = (R_body * (body_vis.verts_m * 1e6).' + X_cm_um(:)).';
        F = body_vis.faces;
        return
    end
    a_um = body_vis.a_um;
    b_um = body_vis.b_um;
    n_circ = 72;
    n_cap = 18;
    n_cyl = 16;
    theta = linspace(0, 2 * pi, n_circ + 1);
    theta = theta(1:end - 1);
    phi_L = linspace(pi, pi / 2, n_cap);
    [TH_L, PH_L] = meshgrid(theta, phi_L);
    XL = b_um * cos(PH_L) - (a_um - b_um);
    YL = b_um * sin(PH_L) .* cos(TH_L);
    ZL = b_um * sin(PH_L) .* sin(TH_L);
    x_cyl = linspace(-(a_um - b_um), (a_um - b_um), n_cyl);
    [TH_C, XC] = meshgrid(theta, x_cyl);
    YC = b_um * cos(TH_C);
    ZC = b_um * sin(TH_C);
    phi_R = linspace(pi / 2, 0, n_cap);
    [TH_R, PH_R] = meshgrid(theta, phi_R);
    XR = b_um * cos(PH_R) + (a_um - b_um);
    YR = b_um * sin(PH_R) .* cos(TH_R);
    ZR = b_um * sin(PH_R) .* sin(TH_R);
    % Drop the duplicated equator shared by a cap and the cylinder.
    Xg = [XL(1:end - 1, :); XC; XR(2:end, :)];
    Yg = [YL(1:end - 1, :); YC; YR(2:end, :)];
    Zg = [ZL(1:end - 1, :); ZC; ZR(2:end, :)];
    Vloc = [Xg(:), Yg(:), Zg(:)];
    V = (R_body * Vloc.' + X_cm_um(:)).';
    nr = size(Xg, 1);
    nc = size(Xg, 2);
    F = zeros(2 * nr * nc, 3);
    nf = 0;
    for i = 1:nr - 1
        for j = 1:nc
            j2 = mod(j, nc) + 1;
            i00 = (j - 1) * nr + i;
            i10 = (j2 - 1) * nr + i;
            i01 = i00 + 1;
            i11 = i10 + 1;
            nf = nf + 1;
            F(nf, :) = [i00, i10, i11];
            nf = nf + 1;
            F(nf, :) = [i00, i11, i01];
        end
    end
    F = F(1:nf, :);
end


function [h, v] = mesh_plane_loops(V, F, fixed_axis, fixed_value, h_axis, v_axis)
%MESH_PLANE_LOOPS  Smooth closed cuts of a triangle mesh by one plane
    nF = size(F, 1);
    P = zeros(nF, 2);
    Q = zeros(nF, 2);
    nseg = 0;
    edges = [1, 2; 2, 3; 3, 1];
    for i = 1:nF
        tri = V(F(i, :), :);
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
                d1 = norm(P(j, :) - tail);
                d2 = norm(Q(j, :) - tail);
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


function [lo_um, hi_um] = grid_window_um(hydro, body_vis, X_cm_um, half_range, spacing_um, L_to_um)
%GRID_WINDOW_UM  Axis-aligned window covering the configured range and the cell
%
%  The configured window is the cube of side grid_range centred on the body.
%  It is expanded, independently on each side, so the body surface and every
%  flagellum node lie inside with a margin of at least 0.5 um.

    n_body = size(hydro.body_data.X, 1);
    n_pili = 0;
    n_mesh = 0;
    has_pili = isfield(hydro, 'geom') && isfield(hydro.geom, 'X_pili');
    if has_pili
        for k = 1:numel(hydro.geom.X_pili)
            n_pili = n_pili + size(hydro.geom.X_pili{k}, 1);
        end
    end
    use_mesh = strcmp(body_vis.type, 'custom') && ~isempty(body_vis.verts_local_dl);
    if use_mesh
        n_mesh = size(body_vis.verts_local_dl, 1);
    end
    pts = zeros(n_body + n_pili + n_mesh, 3);
    pts(1:n_body, :) = hydro.body_data.X;
    row = n_body + 1;
    if has_pili
        for k = 1:numel(hydro.geom.X_pili)
            nk = size(hydro.geom.X_pili{k}, 1);
            pts(row:row + nk - 1, :) = hydro.geom.X_pili{k};
            row = row + nk;
        end
    end
    if use_mesh
        pts(row:end, :) = (hydro.R_body * body_vis.verts_local_dl.' + hydro.X_cm(:)).';
    end
    pts_um = pts * L_to_um;
    bb_min = min(pts_um, [], 1);
    bb_max = max(pts_um, [], 1);
    pad = max(0.5, 4 * spacing_um);
    lo_um = min(X_cm_um(:).' - half_range, bb_min - pad);
    hi_um = max(X_cm_um(:).' + half_range, bb_max + pad);
end


function v = axis_samples(lo, hi, spacing)
%AXIS_SAMPLES  Uniform samples from lo through a point at or beyond hi
    if ~(isfinite(lo) && isfinite(hi) && isfinite(spacing) && spacing > 0)
        error('flowfield_grid_velocity:BadGrid', ...
            'Grid limits and spacing must be finite, and spacing must be positive.');
    end
    if hi < lo
        tmp = lo;
        lo = hi;
        hi = tmp;
    end
    n = max(1, ceil((hi - lo) / spacing));
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
