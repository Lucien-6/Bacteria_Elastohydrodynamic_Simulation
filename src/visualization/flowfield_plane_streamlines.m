function curves = flowfield_plane_streamlines(h_vec, v_vec, Uh, Uv, mask, ...
    origin, spacing_um, length_um, plane_name, skip_check)
%FLOWFIELD_PLANE_STREAMLINES  Evenly spaced integral curves of one plane
%
%  Description:
%    Integrates the unit in-plane velocity (Uh, Uv).  Speed stays in the
%    colour map; these curves carry direction only.  They are integral
%    curves of the planar components, not paths of 3D particles.
%
%    Seeds lie on a lattice of spacing_um, ordered from the origin outward.
%    Each seed is traced forward and backward with steps of 0.2 grid cells.
%    The step is shortened only when 0.2 cells would exceed a quarter of the
%    stop distance, and the last step is shortened to land on the arc cap.
%    A running curve stops at the body, at a stagnation point, at the
%    domain edge, at length_um, when it meets another kept curve, or when
%    it loops back onto itself.  The meeting gap is a few steps, not half
%    a seed spacing.  A new seed inside half a seed spacing of a curve
%    already kept is not started.  Fragments shorter than half a seed
%    spacing are dropped.  A solid arrowhead sits at the arc-length midpoint
%    and points downstream.
%
%  Syntax:
%    curves = flowfield_plane_streamlines(h_vec, v_vec, Uh, Uv, mask, ...
%        origin, spacing_um, length_um)
%    curves = flowfield_plane_streamlines(..., plane_name)
%    curves = flowfield_plane_streamlines(..., plane_name, skip_check)
%
%  Output fields:
%    x, y         - NaN-separated polyline coordinates [um]
%    vertices     - Arrowhead vertices [um], N-by-2
%    faces        - Arrowhead faces, M-by-3, tip then the two base corners
%    polylines    - Cell of kept curves, each K-by-2
%    n_lines      - Number of kept curves
%    spacing_um, length_um, step_um, stop_um
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-13
%  Updated: 2026-11-24 — zero-copy backward trace, endpoint sample reuse,
%           inlined RK2 sampling, dilated fine-grid O(1) collision filter,
%           reusable loop-stamp buffer, and batch bin updates

if nargin < 9 || isempty(plane_name)
    plane_name = 'plane';
end
if nargin < 10 || isempty(skip_check)
    skip_check = false;
end
plane_name = char(plane_name);
skip_check = logical(skip_check);

[h_vec, v_vec, Uh, Uv, mask, origin, spacing_um, length_um, dh, dv] = ...
    validate_inputs(h_vec, v_vec, Uh, Uv, mask, origin, ...
    spacing_um, length_um, plane_name);

curves = empty_curves(spacing_um, length_um, dh);
[grid, seeds] = prepare_grid(h_vec, v_vec, Uh, Uv, mask, origin, ...
    spacing_um, length_um, dh, dv, plane_name);
if isempty(seeds) || ~any(isfinite(grid.U(:)))
    return
end
curves.step_um = grid.ds;
curves.stop_um = grid.d_stop;

[kept, n_keep] = place_curves(grid, seeds);
if n_keep == 0
    return
end
curves = pack_curves(kept, n_keep, spacing_um, length_um, grid.ds, grid.d_stop);
if ~skip_check
    check_result(curves, h_vec, v_vec, mask, plane_name);
end
end


function [h_vec, v_vec, Uh, Uv, mask, origin, spacing_um, length_um, dh, dv] = ...
    validate_inputs(h_vec, v_vec, Uh, Uv, mask, origin, ...
    spacing_um, length_um, plane_name)
%VALIDATE_INPUTS  Grid, mask, and streamline scales for one plane
if ~isnumeric(h_vec) || ~isvector(h_vec) || numel(h_vec) < 2 || ~isreal(h_vec)
    error('flowfield_plane_streamlines:BadGrid', ...
        'Plane %s horizontal coordinates must be a real vector of length >= 2.', ...
        plane_name);
end
if ~isnumeric(v_vec) || ~isvector(v_vec) || numel(v_vec) < 2 || ~isreal(v_vec)
    error('flowfield_plane_streamlines:BadGrid', ...
        'Plane %s vertical coordinates must be a real vector of length >= 2.', ...
        plane_name);
end
h_vec = double(h_vec(:));
v_vec = double(v_vec(:));
dh_all = diff(h_vec);
dv_all = diff(v_vec);
if any(~isfinite(h_vec)) || any(~isfinite(v_vec)) || any(dh_all <= 0) || any(dv_all <= 0)
    error('flowfield_plane_streamlines:BadGrid', ...
        'Plane %s coordinates must be finite and strictly increasing.', plane_name);
end
dh = median(dh_all);
dv = median(dv_all);
if max(abs(dh_all - dh)) > 1e-6 * dh || max(abs(dv_all - dv)) > 1e-6 * dv
    error('flowfield_plane_streamlines:BadGrid', ...
        'Plane %s coordinates must be uniformly spaced.', plane_name);
end

nH = numel(h_vec);
nV = numel(v_vec);
if ~isnumeric(Uh) || ~isreal(Uh) || ~isequal(size(Uh), [nV, nH]) ...
        || ~isnumeric(Uv) || ~isreal(Uv) || ~isequal(size(Uv), [nV, nH])
    error('flowfield_plane_streamlines:BadVelocity', ...
        'Plane %s velocity must be real arrays of size [numel(v), numel(h)].', ...
        plane_name);
end
Uh = double(Uh);
Uv = double(Uv);
if any(isinf(Uh(:))) || any(isinf(Uv(:)))
    error('flowfield_plane_streamlines:BadVelocity', ...
        'Plane %s in-plane velocity contains Inf.', plane_name);
end
if ~(islogical(mask) && isequal(size(mask), [nV, nH]))
    if isnumeric(mask) && isreal(mask) && isequal(size(mask), [nV, nH]) ...
            && all(mask(:) == 0 | mask(:) == 1)
        mask = logical(mask);
    else
        error('flowfield_plane_streamlines:BadGrid', ...
            'Plane %s mask must be a logical array, or 0/1, of size [numel(v), numel(h)].', ...
            plane_name);
    end
end

origin = double(origin(:));
if numel(origin) ~= 2 || any(~isfinite(origin))
    error('flowfield_plane_streamlines:BadGrid', ...
        'Plane %s origin must be a finite 2-vector.', plane_name);
end
spacing_um = require_positive(spacing_um, 'spacing', plane_name, ...
    'flowfield_plane_streamlines:BadSpacing');
length_um = require_positive(length_um, 'length', plane_name, ...
    'flowfield_plane_streamlines:BadLength');

span = min(h_vec(end) - h_vec(1), v_vec(end) - v_vec(1));
if ~(spacing_um < span)
    error('flowfield_plane_streamlines:BadSpacing', ...
        'Plane %s spacing %.4g um must be smaller than the shorter side %.4g um.', ...
        plane_name, spacing_um, span);
end
if length_um < 0.5 * spacing_um
    error('flowfield_plane_streamlines:BadLength', ...
        ['Plane %s length %.4g um must be at least half the spacing ', ...
         '(%.4g um), otherwise every curve is discarded.'], ...
        plane_name, length_um, 0.5 * spacing_um);
end
end


function value = require_positive(value, label, plane_name, err_id)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
    error(err_id, ...
        'Plane %s streamline %s must be a positive finite scalar.', ...
        plane_name, label);
end
value = double(value);
end


function curves = empty_curves(spacing_um, length_um, dh)
curves.x = NaN;
curves.y = NaN;
curves.vertices = [NaN, NaN; NaN, NaN; NaN, NaN];
curves.faces = [1, 1, 1];
curves.polylines = {};
curves.n_lines = 0;
curves.spacing_um = spacing_um;
curves.length_um = length_um;
curves.step_um = 0.2 * dh;
curves.stop_um = 0.5 * spacing_um;
end


function [grid, seeds] = prepare_grid(h_vec, v_vec, Uh, Uv, mask, origin, ...
    spacing_um, length_um, dh, dv, plane_name)
%PREPARE_GRID  Unit field, seed lattice, and arc-length step
spd = hypot(Uh, Uv);
scale = max(spd(:), [], 'omitnan');
bad = mask | ~isfinite(spd);
if isfinite(scale) && scale > 0
    bad = bad | spd <= 1e-8 * scale;
end
Uhat = NaN(size(Uh));
Vhat = NaN(size(Uv));
good = ~bad & spd > 0;
Uhat(good) = Uh(good) ./ spd(good);
Vhat(good) = Uv(good) ./ spd(good);

d_stop = 0.5 * spacing_um;
ds = min([0.2 * dh, 0.2 * dv, 0.25 * d_stop]);
n_max = floor(length_um / ds) + 3;
if n_max > 50000
    error('flowfield_plane_streamlines:BadLength', ...
        ['Plane %s would take %d steps. Increase the grid spacing or ', ...
         'reduce the streamline length.'], plane_name, n_max);
end

span_h = h_vec(end) - h_vec(1);
span_v = v_vec(end) - v_vec(1);
n_h = floor(span_h / spacing_um);
n_v = floor(span_v / spacing_um);
if n_h < 1 || n_v < 1
    error('flowfield_plane_streamlines:BadSpacing', ...
        'Plane %s spacing %.4g um leaves no seed on the grid.', ...
        plane_name, spacing_um);
end
if n_h * n_v > 20000
    error('flowfield_plane_streamlines:BadSpacing', ...
        ['Plane %s spacing %.4g um gives %d seeds (limit 20000). ', ...
         'Increase stream_spacing_um.'], plane_name, spacing_um, n_h * n_v);
end
covered_h = n_h * spacing_um;
covered_v = n_v * spacing_um;
hs = h_vec(1) + (span_h - covered_h) / 2 + spacing_um / 2 ...
    + (0:n_h - 1)' * spacing_um;
vs = v_vec(1) + (span_v - covered_v) / 2 + spacing_um / 2 ...
    + (0:n_v - 1)' * spacing_um;
[Hs, Vs] = ndgrid(hs, vs);
seeds = [Hs(:), Vs(:)];
dil = dilate8(mask);
ih = nearest_index(h_vec, seeds(:, 1), dh);
iv = nearest_index(v_vec, seeds(:, 2), dv);
keep = ~dil(sub2ind(size(dil), iv, ih));
seeds = seeds(keep, :);
if ~isempty(seeds)
    d_origin = [seeds(:, 1) - origin(1), seeds(:, 2) - origin(2)];
    keys = [sum(d_origin.^2, 2), atan2(d_origin(:, 2), d_origin(:, 1)), seeds];
    [~, order] = sortrows(keys);
    seeds = seeds(order, :);
end

grid.U = Uhat;
grid.V = Vhat;
grid.h0 = h_vec(1);
grid.v0 = v_vec(1);
grid.dh = dh;
grid.dv = dv;
grid.nH = numel(h_vec);
grid.nV = numel(v_vec);
grid.ds = ds;
grid.n_max = n_max;
grid.length_um = length_um;
grid.d_stop = d_stop;
grid.bin_h0 = h_vec(1);
grid.bin_v0 = v_vec(1);
grid.nCh = max(1, ceil(span_h / d_stop));
grid.nCv = max(1, ceil(span_v / d_stop));
grid.seed_d2 = d_stop ^ 2;
% Loop test uses a few steps, not the seed spacing, so a bend near the
% flagellum is not treated as a closed orbit.
grid.loop_tol = max(6 * ds, 0.05);
grid.nLoopH = max(1, ceil(span_h / grid.loop_tol)) + 1;
grid.nLoopV = max(1, ceil(span_v / grid.loop_tol)) + 1;
grid.loop_recent = max(4, ceil(0.5 / grid.loop_tol));
end


function [kept, n_keep] = place_curves(grid, seeds)
%PLACE_CURVES  Accept seeds whose curves stay half a spacing apart
n_seeds = size(seeds, 1);
kept = cell(n_seeds, 1);
n_keep = 0;

nCh = grid.nCh;
nCv = grid.nCv;
bins_x = cell(nCv * nCh, 1);
bins_y = cell(nCv * nCh, 1);
% Padded 1-cell-dilated fine occupancy bitmap at loop_tol resolution:
% occ_dilated(iv2 + 1, ih2 + 1) == false mathematically guarantees no kept
% sample lies within < loop_tol of (h2, v2), skipping bin checks in O(1).
occ_dilated = false(grid.nLoopV + 2, grid.nLoopH + 2);
% Reusable stamp array for self-loop detection (avoids per-ray allocation)
visit_stamp = zeros(grid.nLoopV, grid.nLoopH, 'int32');
ray_id = int32(0);

seed_pad = ceil(sqrt(grid.seed_d2) / grid.d_stop) + 1;
inv_dstop = 1 / grid.d_stop;
bin_h0 = grid.bin_h0;
bin_v0 = grid.bin_v0;
h0 = grid.h0;
v0 = grid.v0;
inv_loop = 1 / grid.loop_tol;
nLoopH = grid.nLoopH;
nLoopV = grid.nLoopV;
nLV2 = nLoopV + 2;
empty_blk = zeros(0, 2);

for s = 1:n_seeds
    h = seeds(s, 1);
    v = seeds(s, 2);
    ih = floor((h - bin_h0) * inv_dstop + 1e-12) + 1;
    iv = floor((v - bin_v0) * inv_dstop + 1e-12) + 1;
    if ih < 1, ih = 1; elseif ih > nCh, ih = nCh; end
    if iv < 1, iv = 1; elseif iv > nCv, iv = nCv; end
    if has_point_within(bins_x, bins_y, h, v, ih, iv, nCh, nCv, seed_pad, grid.seed_d2)
        continue
    end

    ray_id = ray_id + int32(1);
    [fwd, visit_stamp] = trace_ray(grid, h, v, 1.0, bins_x, bins_y, ...
        occ_dilated, empty_blk, visit_stamp, ray_id);
    blockers = forward_blockers(fwd, grid.d_stop);
    ray_id = ray_id + int32(1);
    [bwd, visit_stamp] = trace_ray(grid, h, v, -1.0, bins_x, bins_y, ...
        occ_dilated, blockers, visit_stamp, ray_id);

    poly = join_rays(bwd, fwd);
    if size(poly, 1) < 2 || polyline_length(poly) < grid.d_stop
        continue
    end
    n_keep = n_keep + 1;
    kept{n_keep} = poly;

    % Batch-insert polyline points into spatial bins
    px = poly(:, 1);
    py = poly(:, 2);
    ih_b = min(max(floor((px - bin_h0) * inv_dstop + 1e-12) + 1, 1), nCh);
    iv_b = min(max(floor((py - bin_v0) * inv_dstop + 1e-12) + 1, 1), nCv);
    lin_b = iv_b + (ih_b - 1) * nCv;
    [u_lin, ~, ic] = unique(lin_b);
    for m = 1:numel(u_lin)
        li = u_lin(m);
        sel = (ic == m);
        bins_x{li} = [bins_x{li}; px(sel)];
        bins_y{li} = [bins_y{li}; py(sel)];
    end

    % Batch-mark 3x3 dilated fine occupancy cells
    ih_l = min(max(floor((px - h0) * inv_loop + 1e-12) + 1, 1), nLoopH);
    iv_l = min(max(floor((py - v0) * inv_loop + 1e-12) + 1, 1), nLoopV);
    base_lin = (iv_l + 1) + ih_l * nLV2;
    occ_dilated(base_lin)            = true;
    occ_dilated(base_lin - 1)        = true;
    occ_dilated(base_lin + 1)        = true;
    occ_dilated(base_lin - nLV2)     = true;
    occ_dilated(base_lin - nLV2 - 1) = true;
    occ_dilated(base_lin - nLV2 + 1) = true;
    occ_dilated(base_lin + nLV2)     = true;
    occ_dilated(base_lin + nLV2 - 1) = true;
    occ_dilated(base_lin + nLV2 + 1) = true;
end
kept = kept(1:n_keep);
end


function pts = forward_blockers(fwd, d_stop)
%FORWARD_BLOCKERS  Forward samples outside the seed neighbourhood
if size(fwd, 1) < 2
    pts = zeros(0, 2);
    return
end
arc = [0; cumsum(hypot(diff(fwd(:, 1)), diff(fwd(:, 2))))];
pts = fwd(arc >= d_stop, :);
end


function poly = join_rays(bwd, fwd)
%JOIN_RAYS  Backward tail, then the seed, then the forward ray
if size(bwd, 1) <= 1
    poly = fwd;
else
    poly = [bwd(end:-1:2, :); fwd];
end
end


function curves = pack_curves(kept, n_keep, spacing_um, length_um, ds, d_stop)
%PACK_CURVES  NaN-separated coordinates and one arrowhead per curve
counts = zeros(n_keep, 1);
for i = 1:n_keep
    counts(i) = size(kept{i}, 1);
end
n_coord = sum(counts) + n_keep;
x = nan(n_coord, 1);
y = nan(n_coord, 1);
vertices = zeros(n_keep * 3, 2);
pos = 1;
for i = 1:n_keep
    poly = kept{i};
    np = counts(i);
    x(pos:pos + np - 1) = poly(:, 1);
    y(pos:pos + np - 1) = poly(:, 2);
    pos = pos + np + 1;
    vertices(3 * i - 2:3 * i, :) = arrowhead(poly, spacing_um);
end
curves.x = x;
curves.y = y;
curves.vertices = vertices;
curves.faces = reshape(1:n_keep * 3, 3, n_keep).';
curves.polylines = kept;
curves.n_lines = n_keep;
curves.spacing_um = spacing_um;
curves.length_um = length_um;
curves.step_um = ds;
curves.stop_um = d_stop;
end


function tri = arrowhead(poly, spacing_um)
%ARROWHEAD  Triangle centred on the arc-length midpoint, pointing downstream.
[mid, d, arc_len] = curve_midpoint(poly);
head_len = 0.35 * spacing_um;
if head_len > 0.8 * arc_len
    head_len = 0.8 * arc_len;
end
half_w = 0.42 * head_len;
normal = [-d(2), d(1)];
tip = mid + d * (0.5 * head_len);
base = mid - d * (0.5 * head_len);
tri = [tip; base + normal * half_w; base - normal * half_w];
end


function [mid, direction, arc_len] = curve_midpoint(poly)
%CURVE_MIDPOINT  Point and downstream tangent at half the arc length
seg = hypot(diff(poly(:, 1)), diff(poly(:, 2)));
arc = [0; cumsum(seg)];
arc_len = arc(end);
half = 0.5 * arc_len;
k = find(arc >= half - 1e-12, 1, 'first');
if isempty(k) || k < 2
    k = size(poly, 1);
end
span = arc(k) - arc(k - 1);
if span < 1e-12
    t = 0;
else
    t = (half - arc(k - 1)) / span;
    t = min(max(t, 0), 1);
end
mid = poly(k - 1, :) + t * (poly(k, :) - poly(k - 1, :));
d = poly(k, :) - poly(k - 1, :);
dn = hypot(d(1), d(2));
if dn < 1e-12
    direction = [1, 0];
else
    direction = d / dn;
end
end


function check_result(curves, h_vec, v_vec, mask, plane_name)
%CHECK_RESULT  Arc-length cap, body exclusion, and arrowhead indexing
dh = h_vec(2) - h_vec(1);
dv = v_vec(2) - v_vec(1);
cap = 2 * curves.length_um + 1e-6;
for i = 1:curves.n_lines
    poly = curves.polylines{i};
    if size(poly, 1) < 2 || any(~isfinite(poly(:)))
        error('flowfield_plane_streamlines:BadCurve', ...
            'Plane %s kept a degenerate streamline.', plane_name);
    end
    arc = polyline_length(poly);
    if arc > cap || arc + 1e-9 < curves.stop_um
        error('flowfield_plane_streamlines:BadCurve', ...
            'Plane %s streamline %d has arc length %.4g um, outside [%.4g, %.4g].', ...
            plane_name, i, arc, curves.stop_um, curves.length_um * 2);
    end
    ih = nearest_index(h_vec, poly(:, 1), dh);
    iv = nearest_index(v_vec, poly(:, 2), dv);
    if any(mask(sub2ind(size(mask), iv, ih)))
        error('flowfield_plane_streamlines:BadCurve', ...
            'Plane %s streamline %d enters the body mask.', plane_name, i);
    end
    tip = curves.vertices(curves.faces(i, 1), :);
    base = 0.5 * (curves.vertices(curves.faces(i, 2), :) ...
        + curves.vertices(curves.faces(i, 3), :));
    centre = 0.5 * (tip + base);
    [mid, tangent] = curve_midpoint(poly);
    if max(abs(centre - mid)) > 1e-6
        error('flowfield_plane_streamlines:BadCurve', ...
            'Plane %s arrowhead %d is not at the arc-length midpoint.', ...
            plane_name, i);
    end
    direction = tip - base;
    if dot(direction, tangent) <= 0
        error('flowfield_plane_streamlines:BadCurve', ...
            'Plane %s arrowhead %d does not point downstream.', plane_name, i);
    end
end
if curves.n_lines > 0
    n_vert = size(curves.vertices, 1);
    faces = curves.faces;
    if size(faces, 1) ~= curves.n_lines || size(faces, 2) ~= 3 ...
            || any(faces(:) < 1) || any(faces(:) > n_vert)
        error('flowfield_plane_streamlines:BadCurve', ...
            'Plane %s arrowhead faces do not index the vertices.', plane_name);
    end
end
end


function [pts, visit_stamp] = trace_ray(grid, h, v, dir_sign, bins_x, bins_y, ...
    occ_dilated, blockers, visit_stamp, ray_id)
%TRACE_RAY  RK2 until the boundary, the arc cap, a loop, or another curve
%  Inlines bilinear interpolation and reuses endpoint velocity samples on
%  the next step for identical RK2 trajectories with 1/3 fewer lookups.

U = grid.U;
V = grid.V;
h0 = grid.h0;
v0 = grid.v0;
inv_dh = 1 / grid.dh;
inv_dv = 1 / grid.dv;
nH = grid.nH;
nV = grid.nV;
nHm1 = nH - 1;
nVm1 = nV - 1;
ds = grid.ds;
n_max = grid.n_max;
length_um = grid.length_um;
length_lim = length_um - 1e-12;
d_stop = grid.d_stop;
bin_h0 = grid.bin_h0;
bin_v0 = grid.bin_v0;
inv_dstop = 1 / d_stop;
nCh = grid.nCh;
nCv = grid.nCv;
loop_tol = grid.loop_tol;
inv_loop = 1 / loop_tol;
loop2 = loop_tol * loop_tol;
nLoopH = grid.nLoopH;
nLoopV = grid.nLoopV;
recent_n = grid.loop_recent;

pts = zeros(n_max, 2);
pts(1, 1) = h;
pts(1, 2) = v;
n = 1;
arc = 0;

% Initial direction sample at the seed (reused by step 1)
fh = (h - h0) * inv_dh;
fv = (v - v0) * inv_dv;
if fh < 0 || fv < 0 || fh > nHm1 || fv > nVm1
    pts = pts(1, :);
    return
end
ih = floor(fh) + 1;
iv = floor(fv) + 1;
if ih >= nH, ih = nHm1; th = 1; else, th = fh - (ih - 1); end
if iv >= nV, iv = nVm1; tv = 1; else, tv = fv - (iv - 1); end
i00 = iv + (ih - 1) * nV;
i10 = i00 + nV;
om_th = 1 - th;  om_tv = 1 - tv;
u_val = om_tv * (om_th * U(i00) + th * U(i10)) + tv * (om_th * U(i00 + 1) + th * U(i10 + 1));
w_val = om_tv * (om_th * V(i00) + th * V(i10)) + tv * (om_th * V(i00 + 1) + th * V(i10 + 1));
sp = hypot(u_val, w_val);
if ~(sp > 1e-8)
    pts = pts(1, :);
    return
end
inv_sp = dir_sign / sp;
u1 = u_val * inv_sp;
w1 = w_val * inv_sp;

cur_ih = floor((h - h0) * inv_loop + 1e-12) + 1;
cur_iv = floor((v - v0) * inv_loop + 1e-12) + 1;
if cur_ih < 1, cur_ih = 1; elseif cur_ih > nLoopH, cur_ih = nLoopH; end
if cur_iv < 1, cur_iv = 1; elseif cur_iv > nLoopV, cur_iv = nLoopV; end

buf_h = zeros(recent_n, 1);
buf_v = zeros(recent_n, 1);
n_left = 0;

has_block = ~isempty(blockers);
if has_block
    blk_x = blockers(:, 1);
    blk_y = blockers(:, 2);
    blk_xmin = min(blk_x) - loop_tol;
    blk_xmax = max(blk_x) + loop_tol;
    blk_ymin = min(blk_y) - loop_tol;
    blk_ymax = max(blk_y) + loop_tol;
end

while arc < length_lim && n < n_max
    step = length_um - arc;
    if ds < step
        step = ds;
    end

    % RK2 midpoint sample
    hm = h + 0.5 * step * u1;
    vm = v + 0.5 * step * w1;
    fh = (hm - h0) * inv_dh;
    fv = (vm - v0) * inv_dv;
    if fh < 0 || fv < 0 || fh > nHm1 || fv > nVm1
        break
    end
    ih = floor(fh) + 1;
    iv = floor(fv) + 1;
    if ih >= nH, ih = nHm1; th = 1; else, th = fh - (ih - 1); end
    if iv >= nV, iv = nVm1; tv = 1; else, tv = fv - (iv - 1); end
    i00 = iv + (ih - 1) * nV;
    i10 = i00 + nV;
    om_th = 1 - th;  om_tv = 1 - tv;
    u_val = om_tv * (om_th * U(i00) + th * U(i10)) + tv * (om_th * U(i00 + 1) + th * U(i10 + 1));
    w_val = om_tv * (om_th * V(i00) + th * V(i10)) + tv * (om_th * V(i00 + 1) + th * V(i10 + 1));
    sp = hypot(u_val, w_val);
    if ~(sp > 1e-8)
        break
    end
    inv_sp = dir_sign / sp;
    u2 = u_val * inv_sp;
    w2 = w_val * inv_sp;

    % Endpoint sample (reused as u1, w1 on the next iteration)
    h2 = h + step * u2;
    v2 = v + step * w2;
    fh = (h2 - h0) * inv_dh;
    fv = (v2 - v0) * inv_dv;
    if fh < 0 || fv < 0 || fh > nHm1 || fv > nVm1
        break
    end
    ih = floor(fh) + 1;
    iv = floor(fv) + 1;
    if ih >= nH, ih = nHm1; th = 1; else, th = fh - (ih - 1); end
    if iv >= nV, iv = nVm1; tv = 1; else, tv = fv - (iv - 1); end
    i00 = iv + (ih - 1) * nV;
    i10 = i00 + nV;
    om_th = 1 - th;  om_tv = 1 - tv;
    u_val = om_tv * (om_th * U(i00) + th * U(i10)) + tv * (om_th * U(i00 + 1) + th * U(i10 + 1));
    w_val = om_tv * (om_th * V(i00) + th * V(i10)) + tv * (om_th * V(i00 + 1) + th * V(i10 + 1));
    sp = hypot(u_val, w_val);
    if ~(sp > 1e-8)
        break
    end
    inv_sp = dir_sign / sp;
    u_next = u_val * inv_sp;
    w_next = w_val * inv_sp;

    % Fine-grid loop cell (shared by O(1) collision filter and self-loop test)
    ih2 = floor((h2 - h0) * inv_loop + 1e-12) + 1;
    iv2 = floor((v2 - v0) * inv_loop + 1e-12) + 1;
    if ih2 < 1, ih2 = 1; elseif ih2 > nLoopH, ih2 = nLoopH; end
    if iv2 < 1, iv2 = 1; elseif iv2 > nLoopV, iv2 = nLoopV; end

    % Fast O(1) skip unless dilated fine occupancy signals a nearby curve
    if occ_dilated(iv2 + 1, ih2 + 1)
        ihb = floor((h2 - bin_h0) * inv_dstop + 1e-12) + 1;
        ivb = floor((v2 - bin_v0) * inv_dstop + 1e-12) + 1;
        if ihb < 1, ihb = 1; elseif ihb > nCh, ihb = nCh; end
        if ivb < 1, ivb = 1; elseif ivb > nCv, ivb = nCv; end
        if has_point_within(bins_x, bins_y, h2, v2, ihb, ivb, nCh, nCv, 1, loop2)
            break
        end
    end

    if has_block && (arc + step) >= d_stop ...
            && h2 >= blk_xmin && h2 <= blk_xmax ...
            && v2 >= blk_ymin && v2 <= blk_ymax ...
            && any((blk_x - h2).^2 + (blk_y - v2).^2 < loop2)
        break
    end

    if ih2 ~= cur_ih || iv2 ~= cur_iv
        n_left = n_left + 1;
        slot = mod(n_left - 1, recent_n) + 1;
        if n_left > recent_n
            visit_stamp(buf_v(slot), buf_h(slot)) = ray_id;
        end
        buf_h(slot) = cur_ih;
        buf_v(slot) = cur_iv;
        cur_ih = ih2;
        cur_iv = iv2;
        if visit_stamp(iv2, ih2) == ray_id
            break
        end
    end

    h = h2;
    v = v2;
    u1 = u_next;
    w1 = w_next;
    arc = arc + step;
    n = n + 1;
    pts(n, 1) = h;
    pts(n, 2) = v;
end
pts = pts(1:n, :);
end


function hit = has_point_within(bins_x, bins_y, h, v, ih, iv, nCh, nCv, n_pad, thresh2)
%HAS_POINT_WITHIN  True if any stored sample in the (2*n_pad+1)^2 neighbourhood
%  lies strictly within squared distance thresh2 of (h, v).  Checks the
%  home bin first and returns immediately on first match.
hit = false;
li0 = iv + (ih - 1) * nCv;
px = bins_x{li0};
if ~isempty(px)
    py = bins_y{li0};
    if any((px - h).^2 + (py - v).^2 < thresh2)
        hit = true;
        return
    end
end
i_lo = max(1, ih - n_pad);
i_hi = min(nCh, ih + n_pad);
j_lo = max(1, iv - n_pad);
j_hi = min(nCv, iv + n_pad);
for ii = i_lo:i_hi
    col_off = (ii - 1) * nCv;
    for jj = j_lo:j_hi
        if ii == ih && jj == iv
            continue
        end
        li = jj + col_off;
        px = bins_x{li};
        if isempty(px)
            continue
        end
        py = bins_y{li};
        if any((px - h).^2 + (py - v).^2 < thresh2)
            hit = true;
            return
        end
    end
end
end


function idx = nearest_index(vec, x, dx)
idx = round((x - vec(1)) / dx) + 1;
idx = min(max(idx, 1), numel(vec));
end


function dil = dilate8(mask)
dil = mask;
for di = -1:1
    for dj = -1:1
        if di == 0 && dj == 0
            continue
        end
        r1 = (1 + max(di, 0)):(size(mask, 1) + min(di, 0));
        c1 = (1 + max(dj, 0)):(size(mask, 2) + min(dj, 0));
        dil(r1, c1) = dil(r1, c1) | mask(r1 - di, c1 - dj);
    end
end
end


function s = polyline_length(pts)
if size(pts, 1) < 2
    s = 0;
    return
end
s = sum(hypot(diff(pts(:, 1)), diff(pts(:, 2))));
end
