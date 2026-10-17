function add_gravity_indicator_outside(fig_handle, ax_handle, font_size, line_width)
%ADD_GRAVITY_INDICATOR_OUTSIDE  Draw gravity arrow flush outside axis right edge.
%  Uses figure-normalised annotation coordinates so marker location is
%  stable across views and axis limits.
%
%  Root-cause fix (v2.4):
%    With DataAspectRatio=[1 1 1] and a strongly elongated data range
%    (e.g. Z >> X,Y for a sinking bacterium), MATLAB renders a very narrow
%    axis box CENTRED inside the Position rectangle.  Using the Position
%    right edge (the old approach) places the arrow hundreds of pixels to
%    the right of the actual box edge.
%
%    This function instead computes the TRUE rendered right edge by:
%      1. Obtaining the figure size in pixels (aspect ratio).
%      2. Projecting the data bounding box (xlim x ylim x zlim) onto the
%         screen-horizontal direction via the orthographic approximation:
%           right   = [cos(az), sin(az), 0]
%           data_w  = dx|cos az| + dy|sin az|
%           data_h  = (dx|sin az|+dy|cos az|)|sin el| + dz|cos el|
%      3. Applying the equal-aspect scale:
%           px_per_unit = min(ax_w_px/data_w, ax_h_px/data_h)
%      4. Computing the centred right edge:
%           actual_right = ax_pos(1) + ax_pos(3)/2
%                        + data_w*px_per_unit / (2*fig_w_px)
%
%  IMPORTANT:
%    1. Must be called AFTER any sgtitle() on the parent figure.
%    2. Must be called AFTER the axes is fully configured (xlim/ylim/zlim,
%       DataAspectRatio, view, labels, title) so ax.Position is settled.
%
%  Inputs:
%    fig_handle - Parent figure handle
%    ax_handle  - 3-D axes handle whose rendered right edge is the anchor
%    font_size  - Font size of the "G" label [pt]
%    line_width - Line width of the arrow shaft [pt]
%
%  Outputs:
%    None.  Two annotation objects (arrow and "G" label) are added to the figure.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-17  (extracted verbatim from the local function of the same
%           name in run_simulation.m / run_simulation_HPC_Linux.m)
%  Updated: 2026-10-17

    drawnow;   % Flush pending layout so ax.Position is up-to-date

    ax_pos = get(ax_handle, 'Position');   % [left bottom width height] normalised

    % ---- Figure size in pixels (needed for aspect-ratio computation) ----
    old_units = get(fig_handle, 'Units');
    set(fig_handle, 'Units', 'pixels');
    fig_pos_px = get(fig_handle, 'Position');   % [x y w h] in screen pixels
    set(fig_handle, 'Units', old_units);
    fig_w_px = fig_pos_px(3);
    fig_h_px = fig_pos_px(4);

    % ---- Axis Position in pixels ----------------------------------------
    ax_w_px = ax_pos(3) * fig_w_px;
    ax_h_px = ax_pos(4) * fig_h_px;

    % ---- Data bounding-box ranges (equal physical units on all axes) ----
    dx = abs(diff(xlim(ax_handle)));
    dy = abs(diff(ylim(ax_handle)));
    dz = abs(diff(zlim(ax_handle)));

    % ---- Orthographic projection of the data box onto screen axes -------
    %  Screen right-vector:  r = [cos(az), sin(az), 0]
    %  Screen up-vector:     u = [-sin(az)*sin(el), cos(az)*sin(el), cos(el)]
    %
    %  Projected extents of the axis-aligned bounding box:
    %    data_w = dx|cos az| + dy|sin az|
    %    data_h = (dx|sin az|+dy|cos az|)*|sin el| + dz*|cos el|
    [az_deg, el_deg] = view(ax_handle);
    ca = abs(cos(az_deg * pi/180));   sa = abs(sin(az_deg * pi/180));
    ce = abs(cos(el_deg * pi/180));   se = abs(sin(el_deg * pi/180));

    data_w = dx * ca + dy * sa;
    data_h = (dx * sa + dy * ca) * se + dz * ce;
    data_w = max(data_w, 1e-12);   % guard against degenerate limits
    data_h = max(data_h, 1e-12);

    % ---- Equal-aspect scale factor: whichever axis fills first ----------
    px_per_unit = min(ax_w_px / data_w,  ax_h_px / data_h);

    % ---- Actual rendered half-width in normalised figure coords ---------
    % MATLAB centres the equal-aspect box inside the Position rectangle.
    actual_half_w = (data_w * px_per_unit) / (2 * fig_w_px);
    actual_right  = ax_pos(1) + ax_pos(3) / 2 + actual_half_w;

    % Place arrow 1 % of figure width to the right of the rendered box edge
    x_arrow = actual_right + 0.010;
    y_mid   = ax_pos(2) + 0.50 * ax_pos(4);
    y_top   = y_mid + 0.055;
    y_bot   = y_mid - 0.055;

    % Keep annotation within figure bounds
    x_arrow = min(max(x_arrow, 0.02), 0.98);
    y_top   = min(max(y_top,   0.05), 0.98);
    y_bot   = min(max(y_bot,   0.02), 0.95);
    if y_top <= y_bot
        y_top = min(0.98, y_bot + 0.06);
    end

    annotation(fig_handle, 'arrow', [x_arrow, x_arrow], [y_top, y_bot], ...
        'Color', [0.25 0.25 0.25], 'LineWidth', line_width, 'HeadLength', 12, ...
        'HeadWidth', 10);
    % "G" label centred directly below the arrow shaft
    lbl_w = 0.030;
    lbl_h = 0.030;
    annotation(fig_handle, 'textbox', ...
        [x_arrow - lbl_w/2, y_bot - lbl_h - 0.004, lbl_w, lbl_h], ...
        'String', '{\itG}', 'Interpreter', 'tex', 'FontName', 'Arial', ...
        'FontSize', font_size, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');
end
