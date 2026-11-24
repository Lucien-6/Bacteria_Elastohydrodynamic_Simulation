"""Compare the SYM50 simulated gait with the measured gait of cell 1.00 A.

The measurement is the tracked centreline of the left flagellum (side = -1)
of cell 1.00 A in water (Qin et al., Sci. Rep. 5, 9190, 2015; Gutierrez et
al., Dryad doi:10.5061/dryad.sbcc2frj6).  The simulation is the symmetric
closed-loop case output/Chlamy-SB-SYM50, driven at 57.75 Hz by the
two-harmonic fit of the active moment inferred from that recording.

The image axes are not the cell axes: the measured root tangent sits near
81 deg, while the model root is clamped near 13 deg from the anterior axis.
Registration is one rigid placement of the whole simulated flagellum plus
one beat-phase shift, both constant in time and in arc length:

  * translate the root to the origin
  * reflect about the swimming axis if that lowers the error (image y may
    run opposite to the body +Y axis)
  * rotate by one angle so the tangent reference matches the recording
  * shift the beat phase by one fraction of a cycle

Phase 0 of the recording is its first frame, not a stroke event.  The
drive uses that same time origin, and the best phase lag is under one
percent of a cycle, so no shift is applied.  Coloured curves are sampled
at the phases of the measured figure.  Simulated coloured curves are
linearly interpolated between the 1 ms output frames.  Grey curves are
the other frames of the first measured beat, and the simulated
centreline at those same phases, so both shape panels have the same
curves.

Author: Lucien (lucien-6@qq.com)
Created: 2026-11-14
"""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import plot_measured_gait_1p00A as meas

ROOT = Path(__file__).resolve().parents[2]
SIM = ROOT / "output" / "Chlamy-SB-SYM50" / "simulation_results.mat"
OUT = ROOT / "output" / "chlamy_gait" / "cell_1p00A_sym50_gait_comparison"
F_SIM_HZ = 57.75
N_LATE_BEATS = 8
S_GRID = np.linspace(0.02, 0.98, 48)


def quat_to_matrix(q):
    """Active body-to-lab rotation, scalar-first quaternion."""
    q = np.asarray(q, dtype=float)
    q = q / np.linalg.norm(q)
    q0, q1, q2, q3 = q
    return np.array([
        [1 - 2 * (q2 * q2 + q3 * q3), 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
        [2 * (q1 * q2 + q0 * q3), 1 - 2 * (q1 * q1 + q3 * q3), 2 * (q2 * q3 - q0 * q1)],
        [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), 1 - 2 * (q1 * q1 + q2 * q2)],
    ])


def deref(handle, name):
    return [np.array(handle[ref]) for ref in np.array(handle[name][()]).ravel()]


def load_simulation(path):
    """Body-frame centreline of each flagellum, root at the origin.

    Returns X, Y with shape (nodes, times) and the segment tangent angle
    alpha = atan2(dy, dx), in radians, shape (segments, times).
    """
    with h5py.File(path, "r") as handle:
        raw = deref(handle, "results/X_pili_um")
        centre = np.array(handle["results/X_cm_um"])
        quat = np.array(handle["results/q"])
        time = np.array(handle["results/t"]).ravel()
    curves = []
    for flagellum in raw:
        n_time = flagellum.shape[0]
        nodes = flagellum.shape[2]
        x_body = np.empty((nodes, n_time))
        y_body = np.empty((nodes, n_time))
        z_body = np.empty((nodes, n_time))
        for i in range(n_time):
            points = flagellum[i].T
            local = (quat_to_matrix(quat[i]).T @ (points - centre[i]).T).T
            x_body[:, i] = local[:, 0] - local[0, 0]
            y_body[:, i] = local[:, 1] - local[0, 1]
            z_body[:, i] = local[:, 2] - local[0, 2]
        curves.append({
            "x": x_body,
            "y": y_body,
            "z_abs_max": float(np.max(np.abs(z_body))),
            "alpha": np.arctan2(np.diff(y_body, axis=0), np.diff(x_body, axis=0)),
            "s_frac": arc_fraction(x_body, y_body),
            "length_um": float(np.mean(np.sum(
                np.hypot(np.diff(x_body, axis=0), np.diff(y_body, axis=0)), axis=0))),
        })
    return curves, time


def arc_fraction(x_nodes, y_nodes):
    """Segment-midpoint s/L.  x_nodes is (nodes, times)."""
    step = np.hypot(np.diff(x_nodes, axis=0), np.diff(y_nodes, axis=0))
    return (np.cumsum(step, axis=0) - 0.5 * step) / np.sum(step, axis=0)


def sample_angle(alpha, s_frac, grid):
    """Unwrap along the flagellum and interpolate onto a common s/L grid."""
    out = np.empty((len(grid), alpha.shape[1]))
    for i in range(alpha.shape[1]):
        out[:, i] = np.interp(grid, s_frac[:, i], np.unwrap(alpha[:, i]))
    return out


def periodic_interp(phase_src, values, phase_query):
    """Interpolate values[:, time] onto phase_query.  Phase is modulo 1."""
    order = np.argsort(phase_src)
    phase = phase_src[order]
    samples = values[:, order]
    keep = np.concatenate([[True], np.diff(phase) > 1e-8])
    phase = phase[keep]
    samples = samples[:, keep]
    phase_ext = np.concatenate([phase - 1.0, phase, phase + 1.0])
    samples_ext = np.concatenate([samples, samples, samples], axis=1)
    return np.vstack([
        np.interp(phase_query, phase_ext, row) for row in samples_ext
    ])


def register(measured_angle, measured_phase, sim_angle, sim_phase, late):
    """Constant reflection, rotation and phase lag minimising tangent error.

    sim_angle and measured_angle are unwrapped along the arc and already
    interpolated onto the same s/L grid.
    """
    phase = sim_phase[late]
    best = None
    lags = np.linspace(0.0, 1.0, 720, endpoint=False)
    for mirror in (1.0, -1.0):
        source = mirror * sim_angle[:, late]
        for lag in lags:
            predicted = periodic_interp(phase, source, np.mod(measured_phase + lag, 1.0))
            residual = np.angle(np.exp(1j * (measured_angle - predicted)))
            offset = np.angle(np.mean(np.exp(1j * residual)))
            error = np.angle(np.exp(1j * (residual - offset)))
            rms = float(np.sqrt(np.mean(error ** 2)))
            if best is None or rms < best["rms_rad"]:
                best = {
                    "rms_rad": rms,
                    "lag": float(lag),
                    "offset_rad": float(offset),
                    "mirror": float(mirror),
                }
    return best


def place(x_nodes, y_nodes, mirror, offset):
    """Reflect and rotate a root-origin centreline into the image frame."""
    if mirror < 0:
        y_nodes = -y_nodes
    cosine = np.cos(offset)
    sine = np.sin(offset)
    return cosine * x_nodes - sine * y_nodes, sine * x_nodes + cosine * y_nodes


def resample_nodes(x_nodes, y_nodes, grid):
    """Interpolate node coordinates onto fractional arc lengths."""
    step = np.hypot(np.diff(x_nodes, axis=0), np.diff(y_nodes, axis=0))
    length = np.sum(step, axis=0)
    s_nodes = np.vstack([np.zeros((1, x_nodes.shape[1])), np.cumsum(step, axis=0)])
    s_nodes = s_nodes / length
    out_x = np.empty((len(grid), x_nodes.shape[1]))
    out_y = np.empty((len(grid), y_nodes.shape[1]))
    for i in range(x_nodes.shape[1]):
        out_x[:, i] = np.interp(grid, s_nodes[:, i], x_nodes[:, i])
        out_y[:, i] = np.interp(grid, s_nodes[:, i], y_nodes[:, i])
    return out_x, out_y


def shape_rms(measured_x, measured_y, sim_x, sim_y):
    """RMS point distance (µm) at equal fractional arc length."""
    grid = np.linspace(0.0, 1.0, 64)
    mx, my = resample_nodes(measured_x, measured_y, grid)
    sx, sy = resample_nodes(sim_x, sim_y, grid)
    return float(np.sqrt(np.mean((mx - sx) ** 2 + (my - sy) ** 2)))


def draw_shapes(ax, grey, coloured):
    for x_line, y_line in grey:
        ax.plot(x_line, y_line, color="#E4E4E4", lw=0.32,
                solid_capstyle="round", zorder=1)
    for color, linestyle, (x_line, y_line) in zip(meas.PHASE_COLOR, meas.PHASE_LS, coloured):
        ax.plot(x_line, y_line, color="white", ls=linestyle, lw=1.7,
                solid_capstyle="round", zorder=2)
        ax.plot(x_line, y_line, color=color, ls=linestyle, lw=0.9,
                solid_capstyle="round", zorder=3)
    ax.plot(0.0, 0.0, linestyle="none", marker="o", color="black",
            ms=2.2, mew=0.0, zorder=4)


def style_axis(ax):
    ax.tick_params(direction="in", top=True, right=True, length=2.5, width=0.6)


def panel_letter(ax, letter):
    """Upper-left of the subplot, outside the coordinate frame."""
    ax.annotate(
        letter, xy=(0, 1), xycoords="axes fraction",
        xytext=(-14, 3), textcoords="offset points",
        fontsize=8, fontweight="bold", fontfamily="Arial",
        ha="right", va="bottom", annotation_clip=False,
    )


def panel_note(ax, note):
    """Top-right, just above the axes, on the same line as the panel letter."""
    ax.annotate(
        note, xy=(1, 1), xycoords="axes fraction",
        xytext=(0, 3), textcoords="offset points",
        fontsize=7, fontfamily="Arial", color="#222222",
        ha="right", va="bottom", annotation_clip=False,
    )


def arrange_figure(fig, ax_m, ax_s, ax_pm, ax_ps, legend):
    """Pack the four panels and the legend, in millimetres.

    Shape panels keep equal data scales and share their left edge with the
    angle panel below.  The angle panels are only as wide as the block that
    ends at the legend's right edge, so the gap beside panel b stays short.
    """
    dx = abs(ax_m.get_xlim()[1] - ax_m.get_xlim()[0])
    dy = abs(ax_m.get_ylim()[1] - ax_m.get_ylim()[0])
    shape_h = 62.0
    shape_w = shape_h * dx / dy
    gap_ab = 10.0
    gap_leg = 3.5
    # Wide enough for panel letter d, which sits just left of panel d.
    gap_cd = 11.0
    angle_h = 50.0
    left = 17.0
    bottom = 11.0
    top = 7.0
    row_gap = 12.0
    right_pad = 1.5
    fig_h = bottom + angle_h + row_gap + shape_h + top
    # Wide enough for a first measurement of the legend.
    fig_w = left + 2.0 * shape_w + gap_ab + gap_leg + 45.0 + right_pad
    fig.set_size_inches(fig_w / 25.4, fig_h / 25.4)
    fig.set_layout_engine("none")

    def box(x, y, w, h, width_mm=fig_w):
        return [x / width_mm, y / fig_h, w / width_mm, h / fig_h]

    y_shape = bottom + angle_h + row_gap
    ax_m.set_position(box(left, y_shape, shape_w, shape_h))
    b_left = left + shape_w + gap_ab
    ax_s.set_position(box(b_left, y_shape, shape_w, shape_h))
    legend.set_bbox_to_anchor((
        (b_left + shape_w + gap_leg) / fig_w,
        (y_shape + 0.5 * shape_h) / fig_h,
    ))
    fig.canvas.draw()
    extent = legend.get_window_extent(fig.canvas.get_renderer())
    leg_right = extent.transformed(fig.transFigure.inverted()).x1 * fig_w
    span = leg_right - left
    cd_w = (span - gap_cd) / 2.0
    ax_pm.set_position(box(left, bottom, cd_w, angle_h))
    ax_ps.set_position(box(left + cd_w + gap_cd, bottom, cd_w, angle_h))
    new_w = leg_right + right_pad
    scale = fig_w / new_w
    for ax in (ax_m, ax_s, ax_pm, ax_ps):
        pos = ax.get_position()
        ax.set_position([pos.x0 * scale, pos.y0, pos.width * scale, pos.height])
    legend.set_bbox_to_anchor((
        (b_left + shape_w + gap_leg) / new_w,
        (y_shape + 0.5 * shape_h) / fig_h,
    ))
    fig.set_size_inches(new_w / 25.4, fig_h / 25.4)


def build_figure(measured, simulated, frames, sim_coloured, sim_grey, angles, measured_grey):
    meas.apply_style()
    fig, axes = plt.subplots(
        2, 2, figsize=(7.2, 6.2), sharex="row", sharey="row",
    )
    ax_m, ax_s = axes[0]
    ax_pm, ax_ps = axes[1]

    grey_m = [
        (measured["x"][:, i], measured["y"][:, i])
        for i in measured_grey
    ]
    coloured_m = [
        (measured["x"][:, i], measured["y"][:, i]) for i in frames
    ]
    draw_shapes(ax_m, grey_m, coloured_m)
    draw_shapes(ax_s, sim_grey, sim_coloured)

    drawn = grey_m + coloured_m + list(sim_grey) + list(sim_coloured)
    span_x = np.concatenate([np.asarray(xy[0]) for xy in drawn])
    span_y = np.concatenate([np.asarray(xy[1]) for xy in drawn])
    pad = 0.55
    ax_m.set_xlim(np.min(span_x) - pad, np.max(span_x) + pad)
    ax_m.set_ylim(np.min(span_y) - pad, np.max(span_y) + pad)
    ax_m.set_aspect("equal", adjustable="box")
    ax_s.set_aspect("equal", adjustable="box")
    ax_m.set_xlabel(r"$x$ (µm)")
    ax_s.set_xlabel(r"$x$ (µm)")
    ax_m.set_ylabel(r"$y$ (µm)")

    for color, linestyle, frame in zip(meas.PHASE_COLOR, meas.PHASE_LS, frames):
        ax_pm.plot(
            measured["s"][:, frame], np.rad2deg(measured["phi"][:, frame]),
            color=color, ls=linestyle, lw=0.9,
        )
    for color, linestyle, angle, s_line in zip(
        meas.PHASE_COLOR, meas.PHASE_LS, angles["sim"], angles["s_sim"],
    ):
        ax_ps.plot(s_line, np.rad2deg(angle), color=color, ls=linestyle, lw=0.9)

    phi_drawn = [np.rad2deg(measured["phi"][:, frame]) for frame in frames]
    phi_drawn.extend(np.rad2deg(angle) for angle in angles["sim"])
    phi_all = np.concatenate(phi_drawn)
    pad_phi = 0.04 * (np.max(phi_all) - np.min(phi_all))
    ax_pm.set_xlim(0.0, 1.0)
    ax_pm.set_ylim(np.min(phi_all) - pad_phi, np.max(phi_all) + pad_phi)
    ax_pm.set_xlabel(r"Relative arc length $s/L$")
    ax_ps.set_xlabel(r"Relative arc length $s/L$")
    ax_pm.set_ylabel(r"Tangent angle $\it{φ}$ (deg)")

    for ax in axes.ravel():
        style_axis(ax)
    for ax, letter in zip(axes.ravel(), "abcd"):
        panel_letter(ax, letter)

    for ax, note in (
        (ax_m, "Measured"),
        (ax_s, "Simulation"),
        (ax_pm, "Measured"),
        (ax_ps, "Simulation"),
    ):
        panel_note(ax, note)
    handles = []
    for color, linestyle, frame in zip(meas.PHASE_COLOR, meas.PHASE_LS, frames):
        handles.append(Line2D(
            [0], [0], color=color, ls=linestyle, lw=0.9,
            label=f"{measured['phase'][frame]:.2f}",
        ))
    handles.append(Line2D(
        [0], [0], color="#E4E4E4", lw=0.9, label="Other frames",
    ))
    legend = fig.legend(
        handles=handles, title="Beat phase", ncol=1,
        frameon=False, loc="center left",
        bbox_to_anchor=(0.8, 0.7), bbox_transform=fig.transFigure,
        borderaxespad=0.0, handlelength=2.2, labelspacing=0.35,
        title_fontsize=7,
    )
    legend._legend_box.align = "left"
    arrange_figure(fig, ax_m, ax_s, ax_pm, ax_ps, legend)
    return fig


def main():
    x_m, y_m, phi_m, s_m, phase_m, meta = meas.load_gait(meas.MAT)
    pool = meas.first_beat(phase_m)
    _, frames = meas.phase_samples(phase_m, n_phase=6, pool=pool)
    frames = [int(i) for i in frames]
    measured = {
        "x": x_m, "y": y_m, "phi": phi_m, "s": s_m, "phase": phase_m,
    }
    phi_grid = sample_angle(phi_m, s_m, S_GRID)

    curves, time = load_simulation(SIM)
    phase_s = np.mod(time * F_SIM_HZ, 1.0)
    late = time >= N_LATE_BEATS / F_SIM_HZ
    # Fit on the last eight beats, after the start from the resting shape.
    n_keep = int(round(N_LATE_BEATS / F_SIM_HZ / np.median(np.diff(time))))
    late_idx = np.flatnonzero(late)[-n_keep:]
    late_mask = np.zeros(time.shape, dtype=bool)
    late_mask[late_idx] = True

    choice = None
    for index, curve in enumerate(curves):
        alpha_grid = sample_angle(curve["alpha"], curve["s_frac"], S_GRID)
        fit = register(phi_grid, phase_m, alpha_grid, phase_s, late_mask)
        fit["flagellum"] = index + 1
        fit["alpha_grid"] = alpha_grid
        fit["curve"] = curve
        print(
            f"flagellum {index + 1}: mirror {fit['mirror']:+.0f}, "
            f"lag {fit['lag']:.3f} cycles, offset {np.rad2deg(fit['offset_rad']):.1f} deg, "
            f"tangent RMS {np.rad2deg(fit['rms_rad']):.2f} deg, "
            f"L {curve['length_um']:.2f} um, max |z| {curve['z_abs_max']:.2e} um"
        )
        if choice is None or fit["rms_rad"] < choice["rms_rad"]:
            choice = fit

    mirror = choice["mirror"]
    # A lag within 1% of a cycle is the same time origin as the recording.
    if min(choice["lag"], 1.0 - choice["lag"]) < 0.01:
        lag = 0.0
    else:
        lag = choice["lag"]
    source = mirror * choice["alpha_grid"][:, late_mask]
    predicted = periodic_interp(
        np.mod(phase_s[late_mask] - lag, 1.0), source, phase_m,
    )
    residual = np.angle(np.exp(1j * (phi_grid - predicted)))
    offset = float(np.angle(np.mean(np.exp(1j * residual))))
    error = np.angle(np.exp(1j * (residual - offset)))
    rms_tangent = float(np.sqrt(np.mean(error ** 2)))
    rms_along = np.rad2deg(np.sqrt(np.mean(error ** 2, axis=1)))

    curve = choice["curve"]
    x_r, y_r = place(curve["x"], curve["y"], mirror, offset)
    alpha_along = np.unwrap(curve["alpha"], axis=0)
    mid = alpha_along.shape[0] // 2
    branch = np.round((alpha_along[mid] - alpha_along[mid, 0]) / (2.0 * np.pi))
    alpha_along = alpha_along - 2.0 * np.pi * branch
    alpha_r = mirror * alpha_along + offset
    s_r = curve["s_frac"]

    targets = phase_m[frames]
    sim_phase = np.mod(phase_s[late_mask] - lag, 1.0)
    coloured = []
    sim_angles = []
    sim_s = []
    x_at, y_at = [], []
    for target in targets:
        x_line = periodic_interp(sim_phase, x_r[:, late_mask], np.array([target]))[:, 0]
        y_line = periodic_interp(sim_phase, y_r[:, late_mask], np.array([target]))[:, 0]
        a_line = periodic_interp(sim_phase, alpha_r[:, late_mask], np.array([target]))[:, 0]
        s_line = np.mean(s_r[:, late_mask], axis=1)
        coloured.append((x_line, y_line))
        sim_angles.append(a_line)
        sim_s.append(s_line)
        x_at.append(x_line)
        y_at.append(y_line)

    # Phase-matched centreline at every measured frame, for the shape error.
    x_match = periodic_interp(sim_phase, x_r[:, late_mask], phase_m)
    y_match = periodic_interp(sim_phase, y_r[:, late_mask], phase_m)
    # load_gait subtracts the mean base; remove the remaining per-frame base.
    x_meas = x_m - x_m[0]
    y_meas = y_m - y_m[0]
    rms_shape = shape_rms(x_meas, y_meas, x_match, y_match)

    highlight_x = np.column_stack(x_at)
    highlight_y = np.column_stack(y_at)
    rms_highlight = shape_rms(
        x_meas[:, frames], y_meas[:, frames], highlight_x, highlight_y,
    )

    # One beat on each side.  The recording has several beats at 1.65 ms,
    # and drawing all of them made panel a much denser than one simulated
    # beat saved at 1 ms.  The coloured phases already come from the first
    # measured beat; the grey curves are its remaining frames, and the
    # simulation is drawn at those same phases.
    measured_grey = [int(i) for i in pool if int(i) not in set(frames)]
    grey = []
    for target in phase_m[measured_grey]:
        x_line = periodic_interp(sim_phase, x_r[:, late_mask], np.array([target]))[:, 0]
        y_line = periodic_interp(sim_phase, y_r[:, late_mask], np.array([target]))[:, 0]
        grey.append((x_line, y_line))
    envelope = [alpha_r[:, i] for i in late_idx]

    root_m = np.rad2deg(phi_m[0])
    root_s = np.rad2deg(alpha_r[0, late_mask])
    symmetry = float(np.max(np.abs(curves[0]["y"] + curves[1]["y"])))
    print(
        f"selected flagellum {choice['flagellum']}  mirror {mirror:+.0f}  "
        f"lag {lag:.3f} cycles  offset {np.rad2deg(offset):.1f} deg"
    )
    print(
        f"tangent RMS {np.rad2deg(rms_tangent):.2f} deg  "
        f"shape RMS {rms_shape:.2f} um over all frames, "
        f"{rms_highlight:.2f} um at the six phases"
    )
    print(
        "tangent RMS along s/L: "
        + ", ".join(f"{s:.2f}:{r:.1f}" for s, r in zip(S_GRID[::8], rms_along[::8]))
    )
    print(
        f"root angle measured {root_m.mean():.1f} ± {root_m.std():.1f} deg "
        f"[{root_m.min():.1f}, {root_m.max():.1f}], "
        f"simulated {root_s.mean():.1f} ± {root_s.std():.1f} deg"
    )
    print(
        f"mirror residual max |y1+y2| {symmetry:.3f} um  "
        f"measured f {meta['f_Hz']:.2f} Hz, L {meta['L_um']:.2f} um, "
        f"n {meta['n_frame']}, side {meta['side']}, "
        f"grey curves {len(measured_grey)} on each side, "
        f"plus {len(frames)} coloured phases"
    )

    fig = build_figure(
        measured, curve, frames, coloured, grey,
        {"sim": sim_angles, "s_sim": sim_s, "sim_envelope": envelope},
        measured_grey,
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT.with_suffix(".pdf"))
    fig.savefig(OUT.with_suffix(".svg"))
    fig.savefig(OUT.with_suffix(".png"), dpi=600)
    preview = OUT.parent / "_preview_sym50_comparison.png"
    fig.savefig(preview, dpi=150)
    print(f"wrote {OUT.with_suffix('.pdf')}")
    return fig


if __name__ == "__main__":
    main()
