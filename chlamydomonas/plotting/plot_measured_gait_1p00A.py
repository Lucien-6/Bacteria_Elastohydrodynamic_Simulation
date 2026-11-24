"""Measured flagellar gait of cell 1p00 A.

Plots the tracked centerline and tangent angle from the Dryad archive
doi:10.5061/dryad.sbcc2frj6 (Qin et al., Sci. Rep. 5, 9190, 2015;
redistributed by Gutierrez et al.).  Coordinates are the imaging plane
in micrometres, translated so that the time-averaged base is at the
origin.  Phase 0 is the first recorded frame, not a defined stroke event.

Author: Lucien (lucien-6@qq.com)
Created: 2026-11-12
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
MAT = (
    ROOT
    / "chlamydomonas"
    / "Refers"
    / "doi_10_5061_dryad_sbcc2frj6__v20250313"
    / "Quantification_of_Flagellar_Gait_Changes_Gutierrez_Data"
    / "Data"
    / "Newtonian 1p00 A_Viscosity1.00.mat"
)
OUT = ROOT / "output" / "chlamy_gait" / "cell_1p00A_measured_gait"

# Dark Okabe–Ito / Tol colours.  Each phase also has its own dash pattern
# so the curves stay distinct in grayscale.  Light yellow and sky blue are
# omitted: both disappear against the grey beat envelope.
PHASE_COLOR = ["#000000", "#0072B2", "#009E73", "#5E2A84", "#B35C00", "#9A2B0F"]
PHASE_LS = [
    "-",
    (0, (4.2, 1.3)),
    (0, (2.2, 0.95)),
    (0, (5.0, 1.15, 1.05, 1.15)),
    (0, (1.15, 0.75)),
    (0, (3.4, 1.05)),
]
WIDTH_MM = 183.0


def apply_style():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.grid": False,
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "savefig.dpi": 600,
        "figure.dpi": 150,
    })


def beat_frequency(phi, dt):
    """Dominant frequency of the mid-flagellum tangent angle, in Hz."""
    mid = phi.shape[0] // 2
    signal = phi[mid] - np.mean(phi[mid])
    time = np.arange(phi.shape[1]) * dt
    freq = np.arange(30.0, 90.0 + 0.005, 0.01)
    cis = np.exp(-2j * np.pi * freq[:, None] * time[None, :])
    amp = np.abs(cis @ signal)
    return float(freq[int(np.argmax(amp))])


def load_gait(path):
    import scipy.io

    raw = scipy.io.loadmat(path, squeeze_me=True)
    scale = float(raw["space_scale"])
    dt = float(raw["time_scale"])
    xx = np.asarray(raw["XX"], dtype=float)
    yy = np.asarray(raw["YY"], dtype=float)
    phi = np.asarray(raw["phi"], dtype=float)
    xx = (xx - np.mean(xx[0])) * scale
    yy = (yy - np.mean(yy[0])) * scale
    ds = np.hypot(np.diff(xx, axis=0), np.diff(yy, axis=0))
    length = np.sum(ds, axis=0)
    s_mid = np.cumsum(ds, axis=0) - 0.5 * ds
    s_frac = s_mid / length
    freq = beat_frequency(phi, dt)
    phase = np.mod(np.arange(phi.shape[1]) * dt * freq, 1.0)
    meta = {
        "f_Hz": freq,
        "dt_s": dt,
        "L_um": float(raw["flagella_length"]) * scale,
        "L_frame_um": length,
        "atch_deg": float(np.rad2deg(raw["atch_angle"])),
        "side": int(raw["side"]),
        "n_frame": int(phi.shape[1]),
    }
    return xx, yy, phi, s_frac, phase, meta


def first_beat(phase):
    """Frames from phase 0 up to and including the last sample before the wrap."""
    wraps = np.where(np.diff(phase) < -0.5)[0]
    if len(wraps) == 0:
        return np.arange(len(phase))
    return np.arange(0, int(wraps[0]) + 1)


def phase_samples(phase, n_phase, pool):
    """Frame in pool nearest each evenly spaced phase, on the circular phase."""
    target = np.arange(n_phase) / n_phase
    sub = phase[pool]
    delta = np.abs(sub[None, :] - target[:, None])
    delta = np.minimum(delta, 1.0 - delta)
    order = np.argsort(delta, axis=1)
    chosen = []
    used = set()
    for row in order:
        for j in row:
            frame = int(pool[int(j)])
            if frame not in used:
                chosen.append(frame)
                used.add(frame)
                break
    return target, np.asarray(chosen, dtype=int)


def build_figure(xx, yy, phi, s_frac, phase, frames):
    apply_style()
    width = WIDTH_MM / 25.4
    fig, axes = plt.subplots(
        1, 2, figsize=(width, 4.45), layout="constrained",
    )
    ax_xy, ax_phi = axes
    highlight = set(int(i) for i in frames)

    for frame in range(xx.shape[1]):
        if frame in highlight:
            continue
        ax_xy.plot(xx[:, frame], yy[:, frame], color="#E4E4E4",
                   lw=0.32, solid_capstyle="round", zorder=1)

    handles = [
        Line2D([0], [0], color="#C2C2C2", lw=0.8, label="Other frames"),
    ]
    for color, linestyle, frame in zip(PHASE_COLOR, PHASE_LS, frames):
        ax_xy.plot(xx[:, frame], yy[:, frame], color="white",
                   ls=linestyle, lw=1.7, solid_capstyle="round", zorder=2)
        ax_xy.plot(xx[:, frame], yy[:, frame], color=color,
                   ls=linestyle, lw=0.9, solid_capstyle="round", zorder=3)
        ax_phi.plot(s_frac[:, frame], np.rad2deg(phi[:, frame]),
                    color=color, ls=linestyle, lw=0.9, zorder=3)
        handles.append(Line2D(
            [0], [0], color=color, ls=linestyle, lw=0.9,
            label=f"{phase[frame]:.2f}",
        ))

    ax_xy.plot(0.0, 0.0, linestyle="none", marker="o", color="black",
               ms=2.2, mew=0.0, zorder=4)
    pad = 0.55
    ax_xy.set_xlim(np.min(xx) - pad, np.max(xx) + pad)
    ax_xy.set_ylim(np.min(yy) - pad, np.max(yy) + pad)
    ax_xy.set_aspect("equal", adjustable="box")
    ax_xy.set_xlabel(r"$x$ (µm)")
    ax_xy.set_ylabel(r"$y$ (µm)")

    ax_phi.set_xlim(0.0, 1.0)
    phi_deg = np.rad2deg(phi)
    pad_phi = 0.04 * (np.max(phi_deg) - np.min(phi_deg))
    ax_phi.set_ylim(np.min(phi_deg) - pad_phi, np.max(phi_deg) + pad_phi)
    ax_phi.set_xlabel(r"Relative arc length $s/L$")
    ax_phi.set_ylabel(r"Tangent angle $\it{φ}$ (deg)")

    for ax in axes:
        ax.tick_params(direction="in", top=True, right=True, length=2.5, width=0.6)

    for ax, letter in zip(axes, "ab"):
        ax.annotate(
            letter, xy=(0, 1), xycoords="axes fraction",
            xytext=(-14, 3), textcoords="offset points",
            fontsize=8, fontweight="bold", fontfamily="Arial",
            ha="right", va="bottom", annotation_clip=False,
        )

    fig.legend(
        handles=handles, title="Beat phase", ncol=4,
        frameon=False, loc="outside lower center",
        borderaxespad=0.0, handlelength=2.4, columnspacing=1.0,
    )
    return fig


def main():
    xx, yy, phi, s_frac, phase, meta = load_gait(MAT)
    pool = first_beat(phase)
    _, frames = phase_samples(phase, n_phase=6, pool=pool)
    print(
        f"f_Hz={meta['f_Hz']:.2f}  dt_ms={meta['dt_s'] * 1e3:.3f}  "
        f"frames={meta['n_frame']}  L_um={meta['L_um']:.2f}  "
        f"atch_deg={meta['atch_deg']:.2f}  side={meta['side']}"
    )
    print("phase samples:", ", ".join(
        f"frame {int(i)} phase {phase[i]:.3f}" for i in frames
    ))
    fig = build_figure(xx, yy, phi, s_frac, phase, frames)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT.with_suffix(".pdf"))
    fig.savefig(OUT.with_suffix(".svg"))
    fig.savefig(OUT.with_suffix(".png"), dpi=600)
    preview = OUT.parent / "_preview_cell_1p00A.png"
    fig.savefig(preview, dpi=150)
    print(f"wrote {OUT.with_suffix('.pdf')}")
    return fig


if __name__ == "__main__":
    main()
