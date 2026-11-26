"""Redraw result figures at Nature print size.

The tables under ``output/paper/results`` and the saved inverse-fit and
thermal-acceptance records are plotted again at 89, 120 or 183 mm, in Arial,
with 7 pt axis labels, 6.5 pt ticks and 8 pt panel letters.  No new
simulation results are computed.  Where a second series was a closed-form
curve, it is evaluated from the same expression the captions already state.

Author: Lucien (lucien-6@qq.com)
Created: 2026-11-14
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib._mathtext as _mathtext
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "output" / "paper" / "results"
FIG = PAPER / "figures"
QA = ROOT / "output" / "_nature_preview"

MM = 1.0 / 25.4
# Elsevier CPC column widths: single (90 mm), 1.5-column (140 mm), double (190 mm).
W90 = 90.0 * MM
W140 = 140.0 * MM
W190 = 190.0 * MM
# Backward-compatible aliases for scripts importing legacy width names.
W89 = W90
W120 = W140
W183 = W190

# Okabe–Ito, omitting yellow (it disappears on white).
BLUE = "#0072B2"
SKY = "#56B4E9"
ORANGE = "#E69F00"
VERMILION = "#D55E00"
GREEN = "#009E73"
PINK = "#CC79A7"
BLACK = "#000000"
GREY = "#8A8A8A"
REF = "#B3B3B3"

# Dark to light blue for the ordered body-resolution series in panel c.
BLUE_SEQ = ("#003C67", "#0072B2", "#56B4E9")


def apply_style() -> None:
    """Arial, inward ticks, full frame.  Sizes are the printed point sizes.

    Mathtext draws subscripts at SHRINK_FACTOR times the parent size.  At
    0.7 the subscript of a 7 pt label is 4.9 pt, under Nature's 5 pt floor.
    0.82 keeps a 7 pt subscript at 5.7 pt, and every glyph stays in Arial.
    """
    _mathtext.SHRINK_FACTOR = 0.82
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.minor.size": 1.4,
        "ytick.minor.size": 1.4,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.grid": False,
        "axes.axisbelow": True,
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "savefig.dpi": 600,
        "figure.dpi": 120,
        "legend.frameon": False,
    })


def polish(ax: plt.Axes) -> None:
    ax.tick_params(direction="in", top=True, right=True, length=2.4, width=0.5)
    ax.margins(x=0.04, y=0.06)


def panel_letters(axes, letters: str) -> None:
    """Place unified Elsevier CPC panel labels (a), (b), ... at panel top-left."""
    for ax, letter in zip(np.ravel(axes), letters):
        ax.annotate(
            f"({letter.lower()})",
            xy=(0, 1),
            xycoords="axes fraction",
            xytext=(-14, 2),
            textcoords="offset points",
            fontsize=8,
            fontweight="bold",
            fontfamily="Arial",
            ha="right",
            va="bottom",
            annotation_clip=False,
        )


def reference_line(ax: plt.Axes, value: float, axis: str = "h") -> None:
    if axis == "h":
        ax.axhline(value, color=REF, lw=0.5, zorder=0)
    else:
        ax.axvline(value, color=REF, lw=0.5, zorder=0)


def save_figure(fig: plt.Figure, path: Path) -> None:
    """Write PDF and 600 dpi PNG at the figure's own inch size."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".png"), dpi=600)


def preview(fig: plt.Figure, name: str) -> None:
    QA.mkdir(parents=True, exist_ok=True)
    fig.savefig(QA / f"{name}.png", dpi=140)
    # Grayscale copy for the colour-blind check.  The data marks stay.
    fig.savefig(QA / f"{name}_gray.png", dpi=140)
    gray = plt.imread(QA / f"{name}_gray.png")
    if gray.ndim == 3:
        tone = gray[..., :3].mean(axis=2)
        plt.imsave(QA / f"{name}_gray.png", tone, cmap="gray")


def percent(rel: pd.Series | np.ndarray) -> np.ndarray:
    return np.asarray(rel, dtype=float) * 100.0


def load_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(PAPER / name)


# ---------------------------------------------------------------------------
# Paper figures
# ---------------------------------------------------------------------------

def fig_stokes() -> plt.Figure:
    """Sphere sedimentation against Stokes' law.  Single column, 89 mm.

    ``static/S1_sphere_stokes.csv`` is not in the results folder.  On the
    published figure the resistance-solve line passes through the
    time-stepped markers, and the caption states that the two agree, so the
    line is the piecewise-linear connection of the tabulated time-stepped
    deviations.
    """
    raw = load_csv("stokes_sphere.csv")
    table = raw[raw["id"].str.startswith("CONV-Nbody-sphere")].copy()
    table = table.sort_values("N_body")
    n_body = table["N_body"].to_numpy()
    dev = (table["V_over_Stokes"].to_numpy() - 1.0) * 100.0

    apply_style()
    fig, ax = plt.subplots(figsize=(W89, 68 * MM), layout="constrained")
    ax.plot(n_body, dev, color=ORANGE, lw=0.9, zorder=2)
    ax.plot(
        n_body, dev, linestyle="none", marker="o", ms=3.6,
        color=BLUE, zorder=3,
    )
    reference_line(ax, 0.0)
    ax.set_xlabel(r"$N_{\mathrm{body}}$")
    ax.set_ylabel(r"$V_{\mathrm{sink}}\,/\,V_{\mathrm{Stokes}}-1$ (%)")
    ax.set_xlim(120, 1280)
    ax.set_ylim(-0.55, 0.55)
    ax.set_xticks([200, 400, 600, 800, 1000, 1200])
    handles = [
        Line2D([0], [0], color=ORANGE, lw=0.9, label="resistance solve"),
        Line2D([0], [0], color=BLUE, marker="o", lw=0, ms=3.6,
               label="time-stepped runs"),
    ]
    ax.legend(handles=handles, loc="upper left", borderaxespad=0.3)
    polish(ax)
    return fig


def fig_arrangement() -> plt.Figure:
    """Filament arrangement.  The ratio axis starts at 0."""
    table = load_csv("sed_arrangement_length.csv")
    series = [
        ("mono", "monotrichous (1 pole)", BLUE, "o", "-"),
        ("amphi", "amphitrichous (2 poles)", ORANGE, "s", "--"),
        ("lat2", "2 lateral", GREEN, "^", "-."),
        ("lat4", "4 lateral", VERMILION, "D", (0, (1.2, 0.7))),
        ("hemi4", "4 asymmetric", SKY, "v", (0, (3.0, 0.8, 0.7, 0.8))),
        ("peri8", "8 around caps", BLACK, "*", (0, (0.6, 0.7))),
    ]
    apply_style()
    fig, ax = plt.subplots(figsize=(W89, 92 * MM), layout="constrained")
    handles = []
    for key, label, color, marker, linestyle in series:
        part = table[table["arrangement"] == key].sort_values("length_um")
        ax.plot(
            part["length_um"], part["V_over_no_filament"],
            color=color, marker=marker, ms=3.5, lw=0.9, ls=linestyle,
        )
        handles.append(Line2D(
            [0], [0], color=color, marker=marker, ms=3.5, lw=0.9,
            ls=linestyle, label=label,
        ))
    reference_line(ax, 1.0)
    handles.append(Line2D([0], [0], color=REF, lw=0.8, label="no filament"))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 1.12)
    ax.set_xticks([1, 3, 5, 10])
    ax.set_xlabel("filament length (µm)")
    ax.set_ylabel(r"$V_{\mathrm{sink}}\,/\,V_{\mathrm{sink}}$ (no filament)")
    fig.legend(
        handles=handles, loc="outside lower center", ncol=2,
        borderaxespad=0.0, handlelength=2.2, labelspacing=0.3, columnspacing=0.8,
    )
    polish(ax)
    ax.margins(x=0, y=0)
    return fig


def fig_brownian() -> plt.Figure:
    """Ensemble ratios.  Error bars are one standard error of the ratio."""
    table = load_csv("brownian_ensembles.csv").set_index("id")
    order = [
        "VAL-sph-brown",
        "VAL-cap-brown",
        "VAL-cap-brown-sed",
        "VAL-dtB-0.0500",
        "VAL-dtB-0.0250",
        "VAL-dtB-0.0100",
        "VAL-dtB-0.0050",
    ]
    labels = [
        "sphere",
        "capsule",
        "+ gravity",
        "0.050 s",
        "0.025 s",
        "0.010 s",
        "0.005 s",
    ]
    rows = table.loc[order]
    x = np.arange(len(order))
    dt = rows["D_T_ratio"].to_numpy()
    dt_se = (rows["D_T_se"] / rows["D_T_pred"]).to_numpy()
    dr = rows["D_perp_ratio"].to_numpy()
    dr_se = (rows["D_perp_se"] / rows["D_perp_pred"]).to_numpy()
    sink_id = "VAL-cap-brown-sed"
    vs = float(table.loc[sink_id, "V_sink_ratio"])
    vs_se = float(table.loc[sink_id, "V_sink_dZ_se"] / table.loc[sink_id, "V_sink_pred_um_s"])

    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(W183, 78 * MM), layout="constrained")
    specs = [
        (axes[0], dt, dt_se, x, labels, r"$D_{\mathrm{t}}\,/\,D_{\mathrm{t,pred}}$"),
        (axes[1], dr, dr_se, x, labels, r"$D_{\mathrm{r}}\,/\,D_{\mathrm{r,pred}}$"),
    ]
    for ax, y, se, xpos, tick_labels, ylabel in specs:
        reference_line(ax, 1.0)
        ax.errorbar(
            xpos, y, yerr=se, fmt="o", color=BLUE, ms=3.2, lw=0,
            elinewidth=0.6, capsize=1.8, capthick=0.6, zorder=3,
        )
        ax.set_xticks(xpos)
        ax.set_xticklabels(tick_labels, rotation=28, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.6, len(xpos) - 0.4)
        span = np.nanmax(y + se) - np.nanmin(y - se)
        ax.set_ylim(np.nanmin(y - se) - 0.35 * span, np.nanmax(y + se) + 0.35 * span)
        ax.set_ylabel(ylabel)
        polish(ax)
    reference_line(axes[2], 1.0)
    axes[2].errorbar(
        [0], [vs], yerr=[vs_se], fmt="o", color=BLUE, ms=3.2, lw=0,
        elinewidth=0.6, capsize=1.8, capthick=0.6,
    )
    axes[2].set_xticks([0])
    axes[2].set_xticklabels(["capsule\n+ gravity"])
    axes[2].set_xlim(-0.8, 0.8)
    axes[2].set_ylim(vs - vs_se - 0.04, vs + vs_se + 0.04)
    axes[2].set_ylabel(r"$V_{\mathrm{sink}}\,/\,V_{\mathrm{sink,pred}}$")
    polish(axes[2])
    panel_letters(axes, "abc")
    return fig


def _curve(ax, x, y, color, marker, linestyle, label: str | None = None) -> None:
    ax.plot(x, y, color=color, marker=marker, ms=3.2, lw=0.9, ls=linestyle, label=label)


def fig_convergence() -> plt.Figure:
    """Six convergence sweeps.  Panel c uses lightness plus marker shape."""
    sphere = load_csv("convergence_N_body_sphere.csv")
    capsule = load_csv("convergence_N_body_cap.csv")
    one = load_csv("convergence_N_body_one.csv")
    seg_one = load_csv("convergence_seg_per_um_one.csv")
    seg_four = load_csv("convergence_seg_per_um_hemi4.csv")
    grid = load_csv("convergence_grid.csv")
    tol = load_csv("convergence_tol.csv")
    eps = load_csv("convergence_eps_over_h.csv")
    beta_one = load_csv("convergence_beta_one.csv")
    beta_four = load_csv("convergence_beta_hemi4.csv")

    def drop_default(frame: pd.DataFrame) -> pd.DataFrame:
        # The captions omit the default N_body = 300 from panel a.
        return frame[frame["sweep_value"] != 300].sort_values("sweep_value")

    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(W183, 118 * MM), layout="constrained")
    ax_a, ax_b, ax_c, ax_d, ax_e, ax_f = axes.ravel()

    for frame, color, marker, linestyle, label in (
        (drop_default(sphere), BLUE, "o", "-", "sphere"),
        (drop_default(capsule), ORANGE, "s", "--", "capsule"),
        (drop_default(one), GREEN, "^", "-.", "capsule + 1 filament"),
    ):
        _curve(ax_a, frame["sweep_value"], percent(frame["V_rel_to_ref"]),
               color, marker, linestyle, label)
    reference_line(ax_a, 0.0)
    ax_a.set_xlabel(r"$N_{\mathrm{body}}$")
    ax_a.set_ylabel(r"$\Delta V_{\mathrm{sink}}\,/\,V_{\mathrm{ref}}$ (%)")
    ax_a.legend(loc="upper right", borderaxespad=0.15, handlelength=2.4)

    _curve(ax_b, seg_one["sweep_value"], percent(seg_one["V_rel_to_ref"]),
           GREEN, "^", "-", "1 filament")
    _curve(ax_b, seg_four["sweep_value"], percent(seg_four["V_rel_to_ref"]),
           PINK, "D", "--", "4 filaments")
    reference_line(ax_b, 0.0)
    ax_b.set_xlabel("segments per µm")
    ax_b.set_ylabel(r"$\Delta V_{\mathrm{sink}}\,/\,V_{\mathrm{ref}}$ (%)")
    ax_b.set_xticks([5, 10, 15, 20, 25, 30])
    ax_b.legend(loc="center right", borderaxespad=0.15)

    grid = grid.copy()
    grid["n_body"] = (grid["sweep_value"] // 100).astype(int)
    grid["seg"] = (grid["sweep_value"] % 100).astype(int)
    for (n_body, color, marker) in zip((200, 300, 600), BLUE_SEQ, ("o", "s", "^")):
        part = grid[grid["n_body"] == n_body].sort_values("seg")
        _curve(
            ax_c, part["seg"], percent(part["V_rel_to_ref"]),
            color, marker, "-", rf"$N_{{\mathrm{{body}}}}={n_body}$",
        )
    reference_line(ax_c, 0.0)
    ax_c.set_xlabel("segments per µm")
    ax_c.set_ylabel(r"$\Delta V_{\mathrm{sink}}\,/\,V_{\mathrm{ref}}$ (%)")
    ax_c.set_xticks([5, 10, 20])
    ax_c.legend(loc="lower right", borderaxespad=0.15, title=None)

    drawn = tol[tol["V_rel_to_ref"] != 0].sort_values("sweep_value")
    ax_d.plot(
        drawn["sweep_value"], np.abs(percent(drawn["V_rel_to_ref"])),
        color=GREEN, marker="^", ms=3.2, lw=0.9, label="1 filament",
    )
    ax_d.set_xscale("log")
    ax_d.set_yscale("log")
    ax_d.invert_xaxis()
    ax_d.set_xlabel("AbsTol = RelTol (log)")
    ax_d.set_ylabel(r"$|\Delta V_{\mathrm{sink}}|\,/\,V_{\mathrm{ref}}$ (%)")
    ax_d.legend(loc="upper right", borderaxespad=0.15)

    _curve(ax_e, eps["sweep_value"], percent(eps["V_rel_to_ref"]),
           BLUE, "o", "-", "sphere")
    reference_line(ax_e, 0.0)
    ax_e.set_xlabel(r"$\varepsilon/h$")
    ax_e.set_ylabel(r"$\Delta V_{\mathrm{sink}}\,/\,V_{\mathrm{ref}}$ (%)")
    ax_e.legend(loc="upper right", borderaxespad=0.15)

    _curve(ax_f, beta_one["sweep_value"], percent(beta_one["V_rel_to_ref"]),
           GREEN, "^", "-", "1 filament")
    _curve(ax_f, beta_four["sweep_value"], percent(beta_four["V_rel_to_ref"]),
           PINK, "D", "--", "4 filaments")
    reference_line(ax_f, 0.0)
    ax_f.set_xlabel(r"$\beta$")
    ax_f.set_ylabel(r"$\Delta V_{\mathrm{sink}}\,/\,V_{\mathrm{ref}}$ (%)")
    ax_f.legend(loc="upper right", borderaxespad=0.15)

    for ax in axes.ravel():
        polish(ax)
    panel_letters(axes, "abcdef")
    return fig


def _analytic_sink(alpha_deg: np.ndarray, v0: float, v90: float) -> np.ndarray:
    rad = np.deg2rad(alpha_deg)
    return v0 * np.cos(rad) ** 2 + v90 * np.sin(rad) ** 2


def _analytic_drift(alpha_deg: np.ndarray, v0: float, v90: float) -> np.ndarray:
    rad = np.deg2rad(alpha_deg)
    return np.abs(v90 - v0) * np.sin(rad) * np.cos(rad)


def fig_capsule() -> plt.Figure:
    table = load_csv("capsule_orientation.csv").sort_values("alpha_deg")
    v0 = float(table.loc[table["alpha_deg"] == 0, "V_sink_um_s"].iloc[0])
    v90 = float(table.loc[table["alpha_deg"] == 90, "V_sink_um_s"].iloc[0])
    dense = np.linspace(0, 90, 361)

    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(W120, 58 * MM), layout="constrained")
    axes[0].plot(dense, _analytic_sink(dense, v0, v90), color=ORANGE, lw=0.9, label=r"$V_0\cos^2\alpha+V_{90}\sin^2\alpha$")
    axes[0].plot(table["alpha_deg"], table["V_sink_um_s"], linestyle="none",
                 marker="o", ms=3.4, color=BLUE, label="simulation")
    axes[0].set_xlabel(r"tilt $\alpha$ (°)")
    axes[0].set_ylabel(r"$V_{\mathrm{sink}}$ (µm s$^{-1}$)")
    axes[0].set_xlim(-2, 92)
    axes[0].legend(loc="lower right", borderaxespad=0.15, handlelength=2.2)

    axes[1].plot(
        dense, _analytic_drift(dense, v0, v90) * 1e3,
        color=ORANGE, lw=0.9, label=r"$|V_{90}-V_0|\sin\alpha\cos\alpha$",
    )
    axes[1].plot(
        table["alpha_deg"], table["drift_x_um_s"] * 1e3,
        linestyle="none", marker="o", ms=3.4, color=BLUE, label="simulation",
    )
    axes[1].set_xlabel(r"tilt $\alpha$ (°)")
    axes[1].set_ylabel(r"horizontal drift ($10^{-3}$ µm s$^{-1}$)")
    axes[1].set_xlim(-2, 92)
    axes[1].set_ylim(-0.6, 16.2)
    axes[1].legend(loc="upper left", borderaxespad=0.15, handlelength=2.2)
    for ax in axes:
        polish(ax)
    panel_letters(axes, "ab")
    return fig


def _omega_analytic(beta: np.ndarray) -> np.ndarray:
    """Torque balance.  The drag that reproduces the saved predictions is 4 π² μ (β r)²."""
    mu = 1.0e-3
    radius = 3.0e-9
    length = 3.0e-6
    torque = 1.0e-21
    omega_max = 300.0
    drag_length = 4.0 * np.pi ** 2 * mu * (beta * radius) ** 2 * length
    return torque / (drag_length + torque / omega_max)


def fig_motor() -> plt.Figure:
    table = load_csv("validation_motor.csv")
    segments = table[table["id"].str.startswith("VAL-motor-S")].sort_values("seg_per_um")
    betas = table[table["id"].str.startswith("VAL-motor-beta") | table["id"].eq("VAL-motor-S10")]
    betas = betas.sort_values("beta")
    beta_line = np.linspace(0.92, 1.78, 200)
    predicted = _omega_analytic(betas["beta"].to_numpy())
    if np.max(np.abs(predicted - betas["omega_pred_Hz"].to_numpy()) / betas["omega_pred_Hz"]) > 1e-4:
        raise RuntimeError("Analytic rotor speed does not match validation_motor.csv")
    flat = float(segments["omega_pred_Hz"].iloc[0])

    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(W120, 58 * MM), layout="constrained")
    axes[0].axhline(flat, color=ORANGE, lw=0.9, label="analytic")
    axes[0].plot(
        segments["seg_per_um"], segments["omega_rel_Hz"],
        linestyle="none", marker="o", ms=3.4, color=BLUE, label="simulation",
    )
    axes[0].set_xscale("log")
    axes[0].set_xlim(4, 50)
    axes[0].set_ylim(flat * 0.90, flat * 1.10)
    axes[0].set_xlabel("segments per µm (log)")
    axes[0].set_ylabel(r"rotor speed $\omega_{\mathrm{rel}}$ (Hz)")
    axes[0].legend(loc="lower right", borderaxespad=0.2)

    axes[1].plot(beta_line, _omega_analytic(beta_line), color=ORANGE, lw=0.9, label="analytic")
    axes[1].plot(
        betas["beta"], betas["omega_rel_Hz"],
        linestyle="none", marker="o", ms=3.4, color=BLUE, label="simulation (10 µm$^{-1}$)",
    )
    axes[1].set_xlim(0.92, 1.78)
    axes[1].set_xlabel(r"hydrodynamic radius factor $\beta$")
    axes[1].set_ylabel(r"rotor speed $\omega_{\mathrm{rel}}$ (Hz)")
    axes[1].legend(loc="upper right", borderaxespad=0.15)
    for ax in axes:
        polish(ax)
    panel_letters(axes, "ab")
    return fig


def _fdt_display_rows(part: pd.DataFrame) -> pd.DataFrame:
    """Keep every lag of a short series.  Thin a 0.05 s grid to the lags that show the relaxation."""
    if len(part) <= 16:
        return part
    targets = np.array([0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0])
    tau = part["tau_s"].to_numpy()
    chosen = []
    for target in targets:
        j = int(np.argmin(np.abs(tau - target)))
        if abs(tau[j] - target) <= max(0.02, 0.25 * target):
            chosen.append(j)
    if not chosen:
        return part
    return part.iloc[sorted(set(chosen))]


def fig_fdt() -> plt.Figure:
    """Einstein ratio of the body.

    Sphere and bare-capsule series are the body-only ensembles.  When
    ``VAL-ein-th-pili1`` and ``VAL-ein-th-pili4`` are in the table, the
    filament series are those runs: bending orders 1 and 2 thermalised,
    increments taken after the modes have equilibrated, and ``D_t``, ``D_r``
    equal to ``kT`` times the rigid-shape mobility.  Until those rows exist
    the filament series stay the body-only diagnostic.
    """
    table = load_csv("validation_fdt.csv")
    have_thermal = {"VAL-ein-th-pili1", "VAL-ein-th-pili4"}.issubset(set(table["id"]))
    if have_thermal:
        filaments = [
            ("VAL-ein-th-pili1", "capsule + 1 filament", GREEN, "^", "-.", 0.0, False),
            ("VAL-ein-th-pili4", "capsule + 4 filaments", PINK, "D", ":", 0.0, False),
        ]
    else:
        filaments = [
            ("VAL-ein-pili1", "capsule + 1 filament", GREEN, "^", "-.", 0.002, False),
            ("VAL-ein-pili4", "capsule + 4 filaments", PINK, "D", ":", 0.006, False),
        ]
    series = [
        ("VAL-ein-sphere", "sphere", BLUE, "o", "-", -0.006, False),
        ("VAL-ein-sphere-fine", r"sphere, $\Delta t_{\mathrm{B}}=0.01$ s", BLUE, "o", "--", 0.0, True),
        ("VAL-ein-nopili", "capsule", ORANGE, "s", "-", -0.002, False),
        *filaments,
    ]
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(W120, 72 * MM), layout="constrained")
    handles = []
    drawn = [[], []]
    tau_max = 0.0
    for key, label, color, marker, linestyle, shift, open_marker in series:
        part = table[table["id"] == key].sort_values("tau_s")
        part = _fdt_display_rows(part)
        trans = part
        rot = part[part["R_R_valid_window"] == 1]
        x_t = trans["tau_s"].to_numpy() + shift
        x_r = rot["tau_s"].to_numpy() + shift
        y_t = trans["R_T"].to_numpy()
        e_t = trans["R_T_se"].to_numpy()
        y_r = rot["R_R"].to_numpy()
        e_r = rot["R_R_se"].to_numpy()
        face = "none" if open_marker else color
        axes[0].errorbar(
            x_t, y_t, yerr=e_t,
            color=color, marker=marker, mfc=face, mec=color, ms=3.2,
            lw=0.8, ls=linestyle, elinewidth=0.5, capsize=1.4, capthick=0.5,
        )
        axes[1].errorbar(
            x_r, y_r, yerr=e_r,
            color=color, marker=marker, mfc=face, mec=color, ms=3.2,
            lw=0.8, ls=linestyle, elinewidth=0.5, capsize=1.4, capthick=0.5,
        )
        drawn[0].append(np.column_stack([y_t - e_t, y_t + e_t]))
        if y_r.size:
            drawn[1].append(np.column_stack([y_r - e_r, y_r + e_r]))
        tau_max = max(tau_max, float(np.max(trans["tau_s"])) if len(trans) else 0.0)
        handles.append(Line2D(
            [0], [0], color=color, marker=marker, mfc=face, mec=color,
            ms=3.2, lw=0.8, ls=linestyle, label=label,
        ))
    for ax, ylabel, samples in zip(axes, (
        r"Var$|\Delta X|\,/\,(6\,D_{\mathrm{t}}\,\tau)$",
        r"Var$|\Delta\theta|\,/\,(6\,D_{\mathrm{r}}\,\tau)$",
    ), drawn):
        reference_line(ax, 1.0)
        ax.set_xlabel(r"$\tau$ (s)")
        ax.set_ylabel(ylabel)
        polish(ax)
        if have_thermal and samples:
            span = np.concatenate(samples, axis=0)
            lo = min(float(np.min(span)), 1.0)
            hi = max(float(np.max(span)), 1.0)
            pad = 0.08 * (hi - lo)
            ax.set_ylim(lo - pad, hi + pad)
    if have_thermal:
        axes[0].set_xlim(-0.02, tau_max * 1.04)
        axes[1].set_xlim(-0.02, tau_max * 1.04)
    else:
        axes[0].set_xlim(-0.01, 0.33)
        axes[1].set_xlim(-0.01, 0.33)
        axes[0].set_ylim(0.52, 1.18)
        axes[1].set_ylim(0.0, 1.28)
    fig.legend(
        handles=handles, loc="outside upper center", ncol=3,
        borderaxespad=0.0, handlelength=2.6, columnspacing=1.0,
    )
    panel_letters(axes, "ab")
    return fig


# ---------------------------------------------------------------------------
# Inverse fit and thermal acceptance
# ---------------------------------------------------------------------------

def _q2_phase(table_coef: np.ndarray, fit_coef: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Phase of q2 relative to 2 arg(q1), on one common 2π branch.

    The absolute argument of q2 falls along the flagellum because it tracks
    twice the first-harmonic phase.  The published panel is this relative
    phase, which rises from about 2.5 rad to about 6 rad.
    """
    def relative(coef: np.ndarray) -> np.ndarray:
        first = _complex_coeff(coef, 1)
        second = _complex_coeff(coef, 3)
        return np.unwrap(np.angle(second) - 2.0 * np.angle(first))

    table_phase = relative(table_coef)
    fit_phase = relative(fit_coef)
    turns = np.round((np.median(table_phase) - np.median(fit_phase)) / (2.0 * np.pi))
    return table_phase, fit_phase + 2.0 * np.pi * turns


def _harmonic(coef: np.ndarray, phase: float) -> np.ndarray:
    values = np.array(coef[:, 0], dtype=float, copy=True)
    column = 1
    order = 1
    while column + 1 < coef.shape[1]:
        values += coef[:, column] * np.cos(order * phase) + coef[:, column + 1] * np.sin(order * phase)
        column += 2
        order += 1
    return values


def _complex_coeff(coef: np.ndarray, real_col: int) -> np.ndarray:
    # MATLAB stores a = c(:, real) - 1i * c(:, imag), with imag in the next column.
    return coef[:, real_col] - 1j * coef[:, real_col + 1]


def fig_inverse() -> plt.Figure:
    loaded = sio.loadmat(
        ROOT / "output" / "chlamy_inverse" / "fit_1p00A.mat",
        squeeze_me=True, struct_as_record=False,
    )
    tab = sio.loadmat(
        ROOT / "output" / "chlamy_inverse" / "tab_1p00A.mat",
        squeeze_me=True, struct_as_record=False,
    )["tab"]
    fit_coef = np.asarray(loaded["fit"].coef[0], dtype=float)
    table_coef = np.asarray(tab.coef[0], dtype=float)
    samples = np.asarray(tab.q_samples[0], dtype=float)
    nodes = table_coef.shape[0]
    arc = np.arange(nodes) / nodes

    # Sample 10 is the saved frame whose q spans about -6 to +4, matching the
    # previous panel of one instant.
    instant = 10
    phase = float(tab.w_star) * float(tab.t_s[instant]) / float(tab.T_0)
    table_q = samples[:, instant]
    fit_q = _harmonic(fit_coef, phase)

    def amplitude(coef: np.ndarray, real_col: int) -> np.ndarray:
        return np.abs(_complex_coeff(coef, real_col))

    def phase_of(coef: np.ndarray, real_col: int) -> np.ndarray:
        return np.unwrap(np.angle(_complex_coeff(coef, real_col)))

    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(W183, 112 * MM), layout="constrained", sharex=True)
    pairs = (
        (axes[0, 0], amplitude(table_coef, 1), amplitude(fit_coef, 1), r"$|q_1|$ (E/L)"),
        (axes[0, 1], phase_of(table_coef, 1), phase_of(fit_coef, 1), r"phase of $q_1$ (rad)"),
        (axes[0, 2], table_coef[:, 0], fit_coef[:, 0], r"static moment $c_0$ (E/L)"),
        (axes[1, 0], amplitude(table_coef, 3), amplitude(fit_coef, 3), r"$|q_2|$ (E/L)"),
        (axes[1, 1], *_q2_phase(table_coef, fit_coef), r"phase of $q_2$ (rad)"),
        (axes[1, 2], table_q, fit_q, r"$q(s)$ at one instant (E/L)"),
    )
    for ax, table_y, fit_y, ylabel in pairs:
        ax.plot(arc, table_y, color=BLUE, lw=0.8, label="inverse table")
        ax.plot(arc, fit_y, color=ORANGE, lw=0.9, ls="--", label="analytic fit")
        ax.set_ylabel(ylabel)
        polish(ax)
    for ax in axes[1]:
        ax.set_xlabel(r"arc length $s/L$")
    axes[0, 0].set_xlim(0, 1)
    fig.legend(
        handles=[
            Line2D([0], [0], color=BLUE, lw=0.8, label="inverse table"),
            Line2D([0], [0], color=ORANGE, lw=0.9, ls="--", label="analytic fit"),
        ],
        loc="outside upper center", ncol=2, borderaxespad=0.0,
    )
    panel_letters(axes, "abcdef")
    return fig


def _mean_pair(block: np.ndarray) -> np.ndarray:
    return 0.5 * (block[4] + block[5])


def _se_pair(block: np.ndarray) -> np.ndarray:
    return 0.5 * np.sqrt(block[4] ** 2 + block[5] ** 2)


def _case(record, name: str):
    return getattr(record, name)


def fig_thermal() -> plt.Figure:
    """Pilus thermal acceptance.  Colours match across the six panels.

    Blue is pilus noise and vermilion is body-only noise.  In a–c a solid
    line is bending order 1 and a dashed line is order 2.  In d–f the dashed
    curve is the mobility prediction at that lag.  The grey line is the
    model reference stored with the case: the rigid-body value for D⊥ and
    the flexible-body value for D∥.
    """
    record = sio.loadmat(
        ROOT / "output" / "test_pb" / "pba_analysis.mat",
        squeeze_me=True, struct_as_record=False,
    )["R"]
    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(W183, 142 * MM), layout="constrained")

    time_panels = (
        (axes[0, 0], "D1", "D1c", "1 filament, 293.15 K"),
        (axes[0, 1], "C1", None, "1 filament, 2.93 K"),
        (axes[0, 2], "D4", "D4c", "4 filaments, 293.15 K"),
    )
    styles = (
        (0, BLUE, "-", "order 1, pilus noise"),
        (1, ORANGE, "-", "order 2, pilus noise"),
    )
    body_styles = (
        (0, BLUE, (0, (0.7, 0.8)), "order 1, body-only noise"),
        (1, ORANGE, (0, (0.7, 0.8)), "order 2, body-only noise"),
    )
    for ax, pilus_name, body_name, note in time_panels:
        pilus = _case(record, pilus_name)
        series = [(pilus, styles)]
        if body_name is not None:
            series.append((_case(record, body_name), body_styles))
        drawn = []
        for case, spec in series:
            time = np.asarray(case.grow.t, dtype=float)
            sim = np.asarray(case.grow.sim, dtype=float)
            for row, color, linestyle, _label in spec:
                line = sim[row]
                ax.plot(time, line, color=color, lw=0.55, ls=linestyle)
                drawn.append(line)
            if body_name is None or case is not pilus:
                # Linear theory of the noise model drawn in this panel.
                # Body-only panels use the body theory; the cold panel has
                # only the coupled theory, which sits on the equipartition line.
                theory = np.asarray(case.grow.theory, dtype=float)
                for row in (0, 1):
                    ax.plot(time, theory[row], color=BLACK, lw=0.8, ls="--")
                    drawn.append(theory[row])
        reference_line(ax, 1.0)
        stacked = np.concatenate([np.ravel(item) for item in drawn])
        ax.set_ylim(0.0, float(np.nanmax(stacked)) / 0.62)
        ax.set_xlim(time[0], time[-1])
        ax.set_xlabel(r"$t$ (s)")
        ax.set_ylabel(r"$a^{2}\,/\,k_{\mathrm{B}}T$")
        ax.text(
            0.98, 0.98, note, transform=ax.transAxes,
            ha="right", va="top", fontsize=7, fontfamily="Arial",
        )
        polish(ax)

    lag_panels = (
        (axes[1, 0], "D1", "D1c", "1 filament", "perp"),
        (axes[1, 1], "D4", "D4c", "4 filaments", "perp"),
        (axes[1, 2], "D4", "D4c", "4 filaments", "para"),
    )
    for ax, pilus_name, body_name, note, kind in lag_panels:
        lows = []
        highs = []
        for case, color, marker in (
            (_case(record, pilus_name), BLUE, "o"),
            (_case(record, body_name), ORANGE, "s"),
        ):
            lags = np.asarray(case.lags_s, dtype=float)
            sim = np.asarray(case.D_sim, dtype=float)
            se = np.asarray(case.D_se, dtype=float)
            theory = np.asarray(case.D_th, dtype=float)
            if kind == "perp":
                y = _mean_pair(sim)
                err = _se_pair(se)
                y_th = _mean_pair(theory)
                ref = float(_mean_pair(np.asarray(case.D_rigid, dtype=float)))
            else:
                y = sim[3]
                err = se[3]
                y_th = theory[3]
                ref = float(np.asarray(case.D_flex, dtype=float)[3])
            ax.errorbar(
                lags, y, yerr=err, color=color, marker=marker, ms=3.3,
                lw=0, elinewidth=0.6, capsize=1.6, capthick=0.5, zorder=3,
            )
            ax.plot(lags, y_th, color=color, lw=0.8, ls="--", zorder=2)
            ax.axhline(ref, color=REF, lw=0.5, zorder=0)
            lows.append(float(np.min(np.r_[y - err, y_th, ref])))
            highs.append(float(np.max(np.r_[y + err, y_th, ref])))
        lo, hi = min(lows), max(highs)
        ax.set_xscale("log")
        polish(ax)
        if kind == "para":
            ax.set_yscale("log")
            ax.set_ylim(lo / 1.35, hi * 1.35)
            ax.text(
                0.03, 0.06, note, transform=ax.transAxes,
                ha="left", va="bottom", fontsize=7, fontfamily="Arial",
            )
        else:
            pad = 0.16 * (hi - lo)
            ax.set_ylim(lo - pad, hi + pad)
            ax.text(
                0.97, 0.96, note, transform=ax.transAxes,
                ha="right", va="top", fontsize=7, fontfamily="Arial",
            )
        ax.set_xlabel(r"lag $\tau$ (s, log)")
    axes[1, 0].set_ylabel(r"$D_{\mathrm{perp}}$ (rad$^{2}$ s$^{-1}$)")
    axes[1, 1].set_ylabel(r"$D_{\mathrm{perp}}$ (rad$^{2}$ s$^{-1}$)")
    axes[1, 2].set_ylabel(r"$D_{\mathrm{para}}$ (rad$^{2}$ s$^{-1}$)")
    handles = [
        Line2D([0], [0], color=BLUE, lw=0.9, label="order 1, pilus noise"),
        Line2D([0], [0], color=ORANGE, lw=0.9, label="order 2, pilus noise"),
        Line2D([0], [0], color=BLUE, lw=0.9, ls=(0, (0.7, 0.8)), label="order 1, body-only noise"),
        Line2D([0], [0], color=ORANGE, lw=0.9, ls=(0, (0.7, 0.8)), label="order 2, body-only noise"),
        Line2D([0], [0], color=BLACK, lw=0.8, ls="--", label="linear theory"),
    ]
    fig.legend(
        handles=handles, loc="outside upper center", ncol=3,
        borderaxespad=0.0, handlelength=2.8, columnspacing=0.9,
    )
    lag_handles = [
        Line2D([0], [0], color=BLUE, lw=0, marker="o", ms=3.3, label="pilus thermalisation"),
        Line2D([0], [0], color=ORANGE, lw=0, marker="s", ms=3.3, label="body-only noise"),
        Line2D([0], [0], color=BLUE, lw=0.8, ls="--", label="mobility prediction"),
        Line2D([0], [0], color=REF, lw=0.8, label="reference level"),
    ]
    fig.legend(
        handles=lag_handles, loc="outside lower center", ncol=4,
        borderaxespad=0.0, handlelength=2.2, columnspacing=0.9,
    )
    panel_letters(axes, "abcdef")
    return fig


def main() -> None:
    apply_style()
    jobs = [
        ("fig_stokes_sphere", FIG / "fig_stokes_sphere", fig_stokes),
        ("fig_arrangement", FIG / "fig_arrangement", fig_arrangement),
        ("fig_brownian", FIG / "fig_brownian", fig_brownian),
        ("fig_convergence", FIG / "fig_convergence", fig_convergence),
        ("fig_capsule_orientation", FIG / "fig_capsule_orientation", fig_capsule),
        ("fig_motor", FIG / "fig_motor", fig_motor),
        ("fig_fdt", FIG / "fig_fdt", fig_fdt),
        ("fit_1p00A", ROOT / "output" / "chlamy_inverse" / "fit_1p00A", fig_inverse),
        ("fig_pili_thermal_acceptance", ROOT / "output" / "test_pb" / "fig_pili_thermal_acceptance", fig_thermal),
    ]
    qa_dir = os.environ.get("SCIPILOT_SCRIPTS", "")
    if qa_dir and Path(qa_dir).is_dir():
        sys.path.insert(0, qa_dir)
    try:
        from visual_qa import audit_layout, print_report
    except ImportError:
        audit_layout = print_report = None

    for name, path, builder in jobs:
        fig = builder()
        save_figure(fig, path)
        preview(fig, name)
        print(f"\n== {name}  size_in={fig.get_size_inches().tolist()}")
        if audit_layout is not None and print_report is not None:
            print_report(audit_layout(fig))
        plt.close(fig)


if __name__ == "__main__":
    main()
