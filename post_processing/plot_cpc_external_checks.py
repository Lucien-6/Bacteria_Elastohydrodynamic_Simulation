"""CPC validation figures for the external helix and cantilever checks.

The pitch panels compare the RSS propulsive coefficients with Johnson
slender-body theory and with Gray-Hancock resistive-force theory.  The
cantilever panels compare the body-fixed bending time with the local
cantilever formula.  Print size follows Elsevier artwork widths.

Author: Lucien (lucien-6@qq.com)
Created: 2026-11-15
"""

from __future__ import annotations

import os
import string
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "output" / "cpc_checks"
FIG = DATA / "figures"

MM = 1.0 / 25.4
W90 = 90.0 * MM
W140 = 140.0 * MM
W190 = 190.0 * MM

BLUE = "#0072B2"
ORANGE = "#E69F00"
VERMILION = "#D55E00"
GREY = "#7A7A7A"


def apply_style() -> None:
    plt.rcParams.update({
        "figure.figsize": (3.5, 2.625),
        "figure.dpi": 150,
        "figure.constrained_layout.use": True,
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "lines.linewidth": 1.0,
        "lines.markersize": 4,
        "axes.linewidth": 0.6,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
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
    })


def add_panel_labels(
    fig: plt.Figure,
    axes=None,
    fontsize: float = 8,
    x_offset_pt: float = -18.0,
    y_offset_pt: float = 2.0,
) -> None:
    """Place aligned (a), (b), ... labels anchored to each panel's top-left (Elsevier CPC)."""
    if axes is None:
        axes = [ax for ax in fig.axes if ax.get_subplotspec() is not None]
    for ax, letter in zip(axes, string.ascii_lowercase):
        ax.annotate(
            f"({letter})",
            xy=(0.0, 1.0),
            xycoords="axes fraction",
            xytext=(x_offset_pt, y_offset_pt),
            textcoords="offset points",
            fontsize=fontsize,
            fontweight="bold",
            ha="right",
            va="bottom",
            color="black",
            annotation_clip=False,
        )


def polish(ax: plt.Axes) -> None:
    ax.tick_params(direction="in", top=True, right=True, length=2.4, width=0.5)
    ax.margins(x=0.06, y=0.08)


def production_pitch(table: pd.DataFrame) -> pd.DataFrame:
    rows = table[
        (table["is_production"] == 1)
        & (table["with_body"] == 0)
        & np.isclose(table["beta"], np.sqrt(2))
        & table["name"].isin(["pitch_short", "pitch_mid", "pitch_long"])
    ].sort_values("lambda_over_R")
    return rows


def production_all(table: pd.DataFrame) -> pd.DataFrame:
    rows = table[
        (table["is_production"] == 1)
        & (table["with_body"] == 0)
        & np.isclose(table["beta"], np.sqrt(2))
    ].copy()
    order = ["pitch_short", "pitch_mid", "pitch_long", "length_long"]
    rows["order"] = rows["name"].map({name: i for i, name in enumerate(order)})
    return rows.sort_values("order")


def fig_helix(table: pd.DataFrame) -> plt.Figure:
    pitch = production_pitch(table)
    x = pitch["lambda_over_R"].to_numpy()
    series = (
        ("F", r"$|F|\,/\,(\mu\omega R^{2})$"),
        ("T", r"$T\,/\,(\mu\omega R^{3})$"),
        ("D", r"$D\,/\,(\mu U R)$"),
    )
    apply_style()
    fig = plt.figure(figsize=(W190, 118 * MM), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=(1.15, 1.0))
    axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    ax_err = fig.add_subplot(grid[1, 0:2])
    ax_beta = fig.add_subplot(grid[1, 2])

    for ax, (key, ylabel) in zip(axes, series):
        y_sbt = np.abs(pitch[f"{key}_johnson"])
        y_rft = np.abs(pitch[f"{key}_rft"])
        y_rss = np.abs(pitch[f"{key}_rss"])
        ax.plot(x, y_rft, color=VERMILION, lw=0.8, ls=":", marker="^", ms=3.6, zorder=2)
        ax.plot(
            x, y_sbt, color=ORANGE, lw=0.9, ls="--", marker="s", ms=4.6,
            mfc="none", mew=0.8, zorder=3,
        )
        ax.plot(x, y_rss, color=BLUE, lw=0.9, marker="o", ms=3.4, zorder=4)
        ax.set_xlabel(r"pitch $\lambda/R$")
        ax.set_ylabel(ylabel)
        ax.set_xticks([3, 6.5, 12])
        ax.set_xticklabels(["3", "6.5", "12"])
        polish(ax)

    cases = production_all(table)
    xpos = np.arange(len(cases))
    err_handles = []
    for key, marker, label in (
        ("F", "o", "thrust"),
        ("T", "s", "torque"),
        ("D", "^", "drag"),
    ):
        rel = (cases[f"{key}_rss"] - cases[f"{key}_johnson"]) / np.abs(cases[f"{key}_johnson"])
        ax_err.plot(xpos, rel, color=BLUE, marker=marker, ms=4.0, lw=0.8, label=label)
        err_handles.append(Line2D(
            [0], [0], color=BLUE, marker=marker, ms=4.0, lw=0.8, label=label,
        ))
    ax_err.axhline(0.0, color=GREY, lw=0.5, zorder=0)
    ax_err.axhline(0.05, color=GREY, lw=0.5, ls="--", zorder=0)
    ax_err.axhline(-0.05, color=GREY, lw=0.5, ls="--", zorder=0)
    ax_err.set_xticks(xpos)
    ax_err.set_xticklabels(["3\n20", "6.5\n20", "12\n20", "6.5\n40"])
    ax_err.set_xlabel(r"$\lambda/R$ ($L/R$ below each tick)")
    ax_err.set_ylabel(r"$(c_{\mathrm{RSS}}-c_{\mathrm{SBT}})\,/\,|c_{\mathrm{SBT}}|$")
    ax_err.set_ylim(-0.055, 0.105)
    ax_err.legend(
        handles=err_handles, loc="upper center", ncol=3, frameon=False,
        borderaxespad=0.1, handlelength=1.6, columnspacing=0.8,
    )
    polish(ax_err)

    beta_rows = table[
        (table["name"] == "pitch_mid")
        & (table["with_body"] == 0)
        & (table["segs_per_turn"] == 32)
    ].sort_values("beta")
    ax_beta.axhline(abs(float(beta_rows["T_johnson"].iloc[0])), color=ORANGE, lw=0.9, ls="--")
    ax_beta.plot(
        beta_rows["beta"], np.abs(beta_rows["T_rss"]),
        color=BLUE, marker="o", ms=4.0, lw=0.8,
    )
    ax_beta.set_xlabel(r"radius factor $\beta$")
    ax_beta.set_ylabel(r"$T\,/\,(\mu\omega R^{3})$")
    ax_beta.set_xticks([1.0, 1.41, 1.65])
    ax_beta.set_xticklabels(["1", r"$\sqrt{2}$", r"$\sqrt{e}$"])
    polish(ax_beta)

    handles = [
        Line2D([0], [0], color=BLUE, marker="o", lw=0.9, ms=3.4, label="RSS"),
        Line2D(
            [0], [0], color=ORANGE, marker="s", mfc="none", mew=0.8,
            lw=0.9, ls="--", ms=4.6, label="Johnson SBT",
        ),
        Line2D(
            [0], [0], color=VERMILION, marker="^", lw=0.8, ls=":",
            ms=3.6, label="Gray–Hancock RFT",
        ),
    ]
    fig.legend(
        handles=handles, loc="outside upper center", ncol=3,
        borderaxespad=0.0, handlelength=2.4,
    )
    add_panel_labels(fig, axes=axes + [ax_err, ax_beta], fontsize=8)
    return fig


def fig_cantilever(table: pd.DataFrame) -> plt.Figure:
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(W140, 62 * MM), layout="constrained")
    fine = table[table["seg_per_um"] == 20]
    base = table[table["seg_per_um"] == 10].sort_values("length_um")

    for column, color, marker, label in (
        ("tau_ref_s", ORANGE, "s", "cantilever formula"),
        ("tau_fixed_s", BLUE, "o", "body fixed"),
        ("tau_free_s", VERMILION, "^", "body free"),
    ):
        axes[0].plot(
            base["length_um"], base[column], color=color, marker=marker,
            ms=4.0, lw=0.9, label=label,
        )
    axes[0].plot(
        fine["length_um"], fine["tau_fixed_s"], color=BLUE, marker="o",
        ms=4.2, mfc="none", lw=0, label="body fixed, 20 µm$^{-1}$",
    )
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"filament length $L$ (µm)")
    axes[0].set_ylabel(r"relaxation time $\tau$ (s)")
    axes[0].set_xticks([3, 6])
    axes[0].legend(loc="upper left", borderaxespad=0.15, handlelength=1.8)
    polish(axes[0])

    labels = ["3 µm\n10 µm$^{-1}$", "3 µm\n20 µm$^{-1}$", "6 µm\n10 µm$^{-1}$"]
    xpos = np.arange(len(table))
    axes[1].axhspan(0.90, 1.10, color="#E6E6E6", zorder=0)
    axes[1].axhline(1.0, color=GREY, lw=0.6, zorder=1)
    axes[1].plot(xpos, table["ratio"], color=BLUE, marker="o", ms=4.2, lw=0)
    axes[1].set_xticks(xpos)
    axes[1].set_xticklabels(labels)
    axes[1].set_ylabel(r"$\tau_{\mathrm{fixed}}\,/\,\tau_{\mathrm{ref}}$")
    axes[1].set_ylim(0.82, 1.18)
    polish(axes[1])

    add_panel_labels(fig, fontsize=8)
    return fig


def fig_performance(ops: pd.DataFrame, e2e: pd.DataFrame) -> plt.Figure:
    """Two-panel Elsevier CPC performance benchmark figure (140 mm)."""
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(W140, 62 * MM), layout="constrained")

    # (a) Per-call RHS & Jacobian speedup across (N_body, M_pili)
    ax0 = axes[0]
    xpos = np.arange(len(ops))
    w = 0.36
    short_labels = [
        f"M={int(r.M_pili)}\nN={int(r.N_body_actual)}" for _, r in ops.iterrows()
    ]
    b1 = ax0.bar(
        xpos - w / 2, ops["rhs_speedup"], width=w,
        color=BLUE, label="RHS (Block Schur)",
    )
    b2 = ax0.bar(
        xpos + w / 2, ops["jac_speedup"], width=w,
        color=ORANGE, label="Jacobian (3-colour)",
    )
    ax0.axhline(1.0, color=GREY, lw=0.6, ls="--", zorder=0)
    ax0.set_yscale("log")
    ax0.set_xticks(xpos)
    ax0.set_xticklabels(short_labels, fontsize=6.5)
    ax0.set_ylabel("Speedup vs legacy dense + FD")
    ax0.set_ylim(0.8, max(float(ops["jac_speedup"].max()) * 2.5, 100.0))
    ax0.legend(loc="upper left", borderaxespad=0.15, handlelength=1.4)
    for bar in list(b1) + list(b2):
        h = bar.get_height()
        ax0.annotate(
            f"{h:.0f}x" if h >= 10 else f"{h:.1f}x",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 2),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=5.8,
        )
    polish(ax0)

    # (b) End-to-end integration wall time
    ax1 = axes[1]
    xe = np.arange(len(e2e))
    w_e = 0.32
    off_e = 0.20
    e_labels = ["Det. settling\n(1 pilus, 0.5 s)", "Brownian\n(1 pilus, 0.2 s)"]
    if len(e2e) != 2:
        e_labels = [str(x) for x in e2e["case_id"]]
    ax1.bar(xe - off_e, e2e["wall_legacy_s"], width=w_e, color=GREY, label="Legacy (dense + FD)")
    ax1.bar(xe + off_e, e2e["wall_v10_s"], width=w_e, color=BLUE, label="BacEHD v10.2.1")
    ax1.set_xticks(xe)
    ax1.set_xticklabels(e_labels, fontsize=6.5)
    ax1.set_xlim(-0.55, len(e2e) - 0.35)
    ax1.set_ylabel("Wall-clock time (s)")
    ax1.set_ylim(0, float(e2e["wall_legacy_s"].max()) * 1.25)
    ax1.legend(loc="upper left", borderaxespad=0.15, handlelength=1.4)
    for idx, (_, row) in enumerate(e2e.iterrows()):
        ax1.annotate(
            f"{row.speedup:.1f}x\nerr {row.end_pos_rel_diff:.1e}",
            xy=(idx + off_e, row.wall_v10_s),
            xytext=(3, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=5.8,
        )
    polish(ax1)

    add_panel_labels(fig, fontsize=8)
    return fig


def save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".png"), dpi=600)
    fig.savefig(path.with_suffix(".svg"))


def main() -> None:
    helix = pd.read_csv(DATA / "helix_propulsive_matrix.csv")
    cantilever = pd.read_csv(DATA / "cantilever_relaxation.csv")
    jobs = [
        ("fig_helix_propulsion", fig_helix(helix)),
        ("fig_cantilever_relaxation", fig_cantilever(cantilever)),
    ]
    ops_csv = DATA / "performance_benchmark_operators.csv"
    e2e_csv = DATA / "performance_benchmark_e2e.csv"
    if ops_csv.is_file() and e2e_csv.is_file():
        ops = pd.read_csv(ops_csv)
        e2e = pd.read_csv(e2e_csv)
        jobs.append(("fig_performance_benchmark", fig_performance(ops, e2e)))
    qa_dir = os.environ.get("SCIPILOT_SCRIPTS", "")
    if qa_dir and Path(qa_dir).is_dir():
        sys.path.insert(0, qa_dir)
    try:
        from visual_qa import audit_layout, print_report, render_preview
    except ImportError:
        audit_layout = print_report = render_preview = None

    for name, fig in jobs:
        save(fig, FIG / name)
        preview = FIG / f"_{name}_preview.png"
        if render_preview is not None:
            render_preview(fig, str(preview), dpi=150)
        else:
            fig.savefig(preview, dpi=150)
        print(f"\n== {name}")
        if audit_layout is not None and print_report is not None:
            print_report(audit_layout(fig))
        plt.close(fig)


if __name__ == "__main__":
    main()
