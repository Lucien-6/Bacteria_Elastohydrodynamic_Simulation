"""Generate additional Elsevier CPC paper figures for BacEHD.

Produces five publication-ready figures at Elsevier CPC double-column width
(190 mm) using the repository's unified Arial / Okabe-Ito / inward-tick style:
  1. fig_architecture       (Sec. 3.1) — Modular software & data-flow diagram
  2. fig_parallel_and_gpu   (Sec. 5.2-5.3) — parfor strong scaling & GPU flowfield
  3. fig_app_cyanobacteria  (Sec. 6.2) — Pilus arrangement, flexibility & shear
  4. fig_app_ecoli          (Sec. 6.3) — E. coli motor swimming & XZ flow field
  5. fig_app_chlamy         (Sec. 6.4) — Chlamydomonas 3-D beat, spin & superhelix

All numerical plots read directly from verified CSV / MAT outputs; no synthetic
data are generated.

Author: Lucien (lucien-6@qq.com)
Created: 2026-11-25
"""

from __future__ import annotations

import string
from pathlib import Path

import matplotlib._mathtext as _mathtext
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
from matplotlib.collections import PolyCollection
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
CPC_DIR = ROOT / "output" / "cpc_checks"
PAPER_DIR = ROOT / "output" / "paper" / "results"
OUT_DIR = ROOT / "output" / "cpc_paper_figures"

MM = 1.0 / 25.4
W90 = 90.0 * MM
W140 = 140.0 * MM
W190 = 190.0 * MM

# Okabe-Ito colorblind-safe palette
BLUE = "#0072B2"
SKY = "#56B4E9"
ORANGE = "#E69F00"
VERMILION = "#D55E00"
GREEN = "#009E73"
PINK = "#CC79A7"
BLACK = "#1A1A1A"
GREY = "#6E6E6E"
LIGHT_GREY = "#E8E8E8"


def apply_style() -> None:
    _mathtext.SHRINK_FACTOR = 0.82
    plt.rcParams.update({
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.6,
        "lines.linewidth": 1.0,
        "lines.markersize": 4.0,
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
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "savefig.dpi": 600,
        "legend.frameon": False,
    })


def polish(ax: plt.Axes) -> None:
    ax.tick_params(direction="in", top=True, right=True, length=2.4, width=0.5)
    ax.margins(x=0.06, y=0.08)


def add_panel_labels(axes, x_pt: float = -16.0, y_pt: float = 3.0) -> None:
    for ax, letter in zip(axes, string.ascii_lowercase):
        ax.annotate(
            f"({letter})",
            xy=(0.0, 1.0),
            xycoords="axes fraction",
            xytext=(x_pt, y_pt),
            textcoords="offset points",
            fontsize=8,
            fontweight="bold",
            ha="right",
            va="bottom",
            color="black",
            annotation_clip=False,
        )


def save_all(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    base = OUT_DIR / name
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=600)
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(OUT_DIR / f"_{name}_preview.png", dpi=160)
    plt.close(fig)
    print(f"Saved {base.with_suffix('.pdf')}")


# ---------------------------------------------------------------------------
# Figure 1: Modular Software Architecture & Data-Flow Diagram (190 mm)
# ---------------------------------------------------------------------------
def draw_box(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: str,
    edge: str,
    face: str,
    title_bg: str | None = None,
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.008,rounding_size=0.015",
        linewidth=0.85,
        edgecolor=edge,
        facecolor=face,
        zorder=2,
    )
    ax.add_patch(patch)
    if title_bg is not None:
        hdr_h = 0.058
        hdr = FancyBboxPatch(
            (x, y + h - hdr_h),
            w,
            hdr_h,
            boxstyle="round,pad=0.008,rounding_size=0.015",
            linewidth=0.0,
            facecolor=title_bg,
            zorder=2.5,
        )
        ax.add_patch(hdr)
        ax.text(
            x + w / 2,
            y + h - hdr_h / 2 + 0.001,
            title,
            ha="center",
            va="center",
            fontsize=7.0,
            fontweight="bold",
            color="white",
            zorder=4,
        )
        ax.text(
            x + w / 2,
            y + (h - hdr_h) / 2,
            body,
            ha="center",
            va="center",
            fontsize=6.1,
            color=BLACK,
            linespacing=1.28,
            zorder=4,
        )
    else:
        ax.text(
            x + w / 2,
            y + h * 0.78,
            title,
            ha="center",
            va="center",
            fontsize=6.5,
            fontweight="bold",
            color=edge,
            zorder=4,
        )
        ax.text(
            x + w / 2,
            y + h * 0.36,
            body,
            ha="center",
            va="center",
            fontsize=5.4,
            color=BLACK,
            linespacing=1.22,
            zorder=4,
        )


def draw_arrow(
    ax: plt.Axes,
    p0: tuple[float, float],
    p1: tuple[float, float],
    color: str = GREY,
    ls: str = "-",
    rad: float = 0.0,
) -> None:
    arr = FancyArrowPatch(
        p0,
        p1,
        arrowstyle="-|>,head_length=4.5,head_width=2.6",
        connectionstyle=f"arc3,rad={rad}",
        linewidth=0.85,
        linestyle=ls,
        color=color,
        zorder=3,
    )
    ax.add_patch(arr)


def fig_architecture() -> None:
    apply_style()
    fig, ax = plt.subplots(figsize=(W190, 96 * MM))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    # Main 4-stage pipeline across the top/middle band
    # Column 1: Configuration & Validation
    draw_box(
        ax,
        0.015,
        0.46,
        0.205,
        0.48,
        "1. Config & Pre-flight",
        "configs/*.json or inline cfg\n"
        "load_sim_config.m\n"
        "validate_simulation_params.m\n"
        "────────────────────\n"
        "• Parametric or custom mode\n"
        "• Root BC: clamped / motor\n"
        "• Active moment / bg flow\n"
        "• Brownian & pili_thermal",
        edge=BLUE,
        face="#F2F8FC",
        title_bg=BLUE,
    )

    # Column 2: Geometry & Rigid-Body Precomputation
    draw_box(
        ax,
        0.255,
        0.46,
        0.215,
        0.48,
        "2. Geometry & Precomp",
        "init_coupled_params.m\n"
        "capsule_ / custom_body_building\n"
        "canonical_pili_centerline.m\n"
        "────────────────────\n"
        "• S2 Riesz cap mesh cache\n"
        "• Equal-chord PCHIP rods\n"
        "• build_body_precomp.m:\n"
        "  Cholesky(M_BB^loc) once",
        edge=GREEN,
        face="#F1FAF7",
        title_bg=GREEN,
    )

    # Column 3: Core EHD + Brownian Solver Engine
    draw_box(
        ax,
        0.505,
        0.46,
        0.255,
        0.48,
        "3. Coupled EHD & Brownian Core",
        "run_simulation.m  (ode15s BDF/NDF)\n"
        "──────────────────────────\n"
        "• build_coupled_lin_sys.m:\n"
        "  Body-frame Block Schur + weighted\n"
        "  reciprocal RSM / RSS / Rotlet (β=√2)\n"
        "• ode_coupled_rhs & 3-colour Jacobian\n"
        "  + LaggedJacobian + exact motor solve\n"
        "• Body 6-DOF or slow-mode pili thermal\n"
        "  (brown_pili_tau + held-load drift)",
        edge=VERMILION,
        face="#FFF6F2",
        title_bg=VERMILION,
    )

    # Column 4: Diagnostics, Outputs & Post-Processing
    draw_box(
        ax,
        0.795,
        0.46,
        0.190,
        0.48,
        "4. Diagnostics & Output",
        "check_linear_solve / rcond\n"
        "──────────────────\n"
        "• simulation_results.mat\n"
        "• timeseries.csv & log\n"
        "• config.json + RNG seed\n"
        "• Vector PDF / 600 dpi PNG\n"
        "• 4-panel MP4/AVI movie\n"
        "• GPU + Threads flow field",
        edge=ORANGE,
        face="#FFFBF0",
        title_bg="#C98700",
    )

    # Main horizontal arrows
    draw_arrow(ax, (0.228, 0.70), (0.247, 0.70), color=BLACK)
    draw_arrow(ax, (0.478, 0.70), (0.497, 0.70), color=BLACK)
    draw_arrow(ax, (0.768, 0.70), (0.787, 0.70), color=BLACK)

    # Adaptive basis & Brownian restart feedback loop inside/above Stage 3
    ax.text(
        0.632,
        0.972,
        "Event loop: adaptive_basis_check (SO(3) singularity guard)  &  Brownian interval restarts",
        ha="center",
        va="center",
        fontsize=6.1,
        fontstyle="italic",
        color=VERMILION,
    )

    # Bottom row: 3 Pluggable Ecosystem Modules
    draw_box(
        ax,
        0.015,
        0.05,
        0.295,
        0.30,
        "Pluggable Species (chlamydomonas/)",
        "chlamy_inverse_moment -> chlamy_fit_moment\n"
        "(dual-k 1ω+2ω measured-gait active drive)\n"
        "chlamy_make_perp_envelope (zero-twist 3-D OOP)\n"
        "Plugged via active_moment='custom'; SO(3) metrics",
        edge=PINK,
        face="#FCF4F9",
    )

    draw_box(
        ax,
        0.345,
        0.05,
        0.310,
        0.30,
        "Batch & HPC Drivers (batch_sweep/)",
        "batch_repeat_brownian(_HPC_Linux).m\n"
        "batch_run_from_json.m  |  postprocess_cases.m\n"
        "Auto-sized parfor (1 BLAS/worker, Slurm TMPDIR)\n"
        "Reproducible RNG seeds & DataQueue heartbeats",
        edge=BLUE,
        face="#F4F8FB",
    )

    draw_box(
        ax,
        0.690,
        0.05,
        0.295,
        0.30,
        "Verification & Benchmarks (cpc/, validation/)",
        "cpc/run_cpc_tests.m\n"
        "(18 built-in, ~23 s; 2 external, 1-2 min)\n"
        "cpc_helix_propulsive_matrix & cantilever_relaxation\n"
        "cpc_performance_benchmark (Schur + 3-colour Jac)\n"
        "cpc_parallel_and_flowfield_benchmark (parfor/GPU)",
        edge=GREEN,
        face="#F3FAF7",
    )

    # Dashed connector arrows between bottom ecosystem and main pipeline
    draw_arrow(ax, (0.16, 0.36), (0.56, 0.452), color=PINK, ls="--", rad=-0.08)
    draw_arrow(ax, (0.50, 0.36), (0.61, 0.452), color=BLUE, ls="--")
    draw_arrow(ax, (0.82, 0.36), (0.68, 0.452), color=GREEN, ls="--", rad=0.08)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.01)
    save_all(fig, "fig_architecture")


# ---------------------------------------------------------------------------
# Figure 9: Multi-worker parfor Strong Scaling & GPU Flow-Field (190 mm)
# ---------------------------------------------------------------------------
def fig_parallel_and_gpu() -> None:
    par = pd.read_csv(CPC_DIR / "parallel_scaling_benchmark.csv")
    flw = pd.read_csv(CPC_DIR / "flowfield_gpu_benchmark.csv")

    p1t = par[par["blas_threads_per_worker"] == 1].sort_values("n_workers")
    pdef = par[par["mode"] == "serial_default_blas"].iloc[0]

    apply_style()
    fig, axes = plt.subplots(
        1, 3, figsize=(W190, 62 * MM), layout="constrained",
        gridspec_kw={"width_ratios": [1.05, 0.95, 1.15]},
    )

    # (a) parfor strong scaling speedup
    ax = axes[0]
    w = p1t["n_workers"].to_numpy(dtype=float)
    sp = p1t["speedup_vs_1w1t"].to_numpy(dtype=float)
    ax.plot(w, w, color=GREY, ls="--", lw=0.75, label="Ideal linear")
    ax.plot(w, sp, color=BLUE, marker="o", ms=4.2, lw=1.1, label="parfor (1 BLAS/worker)")
    ax.plot(
        [1.0],
        [float(pdef["speedup_vs_1w1t"])],
        color=ORANGE,
        marker="s",
        ms=4.5,
        lw=0,
        label="Serial (24 BLAS threads)",
    )
    for xi, yi in zip(w, sp):
        if xi > 1:
            dx_pt, dy_pt, va_a = (1, -8, "top") if xi == 2 else ((14, -7, "top") if xi == 4 else (0, 5, "bottom"))
            ax.annotate(
                f"{yi:.2f}x",
                xy=(xi, yi),
                xytext=(dx_pt, dy_pt),
                textcoords="offset points",
                ha="center",
                va=va_a,
                fontsize=6.0,
                color=BLUE,
            )
    ax.set_xlabel("Number of parfor workers")
    ax.set_ylabel("Speedup vs 1 worker (1 BLAS)")
    ax.set_xticks([1, 2, 4, 8, 16])
    ax.set_xlim(0.2, 17.2)
    ax.set_ylim(0.0, 17.5)
    ax.legend(loc="upper left", handlelength=1.4, borderaxespad=0.2)
    polish(ax)

    # (b) Parallel efficiency (%)
    ax = axes[1]
    eff = p1t["efficiency_pct"].to_numpy(dtype=float)
    xpos = np.arange(len(w))
    bars = ax.bar(xpos, eff, width=0.56, color=SKY, edgecolor=BLUE, linewidth=0.6)
    ax.axhline(100.0, color=GREY, ls="--", lw=0.6, zorder=0)
    for b, ev in zip(bars, eff):
        ax.annotate(
            f"{ev:.1f}%",
            xy=(b.get_x() + b.get_width() / 2, ev),
            xytext=(0, 2.5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=6.0,
        )
    ax.set_xticks(xpos)
    ax.set_xticklabels([f"{int(x)}" for x in w])
    ax.set_xlabel("Number of parfor workers")
    ax.set_ylabel("Parallel efficiency (%)")
    ax.set_ylim(0, 118)
    polish(ax)

    # (c) 12-frame E. coli 3-plane flow-field reconstruction stacked wall time
    ax = axes[2]
    xf = np.arange(len(flw))
    f_t = flw["force_s"].to_numpy(dtype=float)
    v_t = flw["velocity_s"].to_numpy(dtype=float)
    s_t = flw["streamlines_s"].to_numpy(dtype=float)
    tot = flw["total_s"].to_numpy(dtype=float)
    sp_f = flw["speedup_vs_serial"].to_numpy(dtype=float)

    bw = 0.50
    ax.bar(xf, f_t, width=bw, color=GREEN, label="Force extraction")
    ax.bar(xf, v_t, width=bw, bottom=f_t, color=BLUE, label="3-plane velocity")
    ax.bar(xf, s_t, width=bw, bottom=f_t + v_t, color=ORANGE, label="2-D streamlines")
    for idx, (tv, sv) in enumerate(zip(tot, sp_f)):
        txt = f"{tv:.2f} s" if idx == 0 else f"{tv:.2f} s\n({sv:.2f}x)"
        ax.annotate(
            txt,
            xy=(idx, tv),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=6.0,
        )
    ax.set_xticks(xf)
    ax.set_xticklabels(["CPU\nserial", "CPU\nThreads (8)", "GPU FP64 +\nThreads (8)"], fontsize=6.4)
    ax.set_ylabel("Wall-clock time, 12 frames (s)")
    ax.set_ylim(0, tot.max() * 1.28)
    ax.legend(loc="upper right", handlelength=1.2, borderaxespad=0.2)
    ax.text(
        0.97,
        0.54,
        "GPU vs CPU FP64:\nmax rel err 4.3e-14",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.9,
        color=GREY,
    )
    polish(ax)

    add_panel_labels(axes, x_pt=-15.0, y_pt=3.0)
    save_all(fig, "fig_parallel_and_gpu")


# ---------------------------------------------------------------------------
# Figure 10: Application 1 — Cyanobacterial Pili Sedimentation & Shear (190 mm)
# ---------------------------------------------------------------------------
def fig_app_cyanobacteria() -> None:
    arr = pd.read_csv(PAPER_DIR / "sed_arrangement_length.csv")
    flx = pd.read_csv(PAPER_DIR / "sed_flexibility.csv")
    shr = pd.read_csv(PAPER_DIR / "sed_shear.csv")

    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(W190, 64 * MM), layout="constrained")

    # (a) Arrangement & pilus length drag reduction
    ax = axes[0]
    arr_styles = [
        ("mono", "1 polar", BLUE, "o", "-"),
        ("amphi", "2 polar", SKY, "s", "-"),
        ("lat2", "2 lateral", GREEN, "^", "-"),
        ("lat4", "4 lateral", ORANGE, "D", "-"),
        ("hemi4", "4 asymmetric", VERMILION, "v", "--"),
        ("peri8", "8 around caps", PINK, "P", "-"),
    ]
    ax.axhline(1.0, color=GREY, ls="--", lw=0.6, zorder=0)
    for key, lbl, col, mk, ls in arr_styles:
        sub = arr[arr["arrangement"] == key].sort_values("length_um")
        if len(sub) > 0:
            ax.plot(
                sub["length_um"],
                sub["V_over_no_filament"],
                color=col,
                marker=mk,
                ls=ls,
                ms=3.8,
                lw=0.95,
                label=lbl,
            )
    ax.set_xlabel("Pilus length (µm)")
    ax.set_ylabel(r"$V_{\mathrm{sink}}\,/\,V_{\mathrm{bare}}$")
    ax.set_xlim(0.5, 10.6)
    ax.set_ylim(0.14, 1.07)
    ax.legend(loc="lower left", ncol=2, columnspacing=0.7, handlelength=1.3, fontsize=5.8)
    polish(ax)

    # (b) Flexibility-induced sedimentation speedup
    ax = axes[1]
    ax.axhline(1.0, color=GREY, ls="--", lw=0.6, zorder=0)
    for key, lbl, col, mk in [
        ("hemi4", "4 asymmetric (3 µm)", VERMILION, "o"),
        ("lat2", "2 lateral (3 µm)", GREEN, "s"),
    ]:
        sub = flx[flx["arrangement"] == key].sort_values("E_Pa")
        ax.plot(
            sub["E_Pa"],
            sub["V_over_default_E"],
            color=col,
            marker=mk,
            ms=4.0,
            lw=1.0,
            label=lbl,
        )
    ax.axvline(3.279e8, color=GREY, ls=":", lw=0.7, zorder=0)
    ax.text(
        3.279e8 * 0.72,
        1.53,
        "Default\n328 MPa",
        ha="right",
        va="top",
        fontsize=5.9,
        color=GREY,
    )
    ax.set_xscale("log")
    ax.set_xlabel("Young's modulus $E$ (Pa)")
    ax.set_ylabel(r"$V_{\mathrm{sink}}\,/\,V_{\mathrm{sink}}(E_{\mathrm{def}})$")
    ax.set_ylim(0.94, 1.68)
    ax.legend(loc="center left", bbox_to_anchor=(0.01, 0.25), handlelength=1.4)
    polish(ax)

    # (c) Shear-flow sedimentation
    ax = axes[2]
    for key, lbl, col, mk in [
        ("none", "Bare capsule", GREY, "s"),
        ("mono", "1 polar pilus", BLUE, "o"),
        ("hemi4", "4 asymmetric pili", VERMILION, "^"),
    ]:
        sub = shr[shr["arrangement"] == key].sort_values("shear_rate_1_s")
        ax.plot(
            sub["shear_rate_1_s"],
            sub["V_sink_um_s"],
            color=col,
            marker=mk,
            ms=4.2,
            lw=1.0,
            label=lbl,
        )
    ax.set_xscale("log")
    ax.set_xlabel("Shear rate (s$^{-1}$)")
    ax.set_ylabel(r"Mean $V_{\mathrm{sink}}$ (µm s$^{-1}$)")
    ax.set_ylim(0.058, 0.172)
    ax.legend(loc="lower left", handlelength=1.4)
    ax.text(
        0.97,
        0.94,
        r"$\Delta X / \Delta X_{\mathrm{fluid}} \in [1.00, 1.08]$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.0,
        color=GREY,
    )
    polish(ax)

    add_panel_labels(axes, x_pt=-15.0, y_pt=3.0)
    save_all(fig, "fig_app_cyanobacteria")


# ---------------------------------------------------------------------------
# Figure 11: Application 2 — E. coli Single-Flagellum Motor Swimmer (190 mm)
# ---------------------------------------------------------------------------
def fig_app_ecoli() -> None:
    mat = sio.loadmat(CPC_DIR / "app_ecoli_chlamy_plotdata.mat", squeeze_me=True, struct_as_record=False)
    ec = mat["ec"]
    t = np.asarray(ec.t, dtype=float)
    X = np.asarray(ec.X_um, dtype=float)
    om_b = np.asarray(ec.om_body_Hz, dtype=float)
    om_rel = np.asarray(ec.om_rel_Hz, dtype=float)

    apply_style()
    fig, axes = plt.subplots(
        1, 3, figsize=(W190, 64 * MM), layout="constrained",
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.25]},
    )

    # (a) Forward swimming X(t) + lateral wiggling zoom inset
    ax = axes[0]
    ax.plot(t, X[:, 0], color=BLUE, lw=1.1, label=r"$X(t)$ (17.24 µm s$^{-1}$)")
    ax.set_xlabel("Time $t$ (s)")
    ax.set_ylabel(r"Forward position $X$ (µm)")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.5, 19.0)
    ax.legend(loc="upper left", handlelength=1.3)
    polish(ax)

    ins = ax.inset_axes([0.52, 0.14, 0.44, 0.39])
    m_z = (t >= 0.40) & (t <= 0.46)
    ins.plot(t[m_z] * 1e3, (X[m_z, 1] - X[m_z, 1].mean()) * 1e3, color=VERMILION, lw=0.75, label="Y")
    ins.plot(t[m_z] * 1e3, (X[m_z, 2] - X[m_z, 2].mean()) * 1e3, color=GREEN, lw=0.75, ls="--", label="Z")
    ins.set_xlabel("t (ms)", fontsize=5.5, labelpad=3)
    ins.set_ylabel("ΔY, ΔZ (nm)", fontsize=5.5, labelpad=8)
    ins.tick_params(direction="in", labelsize=5.3, length=1.8, width=0.4, top=True, right=True)
    ins.legend(loc="upper right", fontsize=5.0, handlelength=1.0, borderaxespad=0.1)

    # (b) Self-consistent motor frequency & torque-free body counter-rotation
    ax = axes[1]
    t_mid = 0.5 * (t[:-1] + t[1:])
    ax.plot(t[2:-2], om_rel[2:-2], color=BLUE, lw=1.0, label=r"Motor $\omega_{\mathrm{rel}}$ (170.4 Hz)")
    ax.plot(t_mid, np.abs(om_b), color=VERMILION, lw=1.0, label=r"Body $|\Omega_{\mathrm{b}}|/2\pi$ (14.8 Hz)")
    ax.axhline(300.0, color=GREY, ls=":", lw=0.7, label=r"Zero-load $\omega_{\max}$ (300 Hz)")
    ax.set_xlabel("Time $t$ (s)")
    ax.set_ylabel("Rotation frequency (Hz)")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 330.0)
    ax.legend(loc="upper right", bbox_to_anchor=(0.99, 0.87), handlelength=1.3, fontsize=5.9)
    polish(ax)

    # (c) Reconstructed XZ mid-plane flow field + streamlines
    ax = axes[2]
    hv = np.asarray(ec.xz_h, dtype=float)
    vv = np.asarray(ec.xz_v, dtype=float)
    Umag = np.asarray(ec.xz_Umag, dtype=float)
    cm = np.asarray(ec.Xcm_fin_um, dtype=float)
    # Shift coordinates relative to final CM
    h_rel = hv - cm[0]
    v_rel = vv - cm[2]
    valid = Umag[np.isfinite(Umag)]
    vmax = float(np.percentile(valid, 99.0))
    im = ax.pcolormesh(
        h_rel,
        v_rel,
        Umag,
        cmap="cividis",
        vmin=0.0,
        vmax=vmax,
        shading="auto",
        rasterized=True,
    )
    # Streamlines
    sx = np.asarray(ec.xz_sx, dtype=float) - cm[0]
    sy = np.asarray(ec.xz_sy, dtype=float) - cm[2]
    ax.plot(sx, sy, color="white", lw=0.45, alpha=0.82, zorder=3)
    # Arrowheads
    av = np.asarray(ec.xz_av, dtype=float)
    af = np.asarray(ec.xz_af, dtype=int) - 1
    if av.ndim == 2 and af.ndim == 2 and len(af) > 0:
        av_s = np.column_stack([av[:, 0] - cm[0], av[:, 1] - cm[2]])
        polys = av_s[af]
        ax.add_collection(PolyCollection(polys, facecolor="white", edgecolor="none", alpha=0.85, zorder=3.5))
    # Body section & flagellum projection
    cap_h = np.asarray(ec.xz_cap_h, dtype=float) - cm[0]
    cap_v = np.asarray(ec.xz_cap_v, dtype=float) - cm[2]
    ax.fill(cap_h, cap_v, facecolor="#D9D9D9", edgecolor="black", lw=0.65, zorder=4)
    flg = np.asarray(ec.flag_um, dtype=float)
    ax.plot(flg[:, 0] - cm[0], flg[:, 2] - cm[2], color="#FF5500", lw=1.05, zorder=5)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(float(h_rel.min()), float(h_rel.max()))
    ax.set_ylim(float(v_rel.min()), float(v_rel.max()))
    ax.set_xlabel(r"$x - X_{\mathrm{cm}}$ (µm)")
    ax.set_ylabel(r"$z - Z_{\mathrm{cm}}$ (µm)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.048, pad=0.025)
    cbar.set_label(r"$|\mathbf{u}|$ (µm s$^{-1}$)", fontsize=7.0)
    cbar.ax.tick_params(labelsize=6.2, length=2.0, width=0.45)
    ax.tick_params(direction="in", top=True, right=True, length=2.4, width=0.5)

    add_panel_labels(axes, x_pt=-15.0, y_pt=3.0)
    save_all(fig, "fig_app_ecoli")


# ---------------------------------------------------------------------------
# Figure 12: Application 3 — Chlamydomonas 3-D Helical Swimmer (190 mm)
# ---------------------------------------------------------------------------
def fig_app_chlamy() -> None:
    mat = sio.loadmat(CPC_DIR / "app_ecoli_chlamy_plotdata.mat", squeeze_me=True, struct_as_record=False)
    ch = mat["ch"]
    oop = pd.read_csv(ROOT / "output" / "chlamy_oop" / "oop_summary.csv")

    f1 = np.asarray(ch.f1_sn, dtype=float)  # [n1 x 3 x nsnap]
    f2 = np.asarray(ch.f2_sn, dtype=float)
    tip1 = np.asarray(ch.tip1, dtype=float)  # [3 x nt]
    tip2 = np.asarray(ch.tip2, dtype=float)
    X_tr = np.asarray(ch.X_um, dtype=float)

    apply_style()
    fig = plt.figure(figsize=(W190, 66 * MM), layout="constrained")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 0.95, 1.15])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2], projection="3d")

    # (a) Body-frame 3-D breaststroke snapshots + non-planar tip excursion
    # Draw prolate spheroid body outline (c=5 um, a=4 um)
    ang = np.linspace(0, 2 * np.pi, 120)
    ax0.fill(5.0 * np.cos(ang), 4.0 * np.sin(ang), facecolor="#E5EFE5", edgecolor=GREEN, lw=0.8, zorder=1)
    nsnap = f1.shape[2]
    shades = plt.cm.Blues(np.linspace(0.32, 0.95, nsnap))
    for j in range(nsnap):
        ax0.plot(f1[:, 0, j], f1[:, 1, j], color=shades[j], lw=0.85, zorder=2)
        ax0.plot(f2[:, 0, j], f2[:, 1, j], color=shades[j], lw=0.85, zorder=2)
    # Highlight one clean beat cycle of the closed tip orbit in XY
    sl_beat = slice(480, 500)
    ax0.plot(tip1[0, sl_beat], tip1[1, sl_beat], color=VERMILION, lw=1.0, ls="-", label="Tip orbit (1 beat)", zorder=3)
    ax0.plot(tip2[0, sl_beat], tip2[1, sl_beat], color=VERMILION, lw=1.0, ls="-", zorder=3)
    ax0.set_aspect("equal", adjustable="box")
    ax0.set_xlim(-6.5, 13.8)
    ax0.set_ylim(-12.5, 12.5)
    ax0.set_xlabel(r"Body $x$ (µm)")
    ax0.set_ylabel(r"Body $y$ (µm)")
    ax0.legend(loc="lower left", handlelength=1.2, fontsize=5.8)
    ax0.text(
        0.04,
        0.96,
        "1 beat (57.8 Hz)\ntip Δz = 3.4 µm",
        transform=ax0.transAxes,
        ha="left",
        va="top",
        fontsize=5.8,
        color=BLACK,
    )
    polish(ax0)

    # (b) Chiral 3-D spin superposition vs out-of-plane ratio r_perp
    sweep = oop[oop["tag"].isin(["OOP0", "OOP05", "OOP10", "OOP15", "OOP20"])].sort_values("r_perp")
    ctrl_rev = oop[oop["tag"] == "OOP10M"].iloc[0]
    ctrl_flat = oop[oop["tag"] == "OOP10F"].iloc[0]
    ctrl_sym = oop[oop["tag"] == "OOP10S"].iloc[0]

    ax1.axhspan(1.0, 2.0, color="#E9F5EE", zorder=0, label="Live-cell band (1–2 Hz)")
    ax1.plot(
        sweep["r_perp"],
        np.abs(sweep["spin_Hz"]),
        color=BLUE,
        marker="o",
        ms=4.2,
        lw=1.05,
        label="Tilt 17° + OOP [+1, −1]",
    )
    ax1.plot([0.25], [1.060], color=VERMILION, marker="*", ms=6.8, lw=0, label="HELIX3D_2S (2.0 s)")
    ax1.plot([float(ctrl_rev["r_perp"])], [abs(float(ctrl_rev["spin_Hz"]))],
             color=ORANGE, marker="v", ms=4.5, lw=0, label="Reversed [−1, +1]")
    ax1.plot([float(ctrl_flat["r_perp"])], [abs(float(ctrl_flat["spin_Hz"]))],
             color=GREY, marker="s", ms=4.0, lw=0, label="Zero tilt + OOP")
    ax1.plot([float(ctrl_sym["r_perp"])], [abs(float(ctrl_sym["spin_Hz"]))],
             color=PINK, marker="D", ms=3.8, lw=0, label="Same-sign [+1, +1]")
    ax1.set_xlabel(r"Out-of-plane ratio $r_\perp$")
    ax1.set_ylabel(r"Axial spin rate $|f_{\mathrm{spin}}|$ (Hz)")
    ax1.set_xlim(-0.015, 0.275)
    ax1.set_ylim(0.0, 1.72)
    ax1.legend(loc="upper left", ncol=1, fontsize=5.3, handlelength=1.1, borderaxespad=0.2)
    polish(ax1)

    # (c) 2.0 s 3-D left-handed superhelical trajectory
    ax2.plot(X_tr[:, 0], X_tr[:, 1], X_tr[:, 2], color=BLUE, lw=1.05)
    ax2.scatter([X_tr[0, 0]], [X_tr[0, 1]], [X_tr[0, 2]], color=GREEN, s=16, zorder=5, label="t = 0 s")
    ax2.scatter([X_tr[-1, 0]], [X_tr[-1, 1]], [X_tr[-1, 2]], color=VERMILION, s=18, zorder=5, label="t = 2.0 s")
    # Project trajectory onto bottom plane for crisp depth cue
    z_min = float(X_tr[:, 2].min()) - 2.0
    ax2.plot(X_tr[:, 0], X_tr[:, 1], np.full_like(X_tr[:, 2], z_min), color=GREY, lw=0.5, alpha=0.45)
    ax2.set_xlabel("X (µm)", fontsize=6.6, labelpad=0.5)
    ax2.set_ylabel("Y (µm)", fontsize=6.6, labelpad=0.5)
    ax2.set_zlabel("Z (µm)", fontsize=6.6, labelpad=-3.0)
    ax2.tick_params(labelsize=5.6, pad=0.0)
    ax2.view_init(elev=22, azim=-58)
    ax2.set_box_aspect((1.05, 0.95, 0.78), zoom=0.86)
    ax2.set_zlim(z_min, float(X_tr[:, 2].max()) + 1.5)
    ax2.legend(loc="upper left", fontsize=5.6, handletextpad=0.3)
    ax2.text2D(
        0.98,
        0.97,
        "U = 77.2 µm/s, R = 1.54 µm\nP = −71.5 µm (LH), 2.07 turns",
        transform=ax2.transAxes,
        ha="right",
        va="top",
        fontsize=5.6,
        color=BLACK,
    )

    add_panel_labels([ax0, ax1], x_pt=-15.0, y_pt=3.0)
    ax2.text2D(0.0, 1.01, "(c)", transform=ax2.transAxes, fontsize=8, fontweight="bold", ha="left", va="bottom")
    save_all(fig, "fig_app_chlamy")


def main() -> None:
    fig_architecture()
    fig_parallel_and_gpu()
    fig_app_cyanobacteria()
    fig_app_ecoli()
    fig_app_chlamy()


if __name__ == "__main__":
    main()
