"""Clamped-free bending modes and their share of transverse variance.

The filament is a straight, uniform Euler-Bernoulli beam, clamped at the
root and free at the tip.  Its bending modes are the eigenfunctions of
that beam.  With the elastic energy quadratic in curvature, each mode
holds kT/2, but the mean-square transverse deflection weights mode n by
1/(beta_n L)^4.  In the continuum normalisation the share of the
arc-length integral of <y^2> and the share of the tip variance are the
same number,

    f_n = 12 / (beta_n L)^4,

where beta_n L is the n-th root of cos(beta L) cosh(beta L) = -1.
The two polarisations carry the same split, so the shares also apply to
y^2 + z^2.  Hydrodynamics sets the relaxation time, not this static share.

The curves in panel a are scaled to a common peak deflection so that the
higher orders stay visible.  The variance shares are the bars in panel b.
Orders 1-5 are drawn separately; order 6 and above are one remainder bar.
The vertical axis of panel b is logarithmic because the shares span three
decades.  These are exact modal shares, not a histogram of a sample.

Panel c is not the free-beam scaling.  It is the elastohydrodynamic
relaxation time of one 3 um filament on the example capsule (VAL-ein-pili1,
30 segments, N_body = 300), from init_pili_brownian_modes.  The two
polarisations agree.  The times are read from
output/pilus_bending_modes/tau_pili1_orders5.mat, written by
pilus_mode_relaxation_pili1.m.

Author: Lucien (lucien-6@qq.com)
Created: 2026-11-19
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_nature_figures import (  # noqa: E402
    BLACK,
    BLUE,
    GREEN,
    GREY,
    ORANGE,
    PINK,
    VERMILION,
    W183,
    apply_style,
    panel_letters,
    save_figure,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "pilus_bending_modes"

N_SHOW = 5
N_SUM = 40
MODE_COLORS = (BLUE, ORANGE, GREEN, VERMILION, PINK)
BAR_GREY = "#4D4D4D"


def beta_roots(n_modes: int) -> np.ndarray:
    """First n_modes roots of cos(x) cosh(x) + 1 = 0, x = beta L.

    Odd orders lie just above (2n-1) pi/2 and even orders just below it.
    Newton is started from that shift, sech((2n-1) pi/2), and the residual
    cos(x) + sech(x) is used so the hyperbolic functions do not overflow.
    """
    roots = np.empty(n_modes)
    for n in range(1, n_modes + 1):
        x0 = (2 * n - 1) * np.pi / 2
        sign = 1.0 if n % 2 == 1 else -1.0
        x = x0 + sign * _sech(x0)
        for _ in range(10):
            slope = -np.sin(x) - _sech(x) * _tanh(x)
            step = _residual(x) / slope
            x = x - step
            if abs(step) < 1e-14 * x:
                break
        roots[n - 1] = x
    return roots


def _sech(x: float) -> float:
    if x > 40.0:
        return 0.0
    if x > 20.0:
        return float(2.0 * np.exp(-x) / (1.0 + np.exp(-2.0 * x)))
    return float(1.0 / np.cosh(x))


def _tanh(x: float) -> float:
    if x > 20.0:
        return 1.0
    return float(np.tanh(x))


def _residual(x: float) -> float:
    """cos(x) + sech(x), zero at a clamped-free eigenvalue."""
    return float(np.cos(x) + _sech(x))


def variance_shares(betas: np.ndarray) -> np.ndarray:
    """Share of transverse mean-square deflection, f_n = 12 / (beta L)^4."""
    shares = 12.0 / betas**4
    total = float(np.sum(shares))
    if abs(total - 1.0) > 1e-3:
        raise RuntimeError(f"modal shares sum to {total}, not 1")
    return shares / total


def mode_shape(beta: float, s: np.ndarray) -> np.ndarray:
    """Clamped-free eigenfunction on s/L in [0, 1].

    The shape satisfies y = y' = 0 at the root.  It is not normalised:
    callers scale it for display.
    """
    e = np.exp(-beta)
    sigma = (1.0 - e * e - 2.0 * e * np.sin(beta)) / (
        1.0 + e * e + 2.0 * e * np.cos(beta)
    )
    x = beta * s
    return (np.cosh(x) - np.cos(x)) - sigma * (np.sinh(x) - np.sin(x))


def _percent_label(share: float) -> str:
    percent = 100.0 * share
    if percent >= 1.0:
        return f"{percent:.2f}%"
    return f"{percent:.3f}%"


def load_relaxation(path: Path) -> np.ndarray:
    """Mean relaxation time of each bending order, in seconds.

    tau_by_order is 5 x 2, one column per polarisation.
    """
    record = sio.loadmat(path)
    tau = np.asarray(record["tau_by_order"], dtype=float)
    if tau.shape != (N_SHOW, 2):
        raise ValueError(f"expected tau_by_order shape {(N_SHOW, 2)}, got {tau.shape}")
    return tau.mean(axis=1)


def _time_label(tau: float) -> str:
    if tau >= 0.1:
        return f"{tau:.3f} s"
    milliseconds = tau * 1.0e3
    if milliseconds >= 10.0:
        return f"{milliseconds:.1f} ms"
    if milliseconds >= 1.0:
        return f"{milliseconds:.2f} ms"
    return f"{milliseconds:.3f} ms"


def draw(betas: np.ndarray, shares: np.ndarray, tau_s: np.ndarray) -> plt.Figure:
    """Panel a: mode shapes.  Panel b: variance shares.  Panel c: relaxation times.

    Each mode is drawn in its own row, scaled to unit peak deflection.
    The thickened left spine is the clamp.  The tip, at s/L = 1, is free.
    """
    apply_style()
    fig = plt.figure(figsize=(W183, 118.0 / 25.4), layout="constrained")
    outer = fig.add_gridspec(1, 2, width_ratios=(1.05, 1.15))
    mode_grid = outer[0, 0].subgridspec(N_SHOW, 1)
    bar_grid = outer[0, 1].subgridspec(2, 1)
    mode_axes = []
    s = np.linspace(0.0, 1.0, 500)
    for i in range(N_SHOW):
        ax = fig.add_subplot(mode_grid[i, 0], sharex=mode_axes[0] if mode_axes else None)
        mode_axes.append(ax)
        shape = mode_shape(float(betas[i]), s)
        shape = shape / np.max(np.abs(shape))
        ax.axhline(0.0, color=GREY, lw=0.45, zorder=1)
        ax.plot(s, shape, color=MODE_COLORS[i], lw=0.95, zorder=2)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(-1.35, 1.55)
        ax.set_yticks([])
        ax.spines["left"].set_linewidth(1.15)
        ax.tick_params(direction="in", top=True, right=True, length=2.2, width=0.5,
                       labelbottom=(i == N_SHOW - 1))
        ax.text(
            0.03, 0.78, str(i + 1), transform=ax.transAxes,
            ha="left", va="center", fontsize=7, fontfamily="Arial",
            color=BLACK, clip_on=False,
        )
    mode_axes[-1].set_xlabel("s / L")
    mode_axes[N_SHOW // 2].set_ylabel("mode shape")

    shown = shares[:N_SHOW]
    tail = float(np.sum(shares[N_SHOW:]))
    bar_shares = np.append(shown, tail)
    labels = [str(i) for i in range(1, N_SHOW + 1)] + ["> 5"]
    colors = list(MODE_COLORS) + [BAR_GREY]
    bx = fig.add_subplot(bar_grid[0, 0])
    x = np.arange(len(bar_shares))
    bx.bar(x, bar_shares, width=0.72, color=colors, zorder=2)
    # Orders 5 and >5 differ by about 7 %.  Upright labels would occupy the
    # same height, so those two are rotated to stay on their own bars.
    for i, (xpos, share) in enumerate(zip(x, bar_shares)):
        upright = i < len(bar_shares) - 2
        bx.text(
            xpos, share * (1.45 if upright else 1.9), _percent_label(float(share)),
            ha="center", va="bottom", rotation=0 if upright else 90,
            fontsize=6.5, fontfamily="Arial", color=BLACK, clip_on=False,
        )
    bx.set_yscale("log")
    bx.set_xlim(-0.55, len(bar_shares) - 0.45)
    bx.set_ylim(8.0e-5, 4.0)
    bx.set_xticks(x)
    bx.set_xticklabels(labels)
    bx.set_xlabel("bending order")
    bx.set_ylabel("fraction of transverse variance")
    bx.tick_params(direction="in", top=True, right=True, which="both",
                   length=2.4, width=0.5)
    bx.tick_params(which="minor", length=1.4, width=0.4)

    tx = fig.add_subplot(bar_grid[1, 0])
    xt = np.arange(N_SHOW)
    tx.bar(xt, tau_s, width=0.72, color=MODE_COLORS, zorder=2)
    for xpos, tau in zip(xt, tau_s):
        tx.text(
            xpos, tau * 1.45, _time_label(float(tau)),
            ha="center", va="bottom", fontsize=6.5, fontfamily="Arial",
            color=BLACK, clip_on=False,
        )
    tx.set_yscale("log")
    tx.set_xlim(-0.55, N_SHOW - 0.45)
    tx.set_ylim(float(np.min(tau_s)) / 8.0, float(np.max(tau_s)) * 6.0)
    tx.set_xticks(xt)
    tx.set_xticklabels([str(i) for i in range(1, N_SHOW + 1)])
    tx.set_xlabel("bending order")
    tx.set_ylabel("relaxation time (s)")
    tx.tick_params(direction="in", top=True, right=True, which="both",
                   length=2.4, width=0.5)
    tx.tick_params(which="minor", length=1.4, width=0.4)

    panel_letters([mode_axes[0], bx, tx], "abc")
    return fig


def main() -> None:
    betas = beta_roots(N_SUM)
    shares = variance_shares(betas)
    print("beta L and percent of transverse variance")
    for n in range(N_SHOW):
        print(f"  order {n + 1}: beta L = {betas[n]:.6f}   {100 * shares[n]:.4f}%")
    print(f"  orders > 5: {100 * float(np.sum(shares[N_SHOW:])):.4f}%")
    tau_s = load_relaxation(OUT / "tau_pili1_orders5.mat")
    print("relaxation time, one filament")
    for n, tau in enumerate(tau_s, start=1):
        print(f"  order {n}: {tau:.6g} s")
    fig = draw(betas, shares, tau_s)
    save_figure(fig, OUT / "fig_clamped_bending_modes")
    print(f"wrote {OUT / 'fig_clamped_bending_modes'}.png / .pdf")


if __name__ == "__main__":
    main()
