# Chlamydomonas Biflagellate 3-D Swimming Module

English | [中文](README.zh-CN.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Created:** 2026-11-24  
> **Applies to:** `Bacteria Elastohydrodynamic Simulation` (Kirchhoff rod + RSS/RSM coupled solver)

This module models free-swimming biflagellate *Chlamydomonas reinhardtii* (wild-type CC-124 / CC-125) in low-Reynolds-number Newtonian fluid using the repository's coupled regularized Stokeslet / Kirchhoff elastic rod engine. It upgrades the measured-gait planar breaststroke into a **calibrated, physically grounded 3-D helical swimmer** that reproduces live-cell axial spin (`1–2 Hz`), left-handed superhelical trajectories, 3-D non-planar tip orbits, and proximal flagellar 3-D torsion without body penetration.

---

## 1. Directory Structure & Integration with the Main Solver

The module separates tracked input data, reusable physics engines, one-click case scripts, and plotting tools, matching the top-level layout of the repository (`configs/`, `geometry/`, `src/`, `output/`):

```text
<project_root>/
├── configs/
│   └── config_chlamy_biflagellate_gz.json   # Authoritative biflagellate JSON (5×4 µm body, 13° roots)
├── geometry/
│   ├── chlamy_body_spheroid.csv             # Prolate spheroid body surface (400 Riesz points)
│   ├── chlamy_flagellum_ccw_gz.csv          # +Y flagellum measured mean rest shape (L = 11.70 µm)
│   ├── chlamy_flagellum_cw_gz.csv           # -Y flagellum mirror rest shape
│   └── generate_chlamydomonas_geometry.m    # Geometry generator (measured 1 cP mean curvature)
└── chlamydomonas/
    ├── README.md                            # This document (English)
    ├── README.zh-CN.md                      # Chinese documentation
    ├── chlamy_run_preset.m                  # Main preset entry point (3-D & planar controls)
    ├── data/
    │   └── fit_1p00A.mat                    # Tracked dual-k 1ω+2ω measured-gait moment table (57.75 Hz)
    ├── core/
    │   ├── chlamy_build_params.m            # Parameter & initial-state builder (mirrors run_simulation)
    │   ├── chlamy_inverse_moment.m          # Inverse EHD solver: measured gait -> active moment q(s,t)
    │   ├── chlamy_fit_moment.m              # Analytic dual-k 1ω+2ω harmonic fit of q(s,t)
    │   ├── chlamy_make_perp_envelope.m      # 3-D proximal 1ω+2ω out-of-plane moment builder
    │   ├── chlamy_table_moment_fun.m        # Custom active-moment handle factory (value-captured par)
    │   ├── chlamy_table_moment_eval.m       # Stand-alone moment evaluator (zero parasitic twist)
    │   ├── chlamy_run_case.m                # General headless case runner with parameter overrides
    │   └── chlamy_swim_analysis.m           # Unwrapped SO(3) kinematics, helix & 3-D orbit post-check
    ├── cases/
    │   ├── run_3d_helix_2s.m                # One-click realistic 2.0 s 3-D helical swimming (~116 beats)
    │   ├── run_planar_sym50.m               # One-click coplanar symmetric reference (50 beats)
    │   └── run_planar_a75f75.m              # One-click coplanar amp/freq asymmetry control (1.0 s)
    ├── plotting/
    │   ├── plot_measured_gait_1p00A.py      # Plot measured experimental gait of cell 1p00 A
    │   └── plot_gait_comparison_sym50.py    # Compare SYM50 simulation against measured gait
    └── Refers/                              # Reference papers & Dryad experimental waveforms
```

**How it plugs into `src/`:** The main solver stays species-agnostic. `chlamy_run_case.m` loads `configs/config_chlamy_biflagellate_gz.json`, applies optional beat-plane tilt via `custom_geometry.pili(k).deflection_z`, and replaces the JSON placeholder `active_moment` with `type = 'custom'` handles pointing to `chlamy_table_moment_eval.m`. Saved `simulation_results.mat` files remain self-contained and callable in fresh sessions via `src/visualization/private/restore_saved_moment_handles.m`.

---

## 2. Quick Start

Run all commands from the project root (or open any script in `chlamydomonas/cases/` and click **Run** in MATLAB; scripts resolve the project root automatically).

### Option A: One-Click Case Scripts

```matlab
% Realistic 2.0 s 3-D helical swimming (~116 beats, ~2.2 helical turns, MP4 animation)
% Output: output/Chlamy-3D-Helix-2s/
run('chlamydomonas/cases/run_3d_helix_2s.m')

% Coplanar symmetric 50-beat reference
run('chlamydomonas/cases/run_planar_sym50.m')
```

### Option B: Unified Preset Interface (`chlamy_run_preset`)

```matlab
addpath('chlamydomonas');

% Fast 3-beat high-precision benchmarks (tol = 1e-4)
s_sym = chlamy_run_preset('SYM3D');       % C2-symmetric 3-D breaststroke (straight + ~1.17 Hz spin)
s_hlx = chlamy_run_preset('HELIX3D');     % 5% cis/trans asymmetry (left-handed superhelix)

% Long runs with 3-D MP4 animation
s_20b = chlamy_run_preset('HELIX3D_20');  % 20 beats (~0.35 s)
s_2s  = chlamy_run_preset('HELIX3D_2S');  % 2.0 s (~116 beats, ~2.2 full helical turns)

% Dry-run validation (checks geometry, handles, and initial ODE RHS in ~3 s without integrating)
diag = chlamy_run_preset('HELIX3D_2S', [], struct('dry_run', true));
```

### Option C: Custom Parameter Sweeps (`chlamy_run_case`)

```matlab
addpath('chlamydomonas/core');
S = load('chlamydomonas/data/fit_1p00A.mat', 'fit');
tab3d = chlamy_make_perp_envelope(S.fit, struct( ...
    'r_perp', 0.25, 'r_perp2', 0.05, 's_w', 0.20, ...
    'delta_perp', pi/2, 'delta_perp2', pi/2));

ov = struct( ...
    'case_name',      'My-Chlamy-3D-Test', ...
    'moment_table',   tab3d, ...
    'tilt_deg',       [17, 17], ...
    'perp_sign',      [1, -1], ...
    'drive_amp',      [1, 0.95], ...
    'T_obs',          3 / S.fit.f_Hz, ...
    'rel_tol',        1e-4, ...
    'abs_tol',        1e-4);

[results, params, out_dir] = chlamy_run_case('configs/config_chlamy_biflagellate_gz.json', ov);
s = chlamy_swim_analysis(results, struct('f_Hz', S.fit.f_Hz, 'n_skip', 1));
```

### Monitoring Live Progress During Long Runs

Every run automatically writes a lightweight heartbeat every 10 s to `output/<case_name>/progress.txt`:

```text
pct=42.50 t=1.4536e+01 T=3.4203e+01 elapsed=14820.1 eta=20050.7 seg=1
```

In PowerShell, inspect live percentage and ETA at any time without interrupting MATLAB:

```powershell
Get-Content output\Chlamy-3D-Helix-2s\progress.txt
```

---

## 3. Physical & Mathematical Formulation (3-D Beating Model)

### 3.1 Geometry & Constitutive Law

- **Cell body:** Prolate spheroid, semi-major axis \(c = 5.0\,\mu\text{m}\) along body \(+X\) (swimming direction), semi-minor axis \(a = 4.0\,\mu\text{m}\) along \(Y, Z\), discretized with 400 cached Riesz points.
- **Flagella:** Two clamped Kirchhoff rods of length \(L_k = 11.70\,\mu\text{m}\) (\(N_{\mathrm{seg}} = 71\) at `seg_per_um = 6`), hydrodynamic radius \(r = 0.125\,\mu\text{m}\), bending stiffness \(EI = 8.40 \times 10^{-22}\,\text{N}\cdot\text{m}^2\) (\(EI/L_k = 7.18 \times 10^{-17}\,\text{N}\cdot\text{m}\)). Attachment sites sit at polar angle \(13^\circ\) from \(+X\) (basal-body separation \(1.80\,\mu\text{m}\)) with root tangent half-angle \(13^\circ\).
- **Rest shape:** Time-averaged curvature profile measured in swimming cells in water (Qin et al., *Sci. Rep.* 2015; Gutierrez et al., *Phil. Trans. R. Soc. A* 2025; total static turning \(\approx 124^\circ\)), held as `intrinsic_curvature = 'frozen'`.
- **Constitutive sign convention:** Internal rod residual is \(R_{\mathrm{rod}} = \kappa_{\mathrm{dev}} - \mathbf{M}_{\mathrm{act}}\).

### 3.2 Two-Component Orthogonal Active Moment

<!-- markdownlint-disable MD049 -->
At balance node \(j \in \{2,\dots,N_{\mathrm{seg}}\}\) (dimensionless arc length \(s = (j-1)/N_{\mathrm{seg}}\), root node \(j=1\) zeroed as clamped reaction), the active internal moment vector on flagellum \(k \in \{1, 2\}\) is:

\[
\mathbf{M}_k(s, t) = q_{\parallel,k}(s, t)\,\hat{\mathbf{n}}_k(t) + q_{\perp,k}(s, t)\,\hat{\mathbf{b}}_k(s, t)
\]

1. **Tilted beat-plane normal \(\hat{\mathbf{n}}_k(t)\):** Obtained by rotating the frontal-plane normal \(\hat{\mathbf{z}}\) right-handedly by \(\gamma_k\) about each flagellum's outward root tangent \(\hat{\mathbf{t}}_{0,k}\), stored in the clamped root frame and transported rigidly with the cell.
   - **Sign convention (`tilt_deg = [+17, +17]`):** Because the two flagella curl in opposite directions (`+Y` vs `-Y`), **same-sign** right-handed rotation `[+17, +17]` lifts flagellum 1 (`+Y`) toward `+Z` (`+2.89 µm` tip) and presses flagellum 2 (`-Y`) toward `-Z` (`-2.89 µm` tip), realizing the exact \(C_2\)-chiral orbital tilt (\(\beta \approx 17^\circ\)) reported by Cortese & Wan (*PRL* 2021). Opposite signs `[+17, -17]` lift both flagella to `+Z` (mirror-symmetric, zero net spin).
2. **Node-local out-of-plane bending axis \(\hat{\mathbf{b}}_k(s, t)\):**
   \[
   \hat{\mathbf{b}}_k(s, t) = \frac{\mathbf{d}_{3,k}(s, t) \times \hat{\mathbf{n}}_k(t)}{\|\mathbf{d}_{3,k}(s, t) \times \hat{\mathbf{n}}_k(t)\|}
   \]
   Because \(\hat{\mathbf{b}}_k \perp \mathbf{d}_{3,k}\) at every node regardless of bending angle, the out-of-plane drive produces **pure out-of-plane bending with strictly zero parasitic twist moment** (\(M_{3,k}^\perp \equiv 0\)).

### 3.3 In-Plane Dual-Wavenumber Drive \(q_{\parallel,k}(s,t)\)

Fitted to the exact inverse solution of cell `1p00 A` at \(f = 57.75\,\text{Hz}\) (`chlamydomonas/data/fit_1p00A.mat`, units of \(EI/L_k\)):

\[
q_{\parallel,k}(s, t) = (-1)^{k-1}\Big\{ c_0(s) + \alpha_k\big[A(s)\cos(\omega t + \phi(s)) + B(s)\cos(2\omega t + 2\phi(s) + \delta(s))\big] \Big\}
\]

- \(A(s) = A_0(1 - d_A s) + A_{\mathrm{pk}}\exp\!\big(-((s - s_p)/w_p)^2\big)\)
- \(\phi(s) = \phi_0 - k_1 s - (k_2 - k_1)\tau\ln\!\big(1 + \exp((s - s_b)/\tau)\big)\) *(dual wavenumber: proximal \(k_1 = 7.61\), distal \(k_2 = 3.97\,\text{rad}/L\))*
- \(B(s) = B_0 + B_{\mathrm{pk}}\exp\!\big(-((s - s_p)/w_p)^2\big)\), \(\delta(s) = d_0 + \dfrac{d_1 - d_0}{1 + \exp(-(s - s_d)/w_d)}\)
- \(c_0(s)\): 8-knot piecewise-linear static moment (unscaled by \(\alpha_k\)).

### 3.4 Out-of-Plane Drive \(q_{\perp,k}(s,t)\)

Zero static bias (\(c_{0,\perp} \equiv 0\)), phase-locked to \(\phi(s)\) with a proximal exponential envelope:

\[
q_{\perp,k}(s, t) = \sigma_k\, \mathrm{e}^{-s/s_w}\Big[ A_\perp \cos(\omega t + \phi(s) + \delta_\perp) + B_\perp \cos(2\omega t + 2\phi(s) + \delta_{\perp 2}) \Big]
\]

where \(A_\perp = r_\perp\, q_{\mathrm{rms}}\), \(B_\perp = r_{\perp 2}\, q_{\mathrm{rms}}\), and \(q_{\mathrm{rms}} = 3.421\,EI/L_k\) (\(2.46 \times 10^{-16}\,\text{N}\cdot\text{m}\)) is the oscillating RMS of \(q_{\parallel,1}\).

**Single-multiplication chirality rule:** `tab.coef_perp = {c1, c1}` stores the identical envelope for both flagella; cis/trans chirality \(\sigma = [+1, -1]\) is applied **once** at runtime via `ov.perp_sign`, guaranteeing exact \(C_2\) symmetry \(\mathbf{M}_2 = R_x(\pi)\mathbf{M}_1\) when \(\alpha_1 = \alpha_2\).
<!-- markdownlint-enable MD049 -->

### 3.5 Calibrated 3-D Parameter Summary

| Group | Parameter | Value | Meaning |
| --- | --- | --- | --- |
| **Scales** | \(f\), \(L_k\), \(EI/L_k\) | `57.75 Hz`, `11.70 µm`, `7.18e-17 N·m` | Beat frequency, flagellum length, moment unit |
| **In-plane** | \(A_0, d_A, A_{\mathrm{pk}}, s_p, w_p\) | `2.997, 0.363, 12.578, 0.097, 0.046` | 1st-harmonic amplitude profile (\(EI/L_k\)) |
| (`fit_1p00A`) | \(\phi_0, k_1, k_2, s_b, \tau\) | `1.410, 7.606, 3.975, 0.435, 0.030` | Dual-wavenumber phase profile |
| | \(B_0, B_{\mathrm{pk}}, d_0, d_1, s_d, w_d\) | `1.083, 2.334, 2.583, 5.831, 0.310, 0.048` | 2nd-harmonic amplitude & phase offset |
| **3-D tilt** | `tilt_deg` (\(\gamma_1, \gamma_2\)) | `[+17°, +17°]` | \(C_2\)-chiral beat-plane tilt |
| **3-D OOP** | `r_perp`, `r_perp2` | **`0.25`**, **`0.05`** | 1ω (`6.14e-17 N·m`) & 2ω (`1.23e-17 N·m`) ratios |
| **(1ω+2ω)** | `s_w`, `delta_perp`, `delta_perp2` | **`0.20`**, **`+π/2`**, **`+π/2`** | Proximal decay length & quadrature phase |
| | `perp_sign` (\(\sigma_1, \sigma_2\)) | **`[+1, -1]`** | Cis/trans chiral sign (enhances tilt spin) |
| **Helix** | `drive_amp` (\(\alpha_1, \alpha_2\)) | `[1, 1]` or **`[1, 0.95]`** | 5% physiological asymmetry -> left-handed helix |

---

## 4. Calibrated Benchmarks & Experimental Comparison

Summary tables and the 4-panel calibration figure live in `output/chlamy_oop/` (`oop_summary.csv`, `oop_step1_summary.csv`, `phase2_3d_summary.pdf/.png`).

| Case | Setup | `U_fwd` [µm/s] | `spin` [Hz] | Helix `R` [µm] | Pitch [µm] | Tip `z_pp` [µm] | Centerline 3-D `sv3/sv1` | `f_body_min` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Live cells** *(lit.)* | Wild-type water | ~80–128 (50 Hz) | **1–2 Hz** | ~1–5 (LH) | ~50–100 | 3-D closed loop | Proximal chiral torsion | > 1 |
| `OOP0` *(baseline)* | `tilt=[17,17]`, `r_perp=0` | 80.1 | -0.685 | 0.03 | -116.9 | 2.49 | 0.67% | 1.124 |
| `OOP10F` *(control)* | `tilt=[0,0]`, `r_perp=0.10` | 81.1 | -0.195 | 0.18 | -414.8 | 0.62 | 0.54% | 1.124 |
| **`SYM3D`** *(3-D, `10⁻⁴`)* | `drive_amp=[1, 1]` | **79.6** | **-1.166** | 0.00 | -68.2 | 3.54 | **1.88%** | **1.124** |
| **`HELIX3D`** *(3-D, `10⁻⁴`)* | `drive_amp=[1, 0.95]` | **77.6** | **-1.092** | **1.49** | **-69.9** | 3.41 | **1.82%** | **1.124** |

Key verification findings:

- **ODE tolerance convergence:** Comparing `tol = 1e-3` vs `1e-4` on `HELIX3D` yields `< 0.7%` difference across all metrics (`ΔU = 0.19%`, `Δspin = 0.09%`, `ΔR = 0.67%`).
- **Linear superposition:** Geometric tilt (`-0.685 Hz`) + pure dynamical OOP (`-0.195 Hz`) = `-0.880 Hz`, matching the combined run (`-0.875 Hz`) within `0.6%`. Reversing `perp_sign` to `[-1, +1]` cancels spin to `-0.497 Hz`; same-sign `[+1, +1]` gives zero spin gain (`-0.685 Hz`) and converts OOP torque into body pitching.
- **Role of the 2ω term (`r_perp2 = 0.05`):** Any pure-1ω 3-D periodic curve lies strictly inside a tilted plane; adding 5% 2ω expands genuine tip-orbit non-planarity (`sv3/sv1`) by `+46%` with negligible cost to thrust or clearance.

---

## 5. Rebuilding Geometry or Moment Tables (Advanced)

The tracked geometry CSVs and `chlamydomonas/data/fit_1p00A.mat` are ready to use out of the box. To regenerate them from raw experimental recordings under `chlamydomonas/Refers/`:

```matlab
% 1. Regenerate geometry CSVs (writes into geometry/)
addpath('geometry');
generate_chlamydomonas_geometry();

% 2. Invert measured gait -> tabulate -> compact dual-k analytic fit
addpath('chlamydomonas/core');
raw_mat = ['chlamydomonas/Refers/doi_10_5061_dryad_sbcc2frj6__v20250313/' ...
           'Quantification_of_Flagellar_Gait_Changes_Gutierrez_Data/Data/' ...
           'Newtonian 1p00 A_Viscosity1.00.mat'];
tab = chlamy_inverse_moment(raw_mat);
[fit, P] = chlamy_fit_moment(tab);
save('chlamydomonas/data/fit_1p00A.mat', 'fit', 'P');
```

---

## 6. Outputs & Analysis Notes

Each preset writes standard solver outputs to `output/<case_name>/` (`simulation_results.mat`, `timeseries.csv`, `config.json`, `simulation.log`, `progress.txt`, optional `bacteria_3d_motion.mp4`) plus `swim_analysis.mat` from `chlamy_swim_analysis.m`.

- **Multi-turn accuracy:** `chlamy_swim_analysis.m` accumulates frame-to-frame \(\mathrm{SO}(3)\) rotation logarithms \(\log(R_m^\top R_{m+1})\) and midpoint-projected body displacements, so `spin_Hz`, `total_spin_turns`, `V_body_um_s`, and `helix` parameters remain exact over multi-second, multi-revolution runs (such as `HELIX3D_2S`).
- **Clearance post-check:** The solver has no steric contact model; always confirm `s.f_body_min > 1` (no flagellum-body penetration) and `s.d_flag_min_um > 0.25` (no flagellum-flagellum overlap).
