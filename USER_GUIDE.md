# User Guide — Coupled Body-Pili Simulation

> **Author:** Lucien (lucien-6@qq.com)  
> **Version:** 5.5.0
> **Updated:** 2026-04-18

This guide is a short entry point. **Full theory, parameters, and file map:** see `README.md`.

---

## 1. What you need

- **MATLAB** (R2020b or later recommended). No Python environment is required to run the simulator (`requirements.txt` documents this).
- **MATLAB Add-On / File Exchange — S2 Sampling Toolbox** (*Suite of functions to perform uniform sampling of a sphere*, Anton Semechko): required for **capsule** body-surface generation (`capsule_sphere_mesh_cache.m`, `capsule_sphere_cap_points.m`) — e.g. `ParticleSampleSphere`, `SpiralSampleSphere`, `IcosahedronMesh`, `QuadCubeMesh`, `SubdivideSphericalMesh`. Install from Add-On Explorer so the functions are on the path.

---

## 2. How to run

### Option A — Inline parameters (default)

1. Open `main_coupled_simulation.m`.
2. Edit the `cfg` block (case name, geometry, pili, time span, etc.). The current inline template is a one-pilus high-frequency active-wave example (`OnePili-Active-HighFre-Test`).
3. Run the script.

### Option B — JSON configuration

1. Copy **`config_example.json`** (full template aligned with `load_sim_config.m` / `save_sim_config.m`) or use a file under **`configs/`** (e.g. `OnePili-P1-L3.json`).
2. In `main_coupled_simulation.m`, set `config_file` to your JSON path.
3. Run.

JSON must include the same **top-level sections** the loader expects: `case_name`, `body`, `pili`, `material`, `fluid`, `simulation`, `brownian`, `output`, `flowfield`, and **`background_flow`** (v5.0; use `"type": "off"` to disable). Optional per-pilus fields: **`intrinsic_curvature`**, **`active_moment`**. Optional **`body.eps_body`** for RSM regularisation override. See **`README.md`** §2.2 for the cfg ↔ JSON mapping table.

### Option C — Programmatic / batch

```matlab
cfg = load_sim_config('config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

Use `batch sweep/` scripts or `batch_run_from_json.m` for unattended sweeps. The folder includes **eight** sweep drivers (e.g. `batch_sweep_alpha`, `batch_sweep_eps_body_RSM`, `batch_sweep_N_body_eps_RSM`, …) plus **`batch_run_from_json`** and **`batch_repeat_brownian`**. For **multiple independent runs of the same Brownian-on case** (outputs under `output/<case_name_base>/001/`, `002`, …), run `batch sweep/batch_repeat_brownian.m` (requires `enable_brown = true`; preview off by default). Set **`brown_rng_seed_mode`**: `'default'` (omit per-run seed; same **MMddHHmmss** time default as a single run) or `'fixed'` (**`cfg.brown_rng_seed = rng_seed_base + k`**, reproducible batch).

**v5.1.1:** After each run, MATLAB **closes** the summary and flow-field figure windows once `summary_plot.png` and `flowfield_panels.png` are saved, so long batches do not leave hundreds of figures open.

**v5.2.0:** **Brownian RNG** — default twister seed is **local wall-clock `MMddHHmmss`** as an integer (`brown_rng_seed_from_clock.m`); set **`cfg.brown_rng_seed`** or JSON **`brownian.rng_seed`** to fix a seed. When Brownian is on, the realised seed is saved in **`config.json`**, **`brownian_rng_seed.txt`**, **`simulation.log`**, and **`results.brown_rng_seed`** / **`params.brown_rng_seed`** in `simulation_results.mat`.

**v5.3.0:** **Body RSM default** — if **`body.eps_body`** is omitted, the simulator uses **`eps_body = h_mean / 2.6667`** (see `README.md` §8.5). **`batch_sweep_N_body_eps_RSM.m`** performs a nested **`N_body` × `k`** sweep (`eps_body = k·h_mean`) vs Stokes settling and writes CSV/MAT summaries under `batch sweep/`.

**v5.3.1 / v5.4.4:** **MP4 animation** — the exporter now probes the **actual physical pixel size** returned by **`getframe`**, clamps the encoded frame size to a codec-safe **3840×2160** even-valued box, and resizes frames on write if Windows HiDPI scaling still returns oversized captures. This avoids **`VideoWriter`** failures caused by frame dimensions that remain too large after logical figure resizing.

---

## 3. Main outputs

Under `output/<case_name>/` (or timestamped folder):

| File | Role |
|------|------|
| `simulation_results.mat` | Full state and results |
| `timeseries.csv` | Scalar time series for external tools |
| `config.json` | Snapshot of the run configuration |
| `summary_plot.png` | Overview figure (600 dpi) |
| `flowfield_panels.png` | 2×2 flow-field panels (600 dpi) |
| `bacteria_3d_motion.mp4` | Optional animation |
| `simulation.log` | Text log |
| `brownian_rng_seed.txt` | v5.2.0: Brownian twister seed + one-line reproduce hint (only if `enable_brown`) |

---

## 4. Post-processing batch results (v5.1)

After you have one or more completed runs (each Case folder must contain `timeseries.csv`, and ideally `config.json`), aggregate metrics into one Excel workbook:

1. In MATLAB, `cd` to the project root and run `addpath('post_processing')`.
2. Run `postprocess_cases()`.
3. Choose the **parent directory** that contains your Case subfolders (e.g. `output/`). In the dialog, multi-select the cases to analyse.
4. Pick a save path for the `.xlsx` summary file.

The script computes per-case summaries (e.g. sinking velocity statistics, maximum horizontal displacement, 3-D tortuosity, MSD-based translational diffusion fit) and writes a **Summary** sheet. On **Windows**, Excel COM may apply table formatting if Excel is available; you can decline auto-open to avoid file-lock prompts. See `CHANGELOG.md` for details.

---

## 5. Flow-field figure (`flowfield_panels.png`)

- **Layout:** 2×2 — 3D XZ / 3D YZ / 2D XZ / 2D YZ.
- **2D panels:** Velocity magnitude on a `jet` scale + **normalized** direction arrows (equal length per arrow; magnitude is read from the color field).
- **v3.3.1:** 2D arrows are single dark-gray vectors (no white halo). Arrows are drawn **50% larger** than the previous dark layer for readability.
- **v3.3.2:** Batch sweep scripts under `batch sweep/` were refreshed (new `N_body` / `seg_per_um` drivers, removed `batch_sweep_pili_count.m`); see `CHANGELOG.md`.
- **v3.3.3:** `batch_sweep_alpha.m` sweeps **`alpha_0` for a body-only capsule** (empty `Pili_Matrix`, no pili).
- **v3.3.4:** Removed `batch_sweep_pili_site.m`; use **`batch_sweep_pili_length.m`** with **`site_list`** for attachment-site studies.
- **v3.3.5:** `batch_sweep_major_axis.m` sweeps **semi-major axis** from **0.4 to 1.6 µm** (step **0.2 µm**), body-only (no pili).
- **v4.0.0:** Full-system effective friction Brownian dynamics — diffusion tensor dynamically computed per step from the coupled body-pili linear system via probing method; pili drag naturally suppresses thermal diffusion; optional Fixman midpoint correction (`fixman_midpoint` in JSON config); new `build_coupled_lin_sys.m` and `compute_effective_body_brownian.m`; deleted old body-only `compute_body_brownian_params.m`; fixed missing `fixman_midpoint` in inline config and `OnePili-P1-L3.json`.
- **v4.1.0:** `batch_repeat_brownian.m` — batch **independent repeats** for the same Brownian-on configuration; results in `output/<case_name_base>/001/`, `002`, …; default `preview_initial = false`, `auto_confirm = true`; reproducible RNG seeds per repeat.
- **v5.0.0:** Three new physics features — **intrinsic curvature** (per-pili static, `frozen`/`zero`/`constant`/`helix`), **active moments** (per-pili time-varying, `off`/`traveling_wave`/`custom`), and **background flow** (global steady, `off`/`uniform`/`shear_z`). New files: `coupled_background_flow.m`, `coupled_active_moment.m`. All default to off — fully backward-compatible.
- **v5.1.0:** **`post_processing/postprocess_cases.m`** — batch post-processing of multiple Case folders → Excel summary (metrics from `timeseries.csv` + metadata from `config.json`). English UI.
- **v5.1.1:** Summary and flow-field figures are **`close`d** after PNG export (`run_simulation.m`, `compute_plot_flowfield.m`) to keep batch runs from opening many windows.
- **v5.2.0:** **Brownian twister seed** — default **`MMddHHmmss`** integer; optional **`cfg.brown_rng_seed`** / JSON **`brownian.rng_seed`**; seed persisted in `config.json`, `brownian_rng_seed.txt`, log, and `simulation_results.mat`; `batch_repeat_brownian` uses **`cfg.brown_rng_seed = rng_seed_base + k`**.
- **v5.3.2:** **Node director interpolation** — new shared utility `compute_node_directors.m` (segment directors → node directors via SVD polar decomposition with `det < 0` reflection fix). Fixes a correctness bug in `build_coupled_lin_sys.m` (missing `det < 0` check, single-segment index out-of-bounds). `ode_coupled_rhs.m` elastic/active moment writes vectorised for better performance. `validate_simulation_params.m` adds warnings for negative `wavenumber` and non-zero `d3_amp` in traveling-wave active moments.

- **v5.4.3:** **Header comment cleanup** — all 47 `.m` files: removed `Updated:` history lines; headers now end with `Author:` / `Created:` only. `build_coupled_lin_sys.m` Syntax comment fixed (3-output form). `brown_rng_seed_from_clock.m` and `postprocess_cases.m` header comments standardised.

- **v5.4.4:** **HiDPI animation export hardening** — `run_simulation.m` now clamps the encoded MP4 frame size to a legal even-valued MPEG-4 box after probing `getframe` output, then rescales oversized frames before `writeVideo` if needed. The inline example in `main_coupled_simulation.m` now defaults to a single-pilus high-frequency active case.

**v5.5.0:** **Capsule sphere caps (S2 toolbox)** — hemispherical cap Stokeslet points use Anton Semechko’s toolbox: optional snap of the cap count to Ico/Quad native vertex totals (±5 from area-weighted allocation), else `ParticleSampleSphere` below 1000 points (see `capsule_sphere_cap_points.m` for documented `Etol`/`Dtol`/`Nitr`), else `SpiralSampleSphere`. Equator points with **exact** `x==0` are duplicated onto both caps. Native mesh sizes come from `SubdivideSphericalMesh` via `capsule_sphere_mesh_cache.m` (not a closed-form formula).

Tunable via `cfg`: `flowfield_grid_spacing`, `flowfield_grid_range`, `flowfield_quiver_skip`, `flowfield_chunk_size`.

---

## 6. New in v5.0: Intrinsic Curvature, Active Moments, Background Flow

All three features are configured in `main_coupled_simulation.m` (inline) or in JSON config. **`config_example.json`** includes a `background_flow` block and per-pilus `intrinsic_curvature` / `active_moment` defaults.

### Background flow

```matlab
cfg.background_flow.type           = 'off';       % 'off' | 'uniform' | 'shear_z'
% cfg.background_flow.U_max        = 1e-6;        % [m/s]
% cfg.background_flow.direction    = [1, 0, 0];
% cfg.background_flow.shear_axis   = [0, 0, 1];
% cfg.background_flow.shear_height_m = 3e-6;      % [m]
```

Equivalent JSON (minimal, flow off):

```json
"background_flow": {
    "type": "off",
    "U_max": 0,
    "direction": [1, 0, 0],
    "shear_axis": [0, 0, 1],
    "shear_height_m": 0
}
```

### Per-pili intrinsic curvature

```matlab
pili_config(1).intrinsic_curvature = 'frozen';     % Keep initial geometry (default)
% pili_config(2).intrinsic_curvature = struct('type','constant','kappa',[0,0.5,0]);
```

In JSON, under each `pili[]` entry: `"intrinsic_curvature": "frozen"` or an object with `"type": "helix"`, `kappa_d2`, `kappa_d3`, etc. (see `init_coupled_params`).

### Per-pili active moment

```matlab
pili_config(1).active_moment = struct('type', 'off');     % Default
% pili_config(2).active_moment = struct('type','traveling_wave', ...
%     'amplitude_Nm', 1e-18, 'omega_Hz', 10, 'wavenumber', 2*pi, ...
%     'd1_amp', 1, 'd2_amp', 0, 'd3_amp', 0);
```

In JSON, use `"active_moment": { "type": "off" }` or a `traveling_wave` object with **`amplitude_Nm`**, **`omega_Hz`**, **`wavenumber`**, optional **`phase_rad`** (time-phase offset [rad], default 0), **`d1_amp`** / **`d2_amp`** / **`d3_amp`** (same names as MATLAB).

Physical units are automatically converted to dimensionless at initialisation.

---

## 7. Version and changes

- Current version is printed at the start of `simulation.log` and listed at the top of `README.md` / this file.
- See **`CHANGELOG.md`** for release notes.

---

## 8. License

MIT License — see `README.md` (license section).
