# User Guide — Coupled Body-Pili Simulation

> **Author:** Lucien (lucien-6@qq.com)
> **Version:** 7.3.1
> **Updated:** 2026-05-04

This guide is a short entry point. **Full theory, parameters, and file map:** see `README.md`.

---

## 1. What you need

- **MATLAB** (R2020b or later recommended). No Python environment is required to run the simulator (`requirements.txt` documents this). **Interactive HTML** figure export (`exportgraphics` → `.html`) requires **MATLAB R2026a or later**; PNG and PDF exports work from **R2020a** onward (`exportgraphics`).
- **MATLAB Add-On / File Exchange — S2 Sampling Toolbox** (*Suite of functions to perform uniform sampling of a sphere*, Anton Semechko): required for **capsule** body-surface generation (`capsule_sphere_mesh_cache.m`, `capsule_sphere_cap_points.m`) — e.g. `ParticleSampleSphere`, `SpiralSampleSphere`, `IcosahedronMesh`, `QuadCubeMesh`, `SubdivideSphericalMesh`. Install from Add-On Explorer so the functions are on the path.

---

## 2. How to run

### Option A — Inline parameters (default)

1. Open `main_coupled_simulation.m`.
2. Edit the `cfg` block (case name, geometry, pili, time span, etc.). The current inline template is a **two-pilus** symmetric high-frequency traveling-wave demo (`TwoPili-Active-HighFre-Test-Symmetry` in `main_coupled_simulation.m`).
3. Run the script.

### Option B — JSON configuration

1. Copy **`config_example.json`** (full template aligned with `load_sim_config.m` / `save_sim_config.m`) to a path of your choice. A common layout is a local **`configs/`** folder (the repository **`.gitignore` ignores `configs/`** so machine-specific paths are not committed); then set `config_file` in `main_coupled_simulation.m` (e.g. `'configs/my_case.json'`).
2. In `main_coupled_simulation.m`, set `config_file` to your JSON path.
3. Run.

JSON must include the same **top-level sections** the loader expects: `case_name`, `body`, `pili`, `material`, `fluid`, `simulation`, `brownian`, `output`, `flowfield`, and **`background_flow`** (v5.0; use `"type": "off"` to disable). Optional per-pilus fields: **`intrinsic_curvature`**, **`active_moment`**. Optional **`body.eps_body`** for RSM regularisation override. See **`README.md`** §2.2 for the cfg ↔ JSON mapping table.

### Option D — Custom geometry mode (v6.0, extended in v7.1)

Import arbitrary 3D point clouds for both the body and pili:

1. Copy **`config_custom_example.json`** as a template.
2. Set `"geometry_mode": "custom"` at the top level.
3. Pre-shipped example files live under `geometry/`:
   - `example_body_sphere.csv`, `example_body_capsule.csv`
   - `example_pilus_linear.csv`, `example_pilus_helix.csv`
   - regenerate with `geometry/generate_example_geometry.m`.
4. Under `"custom_geometry"`, specify:
   - `"body_points_file"`: body point cloud (`.csv`, `.txt`, or `.mat`; N×3 or 3×N).
   - `"attachment_sites"`: array of `{ "position": [x,y,z], "direction": [dx,dy,dz], "label": "name" }`.
   - `"pili"`: array of `{ "points_file": "...", "site_index": 1, "deflection_z": 0, "deflection_x": 0, "intrinsic_curvature": ..., "active_moment": {...}, "root_coupling": "...", "motor_torque": {...} }`.
   - (optional) `"initial_orientation"`: exactly one of
     `{ "quaternion": [w,x,y,z] }`, `{ "axis_angle": {"axis":[x,y,z], "angle_rad": t} }`, or `{ "euler_zyx": [yaw,pitch,roll] }`.
   - (optional) `"L_0_suggest_mode"`: `"volume"` (default; equivalent-sphere radius from triangulated volume) or `"extent"` (half the largest bounding-box extent).
5. Set `config_file = 'config_custom_example.json'` in `main_coupled_simulation.m` and run. Alternatively leave `config_file` empty and use the inline custom block inside `main_coupled_simulation.m` — both paths reference the same pre-shipped example files.

**Key constraints and semantics:**
- Custom body and custom pili **must** be used together — mixing with parametric models is not allowed.
- `alpha_0` and `eps_body` for custom mode live under the `custom_geometry` node. `alpha_0` inherits the parametric convention (rotation about +Y assuming body-X is the long axis); for arbitrary custom bodies, prefer `initial_orientation`, which **overrides** `alpha_0` when present.
- **Pili segment uniformity (v7.1):** every user-supplied pilus centerline file **must** be pre-discretised with equal-length segments (relative spread ≤ 1e-3). Non-uniform inputs are rejected.
- **Pili deflection (v7.1):** in custom mode, `deflection_z` / `deflection_x` are now applied in the pilus **local basis** (built via Gram-Schmidt around the aligned root tangent), instead of the lab frame after alignment. The local basis uses `e3` = aligned tangent and `e1` = bending-plane normal, so `deflection_z` rotates about `e3` (a TWIST around the tangent) and `deflection_x` rotates about `e1` (a TILT in the bending plane). **Naming caveat:** in parametric mode the pilus is built along `+X`, so parametric `deflection_x` is a TWIST (rotation about the tangent) and parametric `deflection_z` is a TILT — the geometric roles are **swapped** between the two modes. The `RX` matrix keeps the same left-handed sign convention in both modes; when porting a non-zero deflection from parametric to custom (or vice versa), swap `deflection_x` ↔ `deflection_z`.
- **`root_dir` semantic (v7.1):** stores the **actual** first-segment tangent after alignment + deflection (previously stored the raw `site.direction`).
- **Alignment diagnostic (v7.1):** the rotation angle used to align the centerline's first segment to `site.direction` is logged per pilus; if it exceeds 45°, a warning is printed so users can verify the pilus file orientation.
- **Per-pilus extended physics (v7.1):** `intrinsic_curvature`, `active_moment`, `root_coupling`, and `motor_torque` now flow through the custom JSON/inline path identically to parametric mode.
- **Custom-mode twist (psi0):** the initial twist angle for custom pili is `psi0 = 0` because a generic 3D centerline does not uniquely define a natural binormal. Use `intrinsic_curvature = "constant"` / `"helix"` to request a non-zero intrinsic twist.
- **Duplicate-point handling (v7.1):** both body and pili centerline readers automatically remove consecutive / near-duplicate points (tol = 1e-6 of bounding-box extent) and error out if fewer than the minimum required distinct points remain.
- The body surface is automatically triangulated (`alphaShape` with a logarithmic alpha scan, `convhulln` fallback) for `trisurf` visualisation and divergence-theorem volume computation. A non-closed-surface diagnostic (Euler characteristic check) warns when the triangulation is non-manifold.
- **Flow field in custom mode (v7.1):** `compute_plot_flowfield.m` now builds an alphaShape-based interior mask and renders the body via `trisurf` using the stored triangulation — the previous `major_axis` / `minor_axis` hard dependency has been removed.
- CSV geometry files support comma-separated, space-separated, and tab-separated delimiters, with comment lines starting with `#` or `%`. N×2 inputs are rejected by default; pass `AllowPlanar2D = true` to the reader to permit Z=0 padding.

### Option C — Programmatic / batch

```matlab
cfg = load_sim_config('config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

Use `batch sweep/` scripts or `batch_run_from_json.m` for unattended sweeps. The folder includes **eight** sweep drivers (e.g. `batch_sweep_alpha`, `batch_sweep_eps_body_RSM`, `batch_sweep_N_body_eps_RSM`, …) plus **`batch_run_from_json`** and **`batch_repeat_brownian`**. For **multiple independent runs of the same Brownian-on case** (outputs under `output/<case_name_base>/001/`, `002`, …), run `batch sweep/batch_repeat_brownian.m` (requires `enable_brown = true`; preview off by default). Set **`brown_rng_seed_mode`**: `'timestamp'` (omit per-run seed; same **MMddHHmmss** time default as a single run; `'default'` kept as a silent backward-compatible alias) or `'fixed'` (**`cfg.brown_rng_seed = rng_seed_base + k`**, reproducible batch).

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
| `summary_plot.pdf` | Same overview, vector PDF (`exportgraphics`) |
| `summary_plot.html` | Same overview, interactive HTML (MATLAB **R2026a+** only, if written) |
| `flowfield_panels.png` | 2×2 flow-field panels (600 dpi) |
| `flowfield_panels.pdf` / `.html` | Flow-field panels, vector / interactive (same rules as above) |
| `initial_pose.pdf` / `.html` | Initial pose preview, vector / interactive (when preview is on) |
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

The script computes per-case summaries (e.g. sinking velocity statistics, maximum horizontal displacement, 3-D tortuosity, MSD-based translational diffusion fit) and writes a **Summary** sheet. It also saves a **batch summary figure** next to the workbook: same base name as **`.png`**, **`.pdf`**, and **`.html`** (interactive on R2026a+). On **Windows**, Excel COM may apply table formatting if Excel is available; you can decline auto-open to avoid file-lock prompts. See `CHANGELOG.md` for details.

---

## 5. Flow-field figure (`flowfield_panels.png`)

- **Layout:** 2×2 — 3D XZ / 3D YZ / 2D XZ / 2D YZ.
- **2D panels:** Velocity magnitude on a `jet` scale + **normalized** direction arrows (equal length per arrow; magnitude is read from the color field).
- **v3.3.1:** 2D arrows are single dark-gray vectors (no white halo). Arrows are drawn **50% larger** than the previous dark layer for readability.
- **v3.3.2:** Batch sweep scripts under `batch sweep/` were refreshed (new `N_body` / `seg_per_um` drivers, removed `batch_sweep_pili_count.m`); see `CHANGELOG.md`.
- **v3.3.3:** `batch_sweep_alpha.m` sweeps **`alpha_0` for a body-only capsule** (empty `Pili_Matrix`, no pili).
- **v3.3.4:** Removed `batch_sweep_pili_site.m`; use **`batch_sweep_pili_length.m`** with **`site_list`** for attachment-site studies.
- **v3.3.5:** `batch_sweep_major_axis.m` sweeps **semi-major axis** from **0.4 to 1.6 µm** (step **0.2 µm**), body-only (no pili).
- **v4.0.0:** Full-system effective friction Brownian dynamics — diffusion tensor dynamically computed per step from the coupled body-pili linear system via probing method; pili drag naturally suppresses thermal diffusion; optional Fixman midpoint correction (`fixman_midpoint` in JSON config); new `build_coupled_lin_sys.m` and `compute_effective_body_brownian.m`; deleted old body-only `compute_body_brownian_params.m`; fixed missing `fixman_midpoint` in inline config and sample JSON configs (see `CHANGELOG` [4.0.0]).
- **v4.1.0:** `batch_repeat_brownian.m` — batch **independent repeats** for the same Brownian-on configuration; results in `output/<case_name_base>/001/`, `002`, …; default `preview_initial = false`, `auto_confirm = true`; reproducible RNG seeds per repeat.
- **v5.0.0:** Three new physics features — **intrinsic curvature** (per-pili static, `frozen`/`zero`/`constant`/`helix`), **active moments** (per-pili time-varying, `off`/`traveling_wave`/`custom`), and **background flow** (global steady, `off`/`uniform`/`shear_z`). New files: `coupled_background_flow.m`, `coupled_active_moment.m`. All default to off — fully backward-compatible.
- **v5.1.0:** **`post_processing/postprocess_cases.m`** — batch post-processing of multiple Case folders → Excel summary (metrics from `timeseries.csv` + metadata from `config.json`). English UI.
- **v5.1.1:** Summary and flow-field figures are **`close`d** after PNG export (`run_simulation.m`, `compute_plot_flowfield.m`) to keep batch runs from opening many windows.
- **v5.2.0:** **Brownian twister seed** — default **`MMddHHmmss`** integer; optional **`cfg.brown_rng_seed`** / JSON **`brownian.rng_seed`**; seed persisted in `config.json`, `brownian_rng_seed.txt`, log, and `simulation_results.mat`; `batch_repeat_brownian` uses **`cfg.brown_rng_seed = rng_seed_base + k`**.
- **v5.3.2:** **Node director interpolation** — new shared utility `compute_node_directors.m` (segment directors → node directors via SVD polar decomposition with `det < 0` reflection fix). Fixes a correctness bug in `build_coupled_lin_sys.m` (missing `det < 0` check, single-segment index out-of-bounds). `ode_coupled_rhs.m` elastic/active moment writes vectorised for better performance. `validate_simulation_params.m` adds warnings for negative `wavenumber` and non-zero `d3_amp` in traveling-wave active moments.

- **v5.4.3:** **Header comment cleanup** — all 47 `.m` files: removed `Updated:` history lines; headers now end with `Author:` / `Created:` only. `build_coupled_lin_sys.m` Syntax comment fixed (3-output form). `brown_rng_seed_from_clock.m` and `postprocess_cases.m` header comments standardised.

- **v5.4.4:** **HiDPI animation export hardening** — `run_simulation.m` now clamps the encoded MP4 frame size to a legal even-valued MPEG-4 box after probing `getframe` output, then rescales oversized frames before `writeVideo` if needed. **`main_coupled_simulation.m` inline defaults evolve over releases** — check `cfg.case_name` in the file (currently **`TwoPili-Active-HighFre-Test-Symmetry`**).

**v5.5.0:** **Capsule sphere caps (S2 toolbox)** — hemispherical cap Stokeslet points use Anton Semechko’s toolbox: optional snap of the cap count to Ico/Quad native vertex totals (±5 from area-weighted allocation), else `ParticleSampleSphere` below 1000 points, else `SpiralSampleSphere`. Equator points with **exact** `x==0` are duplicated onto both caps. Native mesh sizes come from `SubdivideSphericalMesh` via `capsule_sphere_mesh_cache.m` (not a closed-form formula).

**v5.5.1:** **Sphere-cap Particle sampling — deterministic init + disk cache** — `ParticleSampleSphere` now uses toolbox default tolerances (`Etol=1e-5`, `Dtol=1e-4 deg`, `Nitr=1e4`; no custom overrides). Initial particle positions are supplied as `SpiralSampleSphere(N, false)` — fully deterministic, no RNG dependency. For `100 < SPN < 1000`, results are cached to `cache/sphere_particle/spherical_points_NXXXX.mat` via the new `sphere_particle_sample_cached.m`; a cache hit skips the iterative Riesz-energy optimisation entirely. Each file stores a `meta` struct (tolerances, `Vo_method`, format version); any mismatch triggers recomputation and overwrites the file. I/O failures fall back silently to direct computation.

**v5.5.2:** **Cylinder circumferential grid — TN multiple of 4** — In `capsule_body_building.m`, after the area-based estimate `TN = ceil(CPN / Turns)`, the code sets **`TN = max(round(TN / 4) * 4, 4)`** before building the azimuthal grid. Total cylinder points and actual `N_body` may change relative to the previous rule.

**v5.5.3:** **Default ODE timeout `tlim`** — Project-wide default is now **3600 s** (was **360 s**): `load_sim_config` JSON fallback, inline `main_coupled_simulation.m`, `config_example.json`, batch sweep scripts, and docs.

**v5.5.4:** **Hot-path vectorisation (performance only, bit-exact)** — `build_coupled_lin_sys.m` (7 sites) and `build_global_resistance.m` (3 sites) now replace per-index `for` loops with strided indexing, triangular masks, block assignments, and a single `blkdiag(...)` call. Reference `ProfileShort` run: **end-to-end `run_simulation` −14.8 %**, `build_coupled_lin_sys` −17.1 %, `ode_coupled_rhs` −16.4 %. New precomputed `params.quad_weights` caches global trapezoidal quadrature weights once (replaces the old per-step `get_trapezoidal_weight` helper). Numerical equivalence was verified by a dedicated regression harness (tolerance `1e-14`; observed max rel. Frobenius error `~1.8 × 10⁻²⁰`, i.e. effectively bit-exact) before the harness was archived. Frozen pre-optimisation copies are retained under `_legacy/` (`build_coupled_lin_sys_v0.m` / `build_global_resistance_v0.m`) as historical reference. **No configuration, API, or physics change** — existing inline cfg and JSON configs reproduce identical trajectories.

**v5.5.5:** **Direction 1 — `M_weighted` LU factorisation reuse (performance only, bit-exact on all D1-owned quantities)** — `build_coupled_lin_sys.m` now builds a single reusable LU via `decomposition(M_weighted, 'lu', 'CheckCondition', false)` and propagates it downstream through `hydro.dM_weighted`. `ode_coupled_rhs.m::apply_Ainv_u` consumes this factorisation for its background-flow solve, eliminating one full `O(N³)` LU factorisation per RHS call on the `bg ≠ off` path. Notable context: `M_weighted` mixes RSM (body) and RSS (pili) kernels and is **non-symmetric** (empirically `≈ 0.56` relative Frobenius asymmetry) and **not positive definite** — so the factorisation must be LU, not Cholesky. `CheckCondition=false` disables MATLAB's automatic `rcond` estimate on each solve (saves `O(N²)` per invocation). Direct A/B wall-clock benchmark at DOFs = 1807: `bg = off` median **767.5 → 731.6 ms (−4.7 %)**; `bg = uniform` median **868.2 → 731.7 ms (−15.7 %)**. All D1-owned quantities are bit-exact against the frozen `_legacy/ode_coupled_rhs_v0.m` / `_legacy/build_coupled_lin_sys_v0.m` reference (`rel_err = 0`); the only non-zero check is the end-to-end `dZ` at a perturbed state (`~4.2 × 10⁻¹¹`, below the condition-number-amplified floor of `~1 × 10⁻⁹` — residual demonstrably inherited from the 5.5.4 `B_global` vectorisation round, not introduced by Direction 1). One-shot developer harnesses (`profile_short_run.m`, `compare_profiles.m`, `test_optimization_invariance.m`, `benchmark_direction1.m`) were removed from the live tree after validation; the `_legacy/_v0` references are kept. **No configuration, API, or physics change.**

**v5.5.6:** **Brownian module polish (zero numerical impact) + config example tweaks** — Brownian: `validate_simulation_params.m` gains two new warnings — (a) `enable_brown + pili present + fixman_midpoint=false` signals that the `kB·T·∇·M_eff` Ito drift is silently dropped (simple Euler-Maruyama) and recommends enabling Fixman midpoint for long-time diffusion statistics; (b) `dt_B > 10% × min(τ_trans, τ_rot)` via Stokes-Einstein spherical back-of-envelope suggests a smaller `dt_B`. Both warnings fire only when `enable_brown = true`, so Brownian-off workflows are unaffected. `compute_effective_body_brownian.m` adds **read-only diagnostic hooks** — pre-symmetrisation relative Frobenius asymmetry `> 1e-6` raises `AsymmetryHigh`; post-dimensionalisation `rcond(M_eff_dim) < 1e-10` raises `IllConditioned`. **Step 1 → Step 5 numerics and the `[L_eff, D_eff, M_eff_dim]` output signature are unchanged.** `batch sweep/batch_repeat_brownian.m`: seed mode `'default'` → `'timestamp'` (self-documenting; `'default'` kept as a silent backward-compatible alias, normalised at runtime). Documentation: `Brownian_Implementation_Plan.md` adds §1.5 (Scope & Simplifying Assumptions, clarifying body-only Brownian DOFs and the omitted pili internal thermal fluctuations) and §2.4.1 (Stochastic Interpretation & Missing Drift, with the exact Itô SDE, drift-omission magnitude estimate, Fixman recovery mechanism, and a practical decision matrix). Config example: `config_example.json` updates pili Young modulus `E = 1.768e7 → 3.279e8` (stiffer demo pili) and `T_obs = 30.0 → 0.5` s (shorter demo run). **No API change; no numerical path change.**

**v5.5.8:** **Publication exports + optional interactive HTML (zero numerical impact)** — New **`export_figure_png_pdf_html.m`** saves each exported figure as **PNG** (600 dpi), **PDF** (vector via `exportgraphics`), and on **MATLAB R2026a+** an **interactive `.html`** web canvas (pan/zoom/rotate per MathWorks docs). Applies to `initial_pose`, `summary_plot`, `flowfield_panels`, and the post-process Excel **same-stem** summary figure. **`main_coupled_simulation.m`** and **`run_simulation.m`** document a path policy: project root `addpath` at the main script plus a defensive `addpath` at the start of `run_simulation`; **`post_processing/`** and **`batch sweep/`** keep local `addpath` for standalone runs. Final *Files saved* / post-process completion text lists **`.html` only if the file exists** (`isfile`). Output-dir overwrite detection whitelists **`.pdf`** and **`.html`**. **`postprocess_cases.m`** now also writes the batch summary **`.png`/`.pdf`/optional `.html`** next to the workbook.

**v6.0.0:** **Custom geometry model import** — New `geometry_mode = 'custom'` allows importing arbitrary 3D point clouds for both the bacterial body and its pili, replacing the parameterized capsule/analytic-pili models. Five new functions: `read_geometry_csv_mat.m` (point cloud file reader, `.csv`/`.txt`/`.mat`, auto-detects orientation), `rotation_between_vectors.m` (Rodrigues' rotation formula for vector alignment), `custom_body_building.m` (body point cloud processing: `alphaShape` + `convhulln` auto-triangulation, volume calculation, attachment sites), `custom_pili_building.m` (pili centerline processing: arc length, PCHIP resampling, Rodrigues direction alignment, additional rotations), `validate_custom_geometry.m` (custom geometry configuration validation). Six core functions modified: `init_coupled_params.m` (geometry mode dispatch), `load_sim_config.m`/`save_sim_config.m` (JSON custom node parsing/persistence), `validate_simulation_params.m` (mutual exclusion validation + custom dispatch), `run_simulation.m` (`trisurf` visualization + conditional logging), `main_coupled_simulation.m` (inline custom example). New `config_custom_example.json` template. **Custom body and custom pili must be used together** (no mixing with parametric models; validated at startup). **Parametric mode is fully backward-compatible** — `geometry_mode` defaults to `'parametric'`, existing configurations require no changes.

**v5.5.7:** **Logging, reporting & persistence polish (zero numerical impact)** — `run_simulation.m` enriches the `[PARAMS]` section with per-pilus `intrinsic_curvature.type` / `active_moment.type` (with `traveling_wave` additionally printing `amplitude_Nm`, `frequency_Hz`, `wavelength`) and a global `Background flow` block covering the `off` / `uniform` / `shear_z` branches (with `U_max`, `direction`, `shear_axis`, `shear_height_m` as applicable). The final *Files saved* manifest now explicitly lists `Config : config.json` and — when `enable_brown = true` — `Seed : brownian_rng_seed.txt`, so the summary matches what is actually written to disk. The output-directory overwrite probe's file-type whitelist is extended with `.json` and `.txt`, and the empty-`case_name` branch in the parameters banner now prints the auto-generated `run-<MMddHHmmss>` identifier instead of the placeholder text. `save_sim_config.m` / `load_sim_config.m` now round-trip `output.compute_flowfield` (default `true`) so that reloading a saved config no longer silently re-enables flow-field computation. **All eight sweep drivers under `batch sweep/`** (`batch_sweep_alpha`, `batch_sweep_N_body`, `batch_sweep_major_axis`, `batch_sweep_pili_E`, `batch_sweep_pili_length`, `batch_sweep_seg_per_um`, `batch_sweep_eps_body_RSM`, `batch_sweep_N_body_eps_RSM`) now wrap each `run_simulation` call in `try/catch`, so a single failure no longer aborts the whole sweep, and each produces a persistent `batch sweep/<script>_summary_<YYYYMMDD_HHmmss>.log` with per-run status (`OK` / `FAILED` / `ABORTED`), wall-clock time, and the exception message on failure, plus a footer with elapsed time, success/failure counts, and CSV/MAT paths. The two RSM scripts additionally save `run_status`, `run_walltime`, and `run_errmsg` arrays inside their `.mat` outputs. A stray `d` character at the start of `batch_sweep_N_body_eps_RSM.m` (which previously made the script unparseable) is removed; `post_processing/postprocess_cases.m` replaces two deprecated `datestr(now, …)` calls with `datetime('now', 'Format', …)`. **No physics, ODE right-hand side, linear-system assembly, Brownian kernel, or time-stepping path is touched**, so any pre-5.5.7 case reproduces bit-exact `results.mat` and `timeseries.csv`; the new logs and summary lines are strictly additive.

**v6.1.0:** **Dual-channel logger + default parameter unification** — `run_simulation.m` introduces two logger closures: `logf` (console + file dual-write, file entries carry `[HH:mm:ss]` timestamps; console remains timestamp-free) and `logf_file` (file-only, for high-frequency diagnostics). Logger handles are injected into `params` via `init_coupled_params` and removed before `.mat` serialisation. All `fprintf` / `warning` calls in `init_coupled_params.m`, `init_coupled_state.m`, `build_coupled_lin_sys.m`, `compute_effective_body_brownian.m`, `adaptive_basis_check.m`, and `validate_simulation_params.m` now write to `simulation.log` via `logf_h` (with `warning` retained for MATLAB's alert system). Animation progress is logged to file every 10%. Default parameters unified: `E` → `3.279e8` Pa, `N_body` → 300 (including `config_example.json`); both JSON templates now include `compute_flowfield: true`. **Fully backward-compatible** — all logger handle parameters are optional with safe fallbacks; no physics/ODE/linear-system/Brownian path modified.

**v7.0.0:** **Free-torsion root coupling + constant motor torque model** — Per-pilus root boundary condition now supports `'free_torsion'` alongside the default `'clamped'`. In free-torsion mode, the root-segment twist angle ψ₁ is free to evolve under elastic + motor moments (θ₁ and φ₁ remain kinematically constrained). Motor torque model (`motor_torque.type ∈ {'off', 'constant'}`) applies a torsional moment at the root segment with optional stall-speed limiting (`max_omega_Hz`). **As of v7.3.0+,** the motor enters the pilus root **d₃ moment-balance row** only; explicit **body-side** reaction terms in the **system** angular-momentum rows `R_global(4:6)` were removed after audit (those rows already represent the full hydrodynamic moment on the coupled body+pilus assembly — see **v7.3.0** below). Releases v7.0–v7.2 experimented with `−τ·J_euler(3,:)'` and `−τ·d₃_lab` on `R_global(4:6)`; neither is retained in production. Eleven core functions modified at introduction: `build_coupled_lin_sys.m` (conditional root BC), `compute_pili_elastic_moments.m` (preserve κ₃(1) for free_torsion), `coupled_active_moment.m` + `ode_coupled_rhs.m` + `enforce_clamped_geometry.m` + `init_coupled_params.m` + `validate_simulation_params.m` + configuration I/O and runners. New `config_motor_flagella_example.json` template. **Clamped mode fully backward-compatible** — default behaviour unchanged.

**v7.3.0 (2026-05-04):** **Second physics audit — system vs. internal motor torque** — Removes any motor **reaction** contribution from `R_global(4:6)`. Those rows are the **full** body+pilus system angular-momentum balance (`B_global` rows 4–6 sum `r×f` with quadrature weights over every hydrodynamic point). The motor coupling between body and pilus is **internal**: action/reaction torques cancel in the system sum, so adding `−τ·d₃` or `−τ·J_euler(3,:)'` over-constrained the body's rotation. Production code applies **only** `+τ_eff` to the pilus segment-1 d₃ moment row; the body's angular response follows from the coupled Stokes solve and root kinematics (θ₁, φ₁ locked to body Ω). Also ships: README §11 κ correction; `pili_building` site-map comments; `custom_pili_building` `DeflectionAxisSwap` warning; `compute_plot_flowfield` shape checks; `validate_simulation_params` stall-curve warning for motors near stall; `coupled_background_flow` hot-path vectorisation; README §A.3 `ℓ_k` derivation for torsional resistance. **`free_torsion + motor`** trajectories differ from v7.0–v7.2; **`clamped`** and **`free_torsion` without motor** remain bit-exact vs. prior baselines where applicable.

**v7.3.1 (2026-05-04):** **Documentation & project website sync** — Version banners aligned across this guide, `README.md`, `CHANGELOG.md`, `requirements.txt`, `run_simulation.m`, `main_coupled_simulation.m`, and `docs/index.html`. §7 motor narrative and the project website are updated to match the v7.3.0 system-angular-momentum correction. No code or numerical changes relative to v7.3.0.

**v7.2.0 (2026-05-04):** **Audit-driven fixes for the v7.0 free-torsion + motor feature, plus a parametric JSON `pili[]` parsing fix.** (i) ~~Reaction torque on the body~~ — **Superseded by v7.3.0:** v7.2 replaced `−τ_eff · J_euler(3,:)'` with `−τ_eff · d₃_lab(seg 1)` on `R_global(4:6)`; v7.3.0 removes **both** patterns (see v7.3.0 note). (ii) Motor speed-limiting `1 − |ω|/ω_max` now uses a **`persistent DPSI1_PREV`** cache in `ode_coupled_rhs.m` (MATLAB structs are pass-by-value, so `params.dpsi1_prev` never persisted across ODE calls). Gated by `params.run_token` from `init_coupled_params.m`; `run_simulation.m` issues `clear ode_coupled_rhs` at integration start. (iii) `enforce_clamped_geometry` optional `force_full_projection`; `init_coupled_state` calls it with `true` so free_torsion ψ₁ matches the body frame at t = 0 when `alpha_0` or `initial_orientation` is non-trivial. (iv) `load_sim_config` parametric `pili[]` routed through `normalise_json_struct_array.m`. Demo timing in `config_motor_flagella_example.json` retuned. **clamped** path bit-exact; historical **free_torsion+motor** results from v7.0–v7.2 should be regenerated **again** after v7.3.0.

**v7.1.0 (2026-05-04):** **Custom-geometry hardening pass (I/O, diagnostics, pili kinematics, JSON configurability, flow-field renderer).** New helpers `parse_initial_orientation.m` (quaternion / axis_angle / euler_zyx → scalar-first unit quaternion; custom-mode only, overrides `alpha_0`) and `normalise_json_struct_array.m` (reconciles `jsondecode` cell-vs-struct dichotomy). New pre-shipped sample geometry under `geometry/` (`example_body_sphere.csv`, `example_body_capsule.csv`, `example_pilus_linear.csv`, `example_pilus_helix.csv`) plus `generate_example_geometry.m` so the inline + JSON custom examples run out of the box. **New strict requirement**: user-supplied pilus centerlines must be pre-discretised with equal-length segments (relative spread ≤ 1e-3), otherwise rejected at validation. `read_geometry_csv_mat.m` rewritten with `readmatrix` + `detectImportOptions` and pre-allocated fallback. `custom_body_building.m` adds `uniquetol` deduplication, a 15-point logarithmic alpha-shape scan with non-manifold rejection, an Euler-characteristic warning for non-closed surfaces, and `L_0_suggested` defaulting to the equivalent-sphere radius (`L_0_suggest_mode = 'volume' | 'extent'`). `custom_pili_building.m` rewritten to deduplicate, enforce segment uniformity, apply `deflection_z`/`deflection_x` in a Gram-Schmidt **local basis** (`e3` = aligned tangent → `deflection_z` is a twist about the tangent; `e1` = bending-plane normal → `deflection_x` is a tilt), store the **actual** post-transform first-segment tangent in `root_dir`, log the alignment angle, and warn if it exceeds 45°. **Naming caveat**: the parametric pili builder binds `deflection_x` to a twist and `deflection_z` to a tilt (axis roles swapped relative to custom); porting a non-zero deflection between modes requires swapping `deflection_x` ↔ `deflection_z`. Per-pilus extended physics (`intrinsic_curvature`, `active_moment`, `root_coupling`, `motor_torque`) now flow identically through the custom JSON/inline path via a shared `validate_pili_ext` helper and `save_sim_config.m` round-tripping. `compute_plot_flowfield.m` dispatches on `params.body_type` through a new `build_body_visualization` helper: custom bodies use an `alphaShape` `inShape` interior mask (bounding-sphere fallback), are rendered with `trisurf(params.body_faces, …)`, and have their 2D silhouettes extracted by a generic SDF-level contour — `major_axis` / `minor_axis` references are gone from the custom path. `main_coupled_simulation.m` inline custom block now explicitly promotes `cfg.custom_geometry.initial_orientation` to the top-level `cfg.initial_orientation`, matching the JSON path. **Compatibility**: parametric configurations are bit-exact; v6.x–v7.0.x custom configurations with non-uniform pili segments or non-zero custom-mode `deflection_*` need to be re-discretised or re-tuned per the above rules.

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

## 7. New in v7.0: Free-Torsion Root Coupling & Motor Torque

Per-pilus root boundary condition and active motor torque are configured in `main_coupled_simulation.m` (inline) or in JSON config. **`config_motor_flagella_example.json`** is a dedicated template for motor-driven flagella.

### Root coupling

```matlab
pili_config(1).root_coupling = 'clamped';         % Default: all 3 Euler angles locked
% pili_config(2).root_coupling = 'free_torsion';   % θ₁, φ₁ locked; ψ₁ free (motor-driven flagellum)
```

In JSON, under each `pili[]` entry: `"root_coupling": "clamped"` (default) or `"root_coupling": "free_torsion"`.

### Motor torque (only meaningful for `free_torsion`)

```matlab
pili_config(2).motor_torque = struct('type', 'off');           % Default: no motor
% pili_config(2).motor_torque = struct('type', 'constant', ...
%     'torque_pN_nm', 2000, 'max_omega_Hz', 300);             % 2000 pN·nm, stall at 300 Hz
```

In JSON: `"motor_torque": { "type": "off" }` or `"motor_torque": { "type": "constant", "torque_pN_nm": 2000, "max_omega_Hz": 300 }`.

Physical units are automatically converted to dimensionless at initialisation. For `'constant'` motor torque, the effective torque is linearly reduced as `|dψ₁/dt|` approaches `max_omega_Hz` (speed-limiting), using the previous-RHS-call `dψ₁/dt` cached in a **persistent** variable inside `ode_coupled_rhs.m` (reset between independent simulations via `params.run_token`). **There is no separate reaction torque added to the system angular-momentum rows `R_global(4:6)`** — those rows enforce the hydrodynamic moment balance for the **entire** coupled body+pilus discretisation; the motor is internal to that assembly, and the body's angular response follows implicitly from the coupled Stokes solve plus kinematic locking of θ₁ and φ₁ to the body frame (v7.3.0+; see *Notes on physics*).

### Interaction with `active_moment` and `intrinsic_curvature` (v7.2.0+)

When `root_coupling = 'free_torsion'` is combined with an `active_moment` whose `d3_amp` is non-zero, the d₃ component of the integrated active moment at segment 1 is **preserved** (it would otherwise be discarded by the clamped BC).  In that case the ψ₁ evolution at the root segment is driven by the **superposition** of three sources:

1. **Elastic torsion**: `κ₃(1) − κ₃₀(1)` (depends on `intrinsic_curvature`).
2. **Active d₃ moment**: the integrated travelling-wave / custom contribution along d₃, evaluated at segment 1.
3. **Motor torque**: `τ_eff` from `motor_torque.type = 'constant'`, with optional speed-limiting.

This is the intended design — it lets a single pilus model both an internal travelling-wave drive and a bulk motor — but users who only want a "pure motor" should set `active_moment.d3_amp = 0` (or leave `active_moment.type = 'off'`).

### Notes on physics (v7.3.0)

- **System angular-momentum rows vs. internal motor torque:** `R_global(4:6)` is paired with `B_global(4:6,:)`, which sums `r×f` (with quadrature weights) over **every** body and pilus hydrodynamic point — the **total** hydrodynamic moment on the coupled system. Motor torque exchanged between body and a pilus is **internal** to that system: it must **not** be added again as an explicit source on those rows. Releases **v7.0.x** used `−τ_eff · J_euler(3,:)'` and **v7.2.x** used `−τ_eff · d₃_lab(seg 1)` on `R_global(4:6)`; **v7.3.0 removes both** — only `+τ_eff` on the pilus segment-1 d₃ moment-balance row remains in production code.
- **What v7.2 still fixed (and remains in place):** persistent `DPSI1_PREV` speed-limiting cache, full ψ₁ projection at t = 0 for free_torsion when the body orientation is non-trivial, and parametric JSON `pili[]` normalisation.
- For a non-trivial initial body orientation (`alpha_0 ≠ 0` or a custom `initial_orientation`), `init_coupled_state.m` calls `enforce_clamped_geometry(..., true)` so that ψ₁ is projected to a self-consistent value at t = 0 even for free_torsion pili. ψ₁ then evolves freely from that starting point during the run.

---

## 8. Version and changes

- Current version is printed at the start of `simulation.log` and listed at the top of `README.md` / this file.
- See **`CHANGELOG.md`** for release notes.

---

## 9. License

MIT License — see `README.md` (license section).
