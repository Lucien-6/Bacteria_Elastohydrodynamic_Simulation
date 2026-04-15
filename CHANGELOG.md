# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**Author:** Lucien (lucien-6@qq.com)

---

## [5.4.4] — 2026-04-15

### Changed

- **`main_coupled_simulation.m`**: Inline example case updated to **`OnePili-Active-HighFre-Test`** with a single pilus at site 1 and a higher-frequency `traveling_wave` active moment, giving a simpler active-motion demonstration template.

### Fixed

- **`run_simulation.m` — MP4 animation export (`bacteria_3d_motion.mp4`)**: Hardened Windows HiDPI / display-scaling handling. After the figure is resized and probed, the encoded frame size is now clamped to a codec-safe, even-valued MPEG-4 box (`<= 3840 x 2160`). If `getframe` still returns oversized physical frames, `frame.cdata` is resized before `writeVideo`, preventing failures where logical figure shrinking alone did not reduce the captured physical pixel dimensions enough.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**: Version **5.4.4**, date **2026-04-15**; animation export notes updated to match the current HiDPI-safe implementation.

---

## [5.4.3] — 2026-04-15

### Changed

- **Header comment cleanup (all 47 `.m` files)**: Removed all `%  Updated:` history lines from file headers. Header comments now end with `Author:` / `Created:` only, keeping the source tree concise and relying on `CHANGELOG.md` for change history.

### Fixed

- **`build_coupled_lin_sys.m`**: Syntax comment now lists both 2-output and 3-output calling forms (`[lin_sys, geom]` and `[lin_sys, geom, hydro]`), matching the actual function signature.
- **`brown_rng_seed_from_clock.m`**: Added standard `Description:` / `Syntax:` / `Inputs:` / `Outputs:` header sections (previously had only a brief one-liner).
- **`post_processing/postprocess_cases.m`**: Standardised header comment to project style (`%POSTPROCESS_CASES`, two-space indent, `Syntax:` / `Inputs:` / `Outputs:` sections, `Author:` / `Created:`; was using inconsistent `% POSTPROCESS_CASES` + `Date:` format).

---

## [5.4.2] — 2026-04-15

### Added

- **`traveling_wave` active moment**: optional **`phase_rad`** [rad] (default 0) — constant phase added to \(\Phi = k_w s - \omega t\) in `coupled_active_moment.m`; set in `init_coupled_params` per pilus. JSON / `validate_simulation_params` accept a finite scalar.

---

## [5.4.1] — 2026-04-15

### Removed

- **`breaststroke` active moment mode** from `coupled_active_moment.m`, `init_coupled_params.m`, and `validate_simulation_params.m`. Per-pili active moments remain: `off`, `traveling_wave`, and `custom`. Second-pilus example in `config_example.json` now uses `traveling_wave`.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`docs/index.html`**, **`run_simulation.m`**: Dropped references to `breaststroke`.

---

## [5.4.0] — 2026-04-14

### Added

- **`coupled_active_moment.m` — `breaststroke` mode** (NEW): Two-harmonic asymmetric beat pattern modelling Chlamydomonas-like breaststroke flagellar motion. Moment density: $m(s,t) = A \cdot E(s) \cdot [\cos\Phi + C_2\cos(2\Phi+\varphi_2) + C_0]$ where $\Phi = k\,s - \omega\,t$ and $E(s) = (1-\beta)+\beta\,s$ is a linear spatial envelope. The linear envelope admits **exact antiderivatives** (same precision as `traveling_wave`). Three mechanisms absent from the single-harmonic `traveling_wave`:
  - **Temporal asymmetry** (`C2_ratio`, `C2_phase`): second harmonic creates distinct power/recovery stroke durations.
  - **Spatial envelope** (`envelope_slope`): amplitude ramps from base to tip.
  - **Mean curvature offset** (`mean_curv`): beat oscillates around a curved rest shape.
  Setting `envelope_slope = 0`, `C2_ratio = 0`, `mean_curv = 0` reduces exactly to `traveling_wave`.
- **`init_coupled_params.m`**: `case 'breaststroke'` — same `amplitude_dimless` / `omega_dimless` conversion as `traveling_wave` plus defaults for `envelope_slope` (0.5), `C2_ratio` (0.25), `C2_phase` (0), `mean_curv` (0).
- **`validate_simulation_params.m`**: Validation rules for breaststroke — `envelope_slope` in [0, 1], `C2_ratio >= 0`; warnings for `C2_ratio > 1`, `|mean_curv| > 10`.
- **`main_coupled_simulation.m`**: Commented inline breaststroke configuration example.
- **`config_example.json`**: Second pili entry with full breaststroke JSON template.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**: Version **5.4.0**, date **2026-04-14**.

---

## [5.3.2] — 2026-04-14

### Added

- **`compute_node_directors.m`** (NEW): Standalone utility that interpolates per-segment director frames to per-node directors. Internal nodes use neighbour averaging; boundary nodes use linear extrapolation with `min(2, N_k)` safety for single-segment pili. All averaged frames are projected onto the nearest proper rotation via **SVD polar decomposition** with an explicit `det < 0` correction to prevent reflections. Shared by `ode_coupled_rhs.m` (active moment evaluation) and `build_coupled_lin_sys.m` (geometry reconstruction).

### Fixed

- **`build_coupled_lin_sys.m` — node director interpolation**: The inline implementation lacked the `det(R) < 0` SVD fix, so polar decomposition could silently return an improper rotation (reflection) instead of the nearest proper rotation. Also had a potential index-out-of-bounds for single-segment pili (`N_k = 1`) in the boundary extrapolation. Both issues fixed by replacing inline code with a call to the new shared `compute_node_directors.m`.

### Changed

- **`ode_coupled_rhs.m`**: (1) Elastic-moment assignment to `R_global` and active-moment subtraction **vectorised** (replaced per-element `for` loops with `reshape` + array indexing). (2) Local `compute_node_directors` helper removed; calls shared `compute_node_directors.m` instead — eliminates code duplication with `build_coupled_lin_sys.m`.
- **`validate_simulation_params.m`**: Added warnings for `traveling_wave` active moment parameters: negative `wavenumber` (wave propagates tip → root) and non-zero `d3_amp` (drives torsion rather than bending, unusual for standard pili).

### Documentation

- **`coupled_active_moment.m`**: Clarified `d3_amp` physical meaning — non-zero value drives twist (torsional component), not bending; typically zero for standard flagellar/pili beating.
- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**: Version **5.3.2**, date **2026-04-14**.

---

## [5.3.1] — 2026-04-14

### Fixed

- **`run_simulation.m` — MP4 animation (`bacteria_3d_motion.mp4`)**: Two failure modes fixed:
  1. **Large Brownian trajectories** push computed figure width beyond Windows MPEG-4 / H.264 encoder limits. Logical figure size is now **proportionally scaled** to fit within a **3840 x 2160** codec-safe box.
  2. **Hi-DPI / display scaling** (e.g. 125%): `getframe` returns **physical pixels** (e.g. 2250 x 1366) larger than the logical `Position` (1800 x 1200). The codec limit was checked against logical pixels only, which was ineffective. Fix: after creating the figure, a **probe `getframe`** measures the actual physical size; if it exceeds codec limits, the figure `Position` is **shrunk** and re-probed. `VideoWriter` is opened **after** the physical frame size is known. In the render loop, `getframe(fig)` is called **without a rect** (avoids logical/physical coordinate mismatch). A fallback **`normalize_anim_frame_cdata`** (`imresize` or `interp2`) catches any residual per-frame size drift.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**: Version **5.3.1**, date **2026-04-14**.

---

## [5.3.0] — 2026-04-13

### Added

- **`batch sweep/batch_sweep_N_body_eps_RSM.m`**: Nested sweep of **`N_body`** and **`eps_body = k · h_mean`** (RSM regularisation) on a sphere with no pili; compares steady **`V_sink`** to analytic Stokes settling speed; writes `batch_sweep_N_body_eps_RSM_grid.csv`, `batch_sweep_N_body_eps_RSM_summary.csv`, and `.mat` results.

### Changed

- **Default body RSM regularisation** (`init_coupled_params.m`): **`eps_body`** default is now **`h_mean / 2.6667`** (was **`h_mean / 3`**), i.e. \(\varepsilon/\Delta s_{\mathrm{avg}} \approx 0.375\), from the above sweep vs Stokes reference.
- **`build_global_resistance.m`**: Header comment aligned with the new default ratio.

### Documentation

- **`README.md`** (§8.5, JSON/`eps_body` notes, feature history), **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**: Version **5.3.0**, date **2026-04-13**.

---

## [5.2.1] — 2026-04-10

### Changed

- **Default `N_body`**: project-wide default for body surface Stokeslet points is now **200** (was **250**) — `load_sim_config.m`, `init_coupled_params.m`, `main_coupled_simulation.m`, `config_example.json`, `batch sweep/*.m`, `docs/index.html`, `README.md`.
- **Default `tlim`**: ODE segment timeout limit default is now **360** s (was **120** s) in the same locations.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**: Version **5.2.1** (log banner).

---

## [5.2.0] — 2026-04-10

### Added

- **Brownian RNG seed policy** (`run_simulation.m`, `brown_rng_seed_from_clock.m`): When `enable_brown` is true, MATLAB's global generator is initialised with **`rng(seed, 'twister')`**. Default seed is **local wall-clock time** encoded as integer **`MMddHHmmss`** (via `brown_rng_seed_from_clock`). Optional **`cfg.brown_rng_seed`** (positive scalar) or JSON **`brownian.rng_seed`** forces a user-defined reproducible seed.
- **Seed persistence**: The realised seed is written to **`output/<case>/config.json`** (`brownian.rng_seed`), **`brownian_rng_seed.txt`**, **`simulation.log`** (`[PARAMS] Brownian RNG`), **`results.brown_rng_seed`**, and **`params.brown_rng_seed`** in `simulation_results.mat`.
- **`load_sim_config.m` / `save_sim_config.m`**: Parse and emit optional **`brownian.rng_seed`** ↔ flat **`cfg.brown_rng_seed`**.

### Changed

- **`batch sweep/batch_repeat_brownian.m`**: Brownian seeds are applied only inside **`run_simulation`** (not via a separate `rng()` call). **`brown_rng_seed_mode`**: `'default'` omits **`cfg.brown_rng_seed`** so each repeat uses the same **MMddHHmmss** time default as a single run (`brown_rng_seed_from_clock`); `'fixed'` sets **`cfg.brown_rng_seed = rng_seed_base + k`** for reproducible batches.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**: Version **5.2.0**, date **2026-04-10**; Brownian seed and output files documented.

---

## [5.1.1] — 2026-04-08

### Changed

- **Default discretisation**: project-wide defaults for **`N_body`** (body surface Stokeslet points) and **`seg_per_um`** (RSS segments per micron along each pilus) are now **250** and **10** respectively (`load_sim_config.m`, `init_coupled_params.m`, `main_coupled_simulation.m`, JSON templates, `batch sweep/*.m`, `docs/index.html`, `README.md`).
- **`run_simulation.m`**: After `exportgraphics` saves `summary_plot.png`, the summary figure is **`close`d** so batch runs do not accumulate visible figure windows.
- **`compute_plot_flowfield.m`**: After `exportgraphics` saves `flowfield_panels.png`, the flow-field figure is **`close`d** for the same reason.
- **`run_simulation.m`**: Log banner reports **Version 5.1.1**.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**: Version **5.1.1**, date **2026-04-08**; §17 notes batch-friendly figure cleanup.
- **JSON / config docs**: **`config_example.json`** and **`configs/OnePili-P1-L3.json`** updated to include **`background_flow`** (v5.0) and per-pilus **`intrinsic_curvature`** / **`active_moment`** defaults; **`README.md`** §2.2 documents nested JSON keys vs flat `cfg`; **`USER_GUIDE.md`** §6 adds JSON fragments; **`docs/index.html`** quick-start JSON aligned with `save_sim_config`, **`traveling_wave`** field names corrected to **`d1_amp`** / **`d2_amp`** / **`d3_amp`**.

---

## [5.1.0] — 2026-04-07

### Added

- **`post_processing/postprocess_cases.m`**: Batch post-processing of **multiple** simulation result folders under a chosen parent directory (non-recursive: only first-level subfolders containing `timeseries.csv`). Reads each case's `timeseries.csv` and `config.json`, computes summary metrics (sinking velocity statistics, maximum horizontal displacement, 3-D trajectory tortuosity, translational diffusion coefficient via MSD fit `6D·τ + V²·τ²`), and exports a **Summary** sheet to Excel (`.xlsx`). On Windows, optional COM formatting (bold headers, borders, row shading) and optional user prompt before opening the file. English UI dialogs; progress via `waitbar` and console output.

### Changed

- **`run_simulation.m`**: Log banner reports **Version 5.1.0**.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**: Version **5.1.0**, date **2026-04-07**; post-processing workflow documented.

---

## [5.0.0] — 2026-04-03

### Added

- **`coupled_background_flow.m`** (NEW): Evaluates a steady (time-invariant) background flow field at all hydrodynamic evaluation points (body surface + pili nodes). Supports three modes: `'off'` (default, zero flow), `'uniform'` (spatially constant velocity), and `'shear_z'` (piecewise-linear shear). Output vector matches the DOF layout of `build_coupled_lin_sys` (translation block + twist block). All values pre-converted to dimensionless during initialisation.
- **`coupled_active_moment.m`** (NEW): Computes per-pilus time-varying active internal moments at the current dimensionless time. Supports `'off'` (default), `'traveling_wave'` (sinusoidal traveling wave with configurable amplitude, frequency, wavenumber, and per-component weights), and `'custom'` (user-provided function handle). Algorithm follows `3D-Flexible-Filaments/integrated_internal_moments.m` (antiderivative → segment integrals → reverse cumsum → director-basis projection).
- **Intrinsic curvature configuration**: Per-pilus static intrinsic curvature (`params.kb0{k}`) can now be set to `'frozen'` (default, keep initial geometry), `'zero'`, `'constant'` (user-specified `[k1, k2, k3]`), or `'helix'` (helical curvature from `kappa_d2`/`kappa_d3`). Configured via `pili_config(k).intrinsic_curvature` in the inline config or per-pilus JSON fields.
- **Active moment configuration**: Per-pilus time-dependent active moments configured via `pili_config(k).active_moment` (struct with `.type`). Pre-converts physical parameters (amplitude [N·m], frequency [Hz]) to dimensionless at initialisation.
- **Background flow configuration**: Global steady background flow configured via `cfg.background_flow` (struct with `.type`, `.U_max`, `.direction`, `.shear_axis`, `.shear_height_m`). Pre-converts physical velocity to dimensionless at initialisation.

### Changed

- **`build_coupled_lin_sys.m`**: Added optional 3rd output `hydro` struct (contains `B_global`, `M_weighted`, `A_twist_total`, `M_rotlet`, `n_trans`) for background flow computation. Backward-compatible: existing 2-output calls are unaffected (`nargout >= 3` guard).
- **`ode_coupled_rhs.m`**: Integrated active moments (Step 8b: subtracted from `R_global` per-pilus before clamped BC) and background flow (Step 8c: `B*Ainv*u_bg` subtracted from `R_global`). Added local helpers `apply_Ainv_u` and `compute_node_directors`.
- **`init_coupled_params.m`**: Added three Name-Value pairs (`'intrinsic_curvature'`, `'active_moment'`, `'background_flow'`). Intrinsic curvature overrides `params.kb0{k}` per mode; active moment configs pre-convert physical units to dimensionless; background flow pre-computes `U_max_dimless` and `shear_height_dimless`. Summary printout extended with new feature status per pilus.
- **`parse_pili_config.m`**: Added optional 2nd output `pili_ext` (struct array with `.intrinsic_curvature` and `.active_moment` per pilus). Backward-compatible: 1-output calls still work; `Pili_Matrix` [8×M] unchanged.
- **`load_sim_config.m`**: Parses `background_flow` JSON section and collects `pili_ext` from `parse_pili_config`.
- **`save_sim_config.m`**: Writes `background_flow` section and per-pilus `intrinsic_curvature`/`active_moment` fields to JSON. Function handles (`'custom'` mode) are stripped before serialisation.
- **`validate_simulation_params.m`**: Added validation rules for background flow type/parameters, per-pilus intrinsic curvature type/kappa, per-pilus active moment type/amplitude/frequency; warnings for background flow + Brownian combination and active moments without pili.
- **`main_coupled_simulation.m`**: Added inline config examples for background flow, per-pilus intrinsic curvature, and per-pilus active moment (all default to `'off'`/`'frozen'`). Uses `[cfg.Pili_Matrix, cfg.pili_ext] = parse_pili_config(pili_config)`.
- **`run_simulation.m`**: Extracts new feature fields from `cfg`, passes `'intrinsic_curvature'`, `'active_moment'`, `'background_flow'` to `init_coupled_params`; adds new fields to `val_cfg` for validation. Log banner reports **Version 5.0.0**.

### Compatibility

- All three new features default to `'off'`/`'frozen'`: existing configurations and sweep scripts produce **identical results** without modification.
- Background flow is deterministic and does not affect the Brownian diffusion tensor (`compute_effective_body_brownian` probes `lin_sys` which excludes flow contributions).
- Active moments are applied before clamped BC rows are zeroed, avoiding any boundary condition conflict.
- `adaptive_basis_check.m` preserves `params.am_cfg{k}` and `params.bg_flow` across basis rotations (static config, not state-dependent).

---

## [4.1.0] — 2026-04-03

### Added

- **`batch sweep/batch_repeat_brownian.m`**: Batch independent repeats of the **same** Brownian-on configuration. User sets repeat count and `case_name_base`; each run writes to `output/<case_name_base>/001/`, `002`, … (zero-padded). Requires `enable_brown = true`; forces `preview_initial = false` and `auto_confirm = true`; sets `rng(rng_seed_base + k)` before each repeat for reproducible independent noise streams. Optional JSON config via `config_json`, else inline `cfg_base` in the script. Writes `batch_repeat_summary_*.txt` under `output/<case_name_base>/`.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**: Version **4.1.0**, date **2026-04-03**; file tree and feature bullets updated for the new batch script.
- **`docs/index.html`**: Version badge/footer aligned with release **4.1.0**.

### Changed

- **`run_simulation.m`**: Log banner reports **Version 4.1.0** (aligned with `README.md`).

---

## [4.0.0] — 2026-03-31

### Added

- **`build_coupled_lin_sys.m`** (NEW): Extracted Steps 1-7 from `ode_coupled_rhs.m` into a reusable subroutine that assembles the coupled body-pili linear system with clamped BC. Shared by the ODE RHS and the Brownian mobility computation.
- **`compute_effective_body_brownian.m`** (NEW): Computes the full-system effective 6×6 body mobility by probing the coupled linear system, dimensionalises it, and returns the Cholesky factor for Euler-Maruyama Brownian dynamics. Implicitly captures pili hydrodynamic drag via Schur complement.
- **`fixman_midpoint`** configuration option: Optional Fixman midpoint correction for configuration-dependent diffusion drift (2× mobility evaluation per Brownian step).

### Changed

- **`ode_coupled_rhs.m`**: Refactored to call `build_coupled_lin_sys` for Steps 1-7, retaining only elastic moment computation and state-rate packing.
- **`run_simulation.m`**: Brownian kick block replaced: uses per-step dynamically computed `L_eff` (6×6 Cholesky) from full-system effective friction; generates 6D correlated noise (translation-rotation coupling included); supports optional Fixman midpoint correction.
- **`init_coupled_params.m`**: Removed one-time body-only `compute_body_brownian_params` call at initialisation; Brownian diffusion is now computed dynamically per step; added `fixman_midpoint` parameter.
- **`load_sim_config.m`**: Parses new `fixman_midpoint` field from the `brownian` JSON block.
- **`config_example.json`**: Added `fixman_midpoint: false` to the `brownian` block.
- **`validate_simulation_params.m`**: Added validation for `fixman_midpoint` parameter.
- **`Brownian_Implementation_Plan.md`**: Added upgrade notice documenting the transition from body-only to full-system effective friction.
- **`README.md`**: Updated section 16, file listings, architecture diagram, and feature descriptions.

### Removed

- **`compute_body_brownian_params.m`**: Deleted. Old body-only RSM friction matrix implementation fully superseded by `compute_effective_body_brownian.m`.
- **`test_effective_brownian.m`**: Deleted. Verification tests passed (7/7); script no longer needed.

### Fixed

- **`main_coupled_simulation.m`**: Added missing `cfg.fixman_midpoint` to inline configuration block (runtime crash when using inline config path).
- **`configs/OnePili-P1-L3.json`**: Added missing `fixman_midpoint` field to the `brownian` block.

---

## [3.3.5] — 2026-03-30

### Changed

- **`batch sweep/batch_sweep_major_axis.m`**: semi-major axis sweep grid **`0.4:0.2:1.6` µm** (7 cases); **`minor_axis`** fixed at **0.4 µm** (first case is a sphere, `major_axis == minor_axis`).

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**: version **3.3.5**, date **2026-03-30**; feature table row 56.
- **`run_simulation.m`**: log banner **Version 3.3.5**.

---

## [3.3.4] — 2026-03-30

### Removed

- **`batch sweep/batch_sweep_pili_site.m`**: removed; attachment-site sweeps use **`batch_sweep_pili_length.m`** with user-defined **`site_list`**.

### Changed

- **`batch sweep/batch_sweep_major_axis.m`**: **body-only** major-axis sweep — **`Pili_Matrix` empty** (no pili).

### Documentation

- **`README.md`**: file counts (38 `.m` files; 6 sweep drivers + `batch_run_from_json`), feature table row 55; `batch_sweep_major_axis` described as no-pili body.

---

## [3.3.3] — 2026-03-30

### Changed

- **`batch sweep/batch_sweep_alpha.m`**: **Body-only** initial tilt sweep — **`Pili_Matrix` empty** (no pili); fixed ellipsoidal capsule geometry unchanged.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**: version **3.3.3**, date **2026-03-30**; feature table row 54; `batch_sweep_alpha` description clarified.
- **`run_simulation.m`**: log banner **Version 3.3.3**.

---

## [3.3.2] — 2026-03-30

### Added

- **`batch sweep/batch_sweep_N_body.m`**: sweep `N_body` for a **spherical** body (radius **0.4 µm**, `major_axis == minor_axis`) with **no pili** (`Pili_Matrix` empty).
- **`batch sweep/batch_sweep_seg_per_um.m`**: sweep `seg_per_um` (default grid 10:5:35).

### Removed

- **`batch sweep/batch_sweep_pili_count.m`**: removed (multi-pili / site control covered by `batch_sweep_pili_length.m` + JSON batching).

### Changed

- **All `batch_sweep_*.m` drivers** (excluding `batch_run_from_json.m`): `case_name` tokens joined with **`-`**; **`enable_brown = false`**; **`T_obs = 200`**; **`generate_animation = true`**, **`animation_skip = 5`**; default pili Young's modulus **`E = 3.279e8` Pa**.
- **`batch_sweep_pili_length.m`**: user-defined **`site_list`** for arbitrary per-pilus attachment sites (equal length per run).
- **`batch_sweep_pili_E.m`**: E sweep **`1e5`–`1e9` Pa**, tenfold per step (`10.^(5:9)`).
- **`validate_simulation_params.m`**: allow **`major_axis == minor_axis`** (sphere); error only when **`major_axis < minor_axis`**.
- **`README.md`**: batch sweep inventory (7 sweep scripts + JSON runner), file counts, `N_body` script description.

---

## [3.3.1] — 2026-03-30

### Added

- **`CHANGELOG.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**: release notes, concise user guide, and MATLAB-only environment note (no pip dependencies).

### Changed

- **`README.md`**: Version/date bump; §18 flow-field description; §19 file tree; new §25 companion docs; feature table rows 51–52.
- **Flow-field 2D quiver** (`compute_plot_flowfield.m`): Removed the white underlay (double-draw halo) on XZ and YZ panels; arrows are drawn once in dark gray for clearer interpretation on the jet colormap. Arrow length, line width, and head size are scaled **+50%** relative to the previous dark layer for improved visibility.
- **`run_simulation.m`**: Log banner reports **Version 3.3.1** (aligned with `README.md`).

---

## [3.3.0] — 2026-03-18

### Added

- Parameter sweep scripts consolidated under `batch sweep/` (six sweep drivers + shared patterns).
- `batch_run_from_json.m`: batch runner over a folder of JSON configs with summary log.

### Fixed

- Minor batch-script hygiene (`batch_run_from_json.m` dead code / redundant assignments).

---

## Earlier releases

Feature and fix history prior to 3.3.0 is summarized in `README.md` (sections *22. 验证总结* — bugfix and feature tables).
