# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**Author:** Lucien (lucien-6@qq.com)

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

- **`compute_body_brownian_params.m`**: Deleted. Old body-only RSM friction matrix implementation fully superseded by `compute_effective_body_brownian.m`. Degradation test in `test_effective_brownian.m` updated to remove this dependency.

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

- **All `batch_sweep_*.m` drivers** (excluding `batch_run_from_json.m`): `case_name` tokens joined with **`-`**; **`enable_brown = false`**; **`T_obs = 200`**; **`generate_animation = true`**, **`animation_skip = 5`**; default pili Young’s modulus **`E = 3.279e8` Pa**.
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
