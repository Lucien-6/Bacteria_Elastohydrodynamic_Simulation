# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**Author:** Lucien (lucien-6@qq.com)

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
