# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**Author:** Lucien (lucien-6@qq.com)

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
