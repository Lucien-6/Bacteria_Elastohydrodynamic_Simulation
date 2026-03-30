# User Guide — Coupled Body-Pili Simulation

> **Author:** Lucien (lucien-6@qq.com)  
> **Version:** 3.3.5  
> **Updated:** 2026-03-30  

This guide is a short entry point. **Full theory, parameters, and file map:** see `README.md`.

---

## 1. What you need

- **MATLAB** (R2020b or later recommended). No Python environment is required to run the simulator (`requirements.txt` documents this).

---

## 2. How to run

### Option A — Inline parameters (default)

1. Open `main_coupled_simulation.m`.
2. Edit the `cfg` block (case name, geometry, pili, time span, etc.).
3. Run the script.

### Option B — JSON configuration

1. Copy `config_example.json` or use a file under `configs/`.
2. In `main_coupled_simulation.m`, set `config_file` to your JSON path.
3. Run.

### Option C — Programmatic / batch

```matlab
cfg = load_sim_config('config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

Use `batch sweep/` scripts or `batch_run_from_json.m` for unattended sweeps.

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

---

## 4. Flow-field figure (`flowfield_panels.png`)

- **Layout:** 2×2 — 3D XZ / 3D YZ / 2D XZ / 2D YZ.
- **2D panels:** Velocity magnitude on a `jet` scale + **normalized** direction arrows (equal length per arrow; magnitude is read from the color field).
- **v3.3.1:** 2D arrows are single dark-gray vectors (no white halo). Arrows are drawn **50% larger** than the previous dark layer for readability.
- **v3.3.2:** Batch sweep scripts under `batch sweep/` were refreshed (new `N_body` / `seg_per_um` drivers, removed `batch_sweep_pili_count.m`); see `CHANGELOG.md`.
- **v3.3.3:** `batch_sweep_alpha.m` sweeps **`alpha_0` for a body-only capsule** (empty `Pili_Matrix`, no pili).
- **v3.3.4:** Removed `batch_sweep_pili_site.m`; use **`batch_sweep_pili_length.m`** with **`site_list`** for attachment-site studies.
- **v3.3.5:** `batch_sweep_major_axis.m` sweeps **semi-major axis** from **0.4 to 1.6 µm** (step **0.2 µm**), body-only (no pili).

Tunable via `cfg`: `flowfield_grid_spacing`, `flowfield_grid_range`, `flowfield_quiver_skip`, `flowfield_chunk_size`.

---

## 5. Version and changes

- Current version is printed at the start of `simulation.log` and listed at the top of `README.md` / this file.
- See **`CHANGELOG.md`** for release notes.

---

## 6. License

MIT License — see `README.md` (license section).
