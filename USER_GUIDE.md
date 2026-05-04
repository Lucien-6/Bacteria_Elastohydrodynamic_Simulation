# User Guide — Coupled Body-Pili Simulation

English | [中文](USER_GUIDE.zh-CN.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Version:** 7.3.1  
> **Updated:** 2026-05-04

This guide is for running simulations, choosing the right configuration path, finding outputs, and post-processing completed cases. For model equations and derivations, see [`docs/theory.md`](docs/theory.md). For release history, see [`CHANGELOG.md`](CHANGELOG.md).

## 1. Before You Run

### Requirements

- MATLAB R2020b or later is recommended.
- PNG and PDF exports use `exportgraphics`; optional interactive HTML exports require MATLAB R2026a or later.
- Capsule body surface generation requires Anton Semechko's **S2 Sampling Toolbox** from MATLAB Add-On Explorer / File Exchange.
- No Python, Conda, or pip environment is required. [`requirements.txt`](requirements.txt) documents this explicitly.

### Recommended Working Folder

Open MATLAB from the repository root or run:

```matlab
cd('path/to/Coupled_Body_Pili')
```

This keeps relative paths for JSON configs, geometry files, and outputs predictable.

## 2. Choose a Run Mode

| Goal | Recommended path |
| --- | --- |
| Quickly test or modify the default demo | Edit [`main_coupled_simulation.m`](main_coupled_simulation.m). |
| Reproduce or share a case | Use a JSON file copied from [`config_example.json`](config_example.json). |
| Import custom body and pili geometry | Use [`config_custom_example.json`](config_custom_example.json). |
| Test a free-torsion motor-driven flagellum | Use [`config_motor_flagella_example.json`](config_motor_flagella_example.json). |
| Run many cases without prompts | Call [`run_simulation.m`](run_simulation.m) with `cfg.auto_confirm = true`. |

## 3. Run with Inline Parameters

1. Open [`main_coupled_simulation.m`](main_coupled_simulation.m).
2. Edit the `cfg` block and `pili_config` entries.
3. Run:

```matlab
main_coupled_simulation
```

Use this mode when you are developing a new case interactively and want the configuration close to the entry script.

## 4. Run with JSON Configuration

1. Copy [`config_example.json`](config_example.json) to a local folder, for example `configs/my_case.json`.
2. Set `config_file` in [`main_coupled_simulation.m`](main_coupled_simulation.m):

```matlab
config_file = 'configs/my_case.json';
main_coupled_simulation
```

The repository ignores `configs/` by default, so local paths and exploratory parameters are not accidentally committed.

### Required JSON Sections

A complete JSON config should include these top-level sections:

- `case_name`
- `geometry_mode`
- `body`
- `pili`
- `material`
- `fluid`
- `simulation`
- `brownian`
- `output`
- `flowfield`
- `background_flow`

The authoritative schemas are the tracked example files and the loader implementation in [`load_sim_config.m`](load_sim_config.m).

## 5. Programmatic and Batch Runs

For scripts or sweeps, call the simulation engine directly:

```matlab
cfg = load_sim_config('config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

Batch sweep drivers live in [`batch sweep/`](batch%20sweep/). The folder includes parameter sweeps for body resolution, body regularization, pilus length, pilus stiffness, segment density, initial angle, and Brownian repeats.

## 6. Custom Geometry Mode

Use custom mode when both the body and pili are imported from point-cloud files.

1. Copy [`config_custom_example.json`](config_custom_example.json).
2. Set `"geometry_mode": "custom"`.
3. Provide body and pili files under `"custom_geometry"`.
4. Use the example files in [`geometry/`](geometry/) as format references.

### Custom Geometry Rules

- Custom body and custom pili must be used together.
- Body point clouds support `.csv`, `.txt`, and `.mat` files with three spatial coordinates.
- Pili centerlines must be pre-discretized into equal-length segments; the accepted relative segment-length spread is at most $10^{-3}$.
- Consecutive duplicate or near-duplicate points are removed before validation.
- For arbitrary body orientation, prefer `initial_orientation` over `alpha_0`.
- Custom-mode `deflection_z` is a twist about the aligned tangent, while custom-mode `deflection_x` is a tilt. In parametric mode these geometric roles are swapped.

## 7. Motor and Free-Torsion Runs

Use [`config_motor_flagella_example.json`](config_motor_flagella_example.json) as the starting template.

Key fields under each pilus entry:

- `root_coupling`: `"clamped"` or `"free_torsion"`.
- `motor_torque.type`: `"off"` or `"constant"`.
- `motor_torque.torque_pN_nm`: motor torque in pN nm.
- `motor_torque.max_omega_Hz`: optional stall-speed limit in Hz.

For `free_torsion`, the root twist angle $\psi_1$ evolves dynamically. The motor torque is applied to the pilus root $d_3$ moment-balance row only; no separate reaction torque is added to `R_global(4:6)` because those rows represent the angular-momentum balance of the full coupled body-pilus system.

## 8. Brownian Dynamics

Set `brownian.enable` in JSON or `cfg.enable_brown` in MATLAB.

Important options:

- `brownian.rng_seed`: fixed positive seed for reproducible runs.
- `brownian.dt_B_ratio`: Brownian time step relative to output step.
- `brownian.fixman_midpoint`: enables the Fixman midpoint correction for configuration-dependent mobility.

When Brownian dynamics are enabled, the realized seed is saved in `config.json`, `brownian_rng_seed.txt`, `simulation.log`, and `simulation_results.mat`.

## 9. Main Outputs

Outputs are written under `output/<case_name>/` or a timestamped run directory.

| File | Purpose |
| --- | --- |
| `simulation_results.mat` | Full `results` and `params` structures. |
| `timeseries.csv` | Scalar time series for external analysis. |
| `config.json` | Saved configuration snapshot for reproduction. |
| `summary_plot.png` / `.pdf` / optional `.html` | Main summary figure. |
| `flowfield_panels.png` / `.pdf` / optional `.html` | XZ and YZ flow-field panels. |
| `initial_pose.png` / `.pdf` / optional `.html` | Initial-pose preview when enabled. |
| `bacteria_3d_motion.mp4` | Optional 3D animation. |
| `simulation.log` | Complete run log. |
| `brownian_rng_seed.txt` | Brownian seed and reproduce hint when Brownian dynamics are enabled. |

## 10. Post-Processing Completed Runs

To aggregate multiple completed case folders:

1. In MATLAB, run `addpath('post_processing')`.
2. Run:

```matlab
postprocess_cases()
```

1. Select the parent folder containing completed case folders.
2. Select the case folders to include.
3. Choose the output `.xlsx` path.

The post-processor reads `timeseries.csv` and `config.json`, computes summary kinematic metrics, writes an Excel workbook, and exports a summary figure as PNG, PDF, and optional HTML.

## 11. Common Checks

| Symptom | Check |
| --- | --- |
| MATLAB cannot find a function | Confirm the repository root is on the MATLAB path. |
| Capsule generation fails | Confirm the S2 Sampling Toolbox is installed and visible on the MATLAB path. |
| JSON case runs without pili | Check that `pili` is a valid JSON array and that each pilus entry has the required fields. |
| Custom geometry is rejected | Check point-cloud orientation, duplicate points, and equal-length pilus segments. |
| Output folder already exists | Use `auto_confirm = true` for unattended rename behavior, or choose interactively. |
| Brownian repeat results differ | Confirm whether the seed mode is fixed or timestamp-based. |

## 12. Documentation Map

- [`README.md`](README.md): project entry page.
- [`README.zh-CN.md`](README.zh-CN.md): Chinese project entry page.
- [`USER_GUIDE.zh-CN.md`](USER_GUIDE.zh-CN.md): Chinese version of this guide.
- [`docs/theory.md`](docs/theory.md): theory and numerical model.
- [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md): Chinese theory document.
- [`CHANGELOG.md`](CHANGELOG.md): release notes and compatibility changes.
- [`docs/index.html`](docs/index.html): project website.

## 13. License

MIT License. See [`LICENSE`](LICENSE).
