# Bacteria Elastohydrodynamic Simulation

English | [中文](README.zh-CN.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Created:** 2026-02-06  
> **Updated:** 2026-11-04
> **Version:** 8.7.2
> **Runtime:** MATLAB-only. No Python or pip dependencies are required.

This project simulates the low-Reynolds-number settling dynamics of a rigid bacterial body coupled to one or more flexible pili or flagella. The body is represented by regularized Stokeslets on a surface point cloud, while each pilus is represented as an inextensible Kirchhoff elastic rod coupled to the surrounding Stokes flow through regularized Stokeslet segments.

## Key Capabilities

- Rigid capsule or imported custom body geometry.
- Parametric or imported custom pili centerlines.
- Coupled RSM/RSS hydrodynamics with non-local segment-integral twist mobility, reciprocal Schur-complement coupling, and physical standard-rotlet scaling, and a calibrated pilus hydrodynamic radius (`material.hydro_radius_factor`).
- Kirchhoff elastic rods with intrinsic curvature and active moment options.
- Optional body-only Brownian dynamics: the 6-DOF body covariance uses the full coupled body-pili effective mobility with dimensionless SPD repair and an isolated reproducible random stream; each interval's increment is realised as a constant thermal generalised force acting through the coupled ODE (implicit Ito drift, no root kink). Pili internal elastic modes receive no thermal noise of their own.
- Optional background flow: off, uniform, or finite-thickness shear.
- Free-torsion root coupling and constant motor torque for motor-driven flagella.
- Inline MATLAB configuration, JSON configuration, and programmatic batch execution.
- Publication-oriented exports: MAT, CSV, log, PNG, vector PDF, optional interactive HTML, and MP4 animation.
- Optional **parallel** batch drivers ([`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m), [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m)): `parfor` with automatic pool sizing (Parallel Computing Toolbox); serial fallback when the toolbox or pool is unavailable.
- **Linux HPC cluster execution** (new in 7.10.0): [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) / [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m) is a maintained fork hardened for headless, single-node-exclusive Linux clusters: platform-aware video export, Linux-native CPU/memory probing, node-local `parpool` job storage, heavy-output sampling to bound per-worker memory, and live sub-case progress reporting. The original Windows/interactive scripts are unaffected.

## Requirements

- MATLAB R2020b or later is recommended.
- PNG and PDF figure export use `exportgraphics` and work from MATLAB R2020a onward.
- Optional interactive HTML figure export requires MATLAB R2026a or later.
- Capsule body-surface generation requires the MATLAB Add-On / File Exchange **S2 Sampling Toolbox** by Anton Semechko. Install it so `ParticleSampleSphere`, `SpiralSampleSphere`, `IcosahedronMesh`, `QuadCubeMesh`, and `SubdivideSphericalMesh` are on the MATLAB path.
- Python / Conda / pip are not required. See [`requirements.txt`](requirements.txt).
- **Parallel Computing Toolbox** is optional: only needed for the default parallel branches of [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m) and [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m) (`use_parallel = true`).

## Quick Start

Before running a new case, confirm that the project root is the current MATLAB folder or is on the MATLAB path. The main script adds the project root defensively, but starting from the repository root keeps relative JSON paths and local output folders predictable.

### Option A: Inline Parameters

1. Open [`main_coupled_simulation.m`](main_coupled_simulation.m).
2. Edit the `cfg` block, including geometry, pili, material, simulation time, Brownian options, and output options.
3. Run:

```matlab
main_coupled_simulation
```

### Option B: JSON Configuration

1. Copy [`configs/config_example.json`](configs/config_example.json) to a local configuration path, for example `configs/my_case.json`.
2. Set `config_file` in [`main_coupled_simulation.m`](main_coupled_simulation.m):

```matlab
config_file = 'configs/my_case.json';
main_coupled_simulation
```

Use [`configs/config_custom_example.json`](configs/config_custom_example.json) for imported body and pili geometry, and [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json) for free-torsion motor examples.  A realistic single-flagellum *E. coli* swimmer (2.5 µm × 0.88 µm body, left-handed helical filament with pitch 2.22 µm, motor torque 1300 pN·nm) is provided in [`configs/config_ecoli_single_flagellum.json`](configs/config_ecoli_single_flagellum.json).

For machine-specific cases, keep copied JSON files under a local `configs/` path (for example `configs/my_case.json`). The repository tracks example JSON files under `configs/`; only `output/` is listed in `.gitignore`.

### Option C: Programmatic Use

```matlab
cfg = load_sim_config('configs/config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

Batch sweep **driver scripts** live in [`batch_sweep/`](batch_sweep/). Parameter-sweep CSV/MAT products are written under `batch sweep/` (with a space) when those drivers run. To batch-run every `*.json` in a folder (optional parallel, duplicate `case_name` skip, real-time `batch_summary_*.log`), use [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m). Multi-case post-processing is available in [`post_processing/postprocess_cases.m`](post_processing/postprocess_cases.m): 3-D translational MSD fit (\(\mathrm{MSD}=6D\tau+V^2\tau^2\) on \(0<\tau\le T/2\)) and XZ-plane rotational diffusion from MSAD (\(\mathrm{MSAD}=2D_\mathrm{rot}\tau\)), exported to Excel plus a PNG/PDF/(optional HTML) summary figure. To re-render only the motion animation from a saved `simulation_results.mat` without re-running the ODE (different frame stride or frame rate), use [`post_processing/render_animation_from_results.m`](post_processing/render_animation_from_results.m).

## Main Outputs

Results are written under `output/<case_name>/` or an automatically timestamped run folder.

| File | Purpose |
| --- | --- |
| `simulation_results.mat` | Full `results` and `params` structures. |
| `timeseries.csv` | Scalar time series for external analysis. |
| `config.json` | Saved run configuration for reproduction. |
| `summary_plot.png` / `.pdf` / optional `.html` | Main trajectory and summary figure. |
| `flowfield_panels.png` / `.pdf` / optional `.html` | XZ/YZ flow-field panels. |
| `initial_pose.png` / `.pdf` / optional `.html` | Initial pose preview, when enabled. |
| `bacteria_3d_motion.mp4` | Optional 3D animation (`bacteria_3d_motion.avi` on Linux). |
| `bacteria_3d_motion_rerender.mp4` | Optional animation re-rendered from saved results by `render_animation_from_results`. |
| `simulation.log` | Console and file log of the run. |
| `brownian_rng_seed.txt` | Brownian RNG seed, only when Brownian dynamics are enabled. |

`results.completed` and `results.termination_reason` distinguish complete
runs from timeout-truncated runs.  A segment timeout saves the available
MAT/CSV/config/log data and returns normally; batch and post-processing
drivers exclude that partial trajectory from steady-state statistics.

## Core Files

| Path | Role |
| --- | --- |
| [`main_coupled_simulation.m`](main_coupled_simulation.m) | Main script and inline configuration entry point. |
| [`run_simulation.m`](run_simulation.m) | Simulation engine, integration loop, exports, and logging. |
| [`src/config/load_sim_config.m`](src/config/load_sim_config.m) | JSON configuration loader. |
| [`src/config/save_sim_config.m`](src/config/save_sim_config.m) | JSON configuration snapshot writer. |
| [`src/config/validate_simulation_params.m`](src/config/validate_simulation_params.m) | Pre-run parameter validation. |
| [`src/config/init_coupled_params.m`](src/config/init_coupled_params.m) | Geometry, scaling, Brownian, and model parameter setup. |
| [`src/config/init_coupled_state.m`](src/config/init_coupled_state.m) | Initial state vector construction. |
| [`src/core/ode_coupled_rhs.m`](src/core/ode_coupled_rhs.m) | Coupled elastohydrodynamic ODE right-hand side. |
| [`src/core/LaggedJacobian.m`](src/core/LaggedJacobian.m) | Lagged finite-difference Jacobian handle for the Brownian `ode15s` restarts (v8.5.1). |
| [`src/core/build_coupled_lin_sys.m`](src/core/build_coupled_lin_sys.m) | Coupled linear system assembly. |
| [`src/physics/compute_pili_curvature.m`](src/physics/compute_pili_curvature.m) | Shared node-collocated, second-order curvature/twist operator (moments and `kb0`). |
| [`src/physics/compute_pili_elastic_moments.m`](src/physics/compute_pili_elastic_moments.m) | Kirchhoff restoring moments with root boundary rules. |
| [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) | Full-coupled effective body mobility, Brownian covariance factor, and dimensionless mobility. |
| [`src/physics/sample_brownian_body_force.m`](src/physics/sample_brownian_body_force.m) | Per-interval thermal generalised body force (force realisation of the Brownian increment). |
| [`src/physics/background_flow_velocity.m`](src/physics/background_flow_velocity.m) | Shared background-flow velocity evaluator (solver and flow-field plots). |
| [`src/visualization/compute_plot_flowfield.m`](src/visualization/compute_plot_flowfield.m) | End-state flow-field evaluation at the achieved final time and figure export. |
| [`post_processing/render_flowfield_from_results.m`](post_processing/render_flowfield_from_results.m) | Recompute final, selected-time, or animated flow fields from a v8.7.0 `simulation_results.mat`. |
| [`src/visualization/render_motion_animation.m`](src/visualization/render_motion_animation.m) | Shared 4-panel motion-animation renderer (3-D and XZ overviews plus co-moving 3-D and XZ windows, legend, scale bar, a short dark-green waist patch as the roll mark, time stamp with playback speed): graphics objects are created once and only their data are updated per frame; logs per-stage timing. Co-moving ticks are never finer than 0.5 μm. |
| [`src/visualization/add_gravity_indicator_outside.m`](src/visualization/add_gravity_indicator_outside.m) | Gravity-arrow annotation shared by the preview, summary and animation figures (optional `head_scale` argument). |

## Documentation

- [`README.zh-CN.md`](README.zh-CN.md): Chinese version of this entry page.
- [`USER_GUIDE.md`](USER_GUIDE.md): practical run guide, output notes, custom geometry notes, and batch workflow notes.
- [`USER_GUIDE.zh-CN.md`](USER_GUIDE.zh-CN.md): Chinese version of the user guide.
- [`docs/theory.md`](docs/theory.md): theory, numerical model, nondimensionalization, linear system, Brownian dynamics, and validation notes.
- [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md): Chinese version of the theory document.
- [`CHANGELOG.md`](CHANGELOG.md): release history and compatibility notes.
- [`docs/index.html`](docs/index.html): project website.

## License

MIT License. See [`LICENSE`](LICENSE).
Copyright (c) 2026 Lucien.
