# Bacteria Elastohydrodynamic Simulation

English | [中文](README.zh-CN.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Created:** 2026-02-06  
> **Updated:** 2026-09-04
> **Version:** 8.5.0
> **Runtime:** MATLAB-only. No Python or pip dependencies are required.

This project simulates the low-Reynolds-number settling dynamics of a rigid bacterial body coupled to one or more flexible pili or flagella. The body is represented by regularized Stokeslets on a surface point cloud, while each pilus is represented as an inextensible Kirchhoff elastic rod coupled to the surrounding Stokes flow through regularized Stokeslet segments.

## Key Capabilities

- Rigid capsule or imported custom body geometry.
- Parametric or imported custom pili centerlines.
- Coupled RSM/RSS hydrodynamics with non-local segment-integral twist mobility, reciprocal Schur-complement coupling, and physical standard-rotlet scaling.
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

Use [`configs/config_custom_example.json`](configs/config_custom_example.json) for imported body and pili geometry, and [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json) for free-torsion motor examples.

For machine-specific cases, keep copied JSON files under a local `configs/` path (for example `configs/my_case.json`). The repository tracks example JSON files under `configs/`; only `output/` is listed in `.gitignore`.

### Option C: Programmatic Use

```matlab
cfg = load_sim_config('configs/config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

Batch sweep **driver scripts** live in [`batch_sweep/`](batch_sweep/). Parameter-sweep CSV/MAT products are written under `batch sweep/` (with a space) when those drivers run. To batch-run every `*.json` in a folder (optional parallel, duplicate `case_name` skip, real-time `batch_summary_*.log`), use [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m). Multi-case post-processing is available in [`post_processing/postprocess_cases.m`](post_processing/postprocess_cases.m): 3-D translational MSD fit (\(\mathrm{MSD}=6D\tau+V^2\tau^2\) on \(0<\tau\le T/2\)) and XZ-plane rotational diffusion from MSAD (\(\mathrm{MSAD}=2D_\mathrm{rot}\tau\)), exported to Excel plus a PNG/PDF/(optional HTML) summary figure.

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
| `bacteria_3d_motion.mp4` | Optional 3D animation. |
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
| [`src/core/build_coupled_lin_sys.m`](src/core/build_coupled_lin_sys.m) | Coupled linear system assembly. |
| [`src/physics/compute_pili_curvature.m`](src/physics/compute_pili_curvature.m) | Shared node-collocated, second-order curvature/twist operator (moments and `kb0`). |
| [`src/physics/compute_pili_elastic_moments.m`](src/physics/compute_pili_elastic_moments.m) | Kirchhoff restoring moments with root boundary rules. |
| [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) | Full-coupled effective body mobility, Brownian covariance factor, and dimensionless mobility. |
| [`src/physics/sample_brownian_body_force.m`](src/physics/sample_brownian_body_force.m) | Per-interval thermal generalised body force (force realisation of the Brownian increment). |
| [`src/physics/background_flow_velocity.m`](src/physics/background_flow_velocity.m) | Shared background-flow velocity evaluator (solver and flow-field plots). |
| [`src/visualization/compute_plot_flowfield.m`](src/visualization/compute_plot_flowfield.m) | End-state flow-field evaluation at the achieved final time and figure export. |

## Documentation

- [`README.zh-CN.md`](README.zh-CN.md): Chinese version of this entry page.
- [`USER_GUIDE.md`](USER_GUIDE.md): practical run guide, output notes, custom geometry notes, and batch workflow notes.
- [`USER_GUIDE.zh-CN.md`](USER_GUIDE.zh-CN.md): Chinese version of the user guide.
- [`docs/theory.md`](docs/theory.md): theory, numerical model, nondimensionalization, linear system, Brownian dynamics, and validation notes.
- [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md): Chinese version of the theory document.
- [`CHANGELOG.md`](CHANGELOG.md): release history and compatibility notes.
- [`docs/index.html`](docs/index.html): project website.

## Version Notes

The current release is **8.5.0 (2026-09-04)**, which changes how body-only Brownian motion is applied: each Brownian interval's thermal increment is realised as a constant generalised force `f_th = M_eff \ (dU/dt_B)` acting on the body through the coupled ODE (`sample_brownian_body_force.m`, `ode_coupled_rhs` Step 8a) instead of an instantaneous displacement kick followed by root re-projection. Pili now follow the body's thermal motion continuously, so the spurious curvature kink between segments 1 and 2 that every kick used to inject — and the stiff high-wavenumber relaxation it forced onto `ode15s` — is gone; the Ito drift `kBT div(M)` is generated implicitly by the mobility gradient along the trajectory, so the explicit finite-difference drift (12 extra coupled solves per node) was removed. Brownian trajectories differ from earlier versions even for the same seed; deterministic runs are unchanged. The type-2 parabolic pilus template is also anchored exactly at the origin. Release **8.4.0 (2026-09-04)** corrected the `free_torsion` root twist equation: the root $\psi_1$ row is the whole-filament axial moment balance and now carries only the boundary (motor) torque, whereas the interior Kirchhoff twist moment used to be added as well, making the motor balance the drag of the first segment only and inflating the steady rotor speed by the segment count $N$. Release **8.3.0 (2026-09-04)** was a discretisation-accuracy and geometry-consistency update: the Kirchhoff curvature stencil is now node-collocated and second-order accurate (shared operator `compute_pili_curvature.m`, also used for `kb0`), and the 14 attachment-site maps in `pili_building.m` are proper rotations (sites 1, 2, 5, 7, 10, 11, 14 previously reflected chiral pili) with a right-handed `deflection_x`; planar pili are placed identically to before. The preceding release, **8.2.0 (2026-09-04)**, was a physics-correctness update from the full mathematical review: the flagellar motor now follows a signed (asymmetric) linear torque-speed law in the rotor speed *relative to the cell body* (full axial rate `dpsi_1 + cos(theta_1) dphi_1 - Omega_body . d3`), with a documented sign convention (positive torque = right-handed about the outward root tangent); the flow-field post-processor solves force densities from the velocity relative to the ambient flow and adds the background flow to the plotted field; `validate_simulation_params` rejects missing/invalid `shape_params` for pilus types 2-5 with a clear message; and the user guides document the `1/L_pili` and `EI/L_pili` units of all curvature and moment inputs.

**8.1.0 (2026-07-25)** was a batch post-processing feature update: [`post_processing/postprocess_cases.m`](post_processing/postprocess_cases.m) fits 3-D translational diffusion on \(0<\tau\le T/2\) and adds XZ-plane rotational diffusion \(D_\mathrm{rot}\) from MSAD (prefer body quaternion; fall back to `tilt_rad`). The hydrodynamic / elastic / Brownian simulation core is unchanged from 8.0.1.

**8.0.1 (2026-07-21)** was a numerical-diagnostics hardening patch on top of 8.0.0. Dense weighted-mobility dissipation checks use full `eig` (iterative `eigs` only for sparse matrices, with fallback); isolated `rcond < 1e-14` warns as `NearSingularMatrix` instead of aborting; `SingularMatrix` requires the joint gate of a large backward residual and near-singular `rcond`; health checks probe solve quality through `check_linear_solve`. Diagnostic snapshots expose `rcond_near_singular_threshold`.

**8.0.0 (2026-07-13)** was the major mathematical-physics correctness and runtime-contract release: weighted-reciprocal translation/twist hydrodynamics for unequal segment lengths, equal-chord canonical custom-pilus geometry, objective full-director initialisation, dimensionless body-only Brownian SPD repair with remainder kicks and isolated RNG, strict incompressible `shear_z` orthogonality, timeout-safe partial returns, and auditable numerical guards via `simulation.log` / `results.diagnostics`.

**Breaking changes retained from 8.0.0:** `Pili_Matrix` must have exactly eight rows; custom or named `bending_stiffness` fields are rejected in favour of global `material.E` and `r_pili`; non-orthogonal `shear_z` configurations are rejected; and a timed-out run returns `results.completed=false` with `termination_reason='timeout'` after saving its partial trajectory.

**7.10.0 (2026-07-03)** added Linux HPC cluster execution and fixed two cross-platform batch-driver bugs. It introduced the maintained **[`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m)** / **[`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m)** fork with platform-aware video export, Linux-native resource probing, node-local `parpool` job storage, heavy-output sampling, and live progress heartbeats.

**7.9.6 (2026-07-03)** was a console-output and error-resilience correctness release: a review of `run_simulation.m`'s logging pipeline found that `dual_log`'s console-echo flag had been inverted since v7.9.3, so every normal (non-batch) run echoed **nothing** through `logf` while batch/parallel workers echoed everything — the opposite of the documented `cfg.silent` contract; `simulation.log` itself was unaffected. It also found that custom-geometry (`geometry_mode='custom'`) runs crashed deterministically at the (non-optional) `[PLOT] Generating summary figure` step because `trisurf(ax, F, x, y, z, ...)` does not support an axes handle as its first argument the way `surf`/`plot3` do. Both were fixed, plus a run-wide `try/catch` so `simulation.log` always records a `[FATAL]` block instead of truncating silently on error.

**7.9.5 (2026-07-03)** extended silent-mode coverage to **custom geometry**: `custom_body_building` and `custom_pili_building` accepted the injected `logf` handle (so `[Custom Body]` / `[Custom Pili]` diagnostics respect `cfg.silent`), and the summary figure used `Visible='off'` under silent (matching the animation figure). **`cfg.silent`** is a **runtime-only** batch flag (not stored in JSON): when true, console echo of `logf`, direct `fprintf`, ODE/animation progress bars, and warnings are suppressed for the run; full text still goes to `simulation.log`; errors are never suppressed. Parallel batch drivers (e.g. `batch_repeat_brownian`) set `cfg.silent = true` on workers.

**7.9.4 (2026-07-02)** closed the remaining `cfg.silent` console-output gaps (direct `fprintf`, `odetpbar`/`render_progress_bar` progress bars, ~24 `warning` calls across 16 `src/` files). 7.9.3 added `cfg.silent` itself plus incremental batch-summary writes; 7.9.2 was a numerical-robustness pass; 7.9.1 was a post-review bug-fix release; 7.9.0 was a parameter-validation deep refactor; 7.8.0 corrected the reference-curvature `kb0` discretisation and the constant motor-torque nondimensionalisation. For the full history, see [`CHANGELOG.md`](CHANGELOG.md).

## License

MIT License. See [`LICENSE`](LICENSE).
Copyright (c) 2026 Lucien.
