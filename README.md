# Bacteria Elastohydrodynamic Simulation

English | [中文](README.zh-CN.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Created:** 2026-02-06  
> **Updated:** 2026-10-20
> **Version:** 8.6.5
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
| [`src/visualization/render_motion_animation.m`](src/visualization/render_motion_animation.m) | Shared 4-panel motion-animation renderer (3-D and XZ overviews plus co-moving 3-D and XZ windows, legend, scale bar, roll-marker stripes on the body, time stamp with playback speed): graphics objects are created once and only their data are updated per frame; logs per-stage timing. |
| [`src/visualization/add_gravity_indicator_outside.m`](src/visualization/add_gravity_indicator_outside.m) | Gravity-arrow annotation shared by the preview, summary and animation figures (optional `head_scale` argument). |

## Documentation

- [`README.zh-CN.md`](README.zh-CN.md): Chinese version of this entry page.
- [`USER_GUIDE.md`](USER_GUIDE.md): practical run guide, output notes, custom geometry notes, and batch workflow notes.
- [`USER_GUIDE.zh-CN.md`](USER_GUIDE.zh-CN.md): Chinese version of the user guide.
- [`docs/theory.md`](docs/theory.md): theory, numerical model, nondimensionalization, linear system, Brownian dynamics, and validation notes.
- [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md): Chinese version of the theory document.
- [`CHANGELOG.md`](CHANGELOG.md): release history and compatibility notes.
- [`docs/index.html`](docs/index.html): project website.

## Version Notes

The current release is **8.6.5 (2026-10-20)**, an example-case release. It adds a realistic single-polar-flagellum *E. coli* swimmer as a custom-geometry configuration ([`configs/config_ecoli_single_flagellum.json`](configs/config_ecoli_single_flagellum.json), geometry from [`geometry/generate_ecoli_geometry.m`](geometry/generate_ecoli_geometry.m)): a 2.5 µm × 0.88 µm sphero-cylinder (`N_body = 300` requested, 326 points from `capsule_body_building`), one left-handed helical filament at the tail pole (pitch 2.22 µm, helix diameter 0.41 µm, contour about 7.6 µm, filament radius 10 nm, EI = 3.5 pN·µm²) driven by a free-torsion motor (1300 pN·nm, 300 Hz), neutrally buoyant, with `T_obs = 1 s`. **The solver is unchanged, so the results of existing cases are unchanged.**

**8.6.3 (2026-10-19)** was an animation-layout refinement release. The co-moving windows c and d ceased to be cubes / squares: along every lab axis their extent became the reach of the body and all pili over all frames plus a margin, at equal scale on X, Y and Z, so the bacterium filled the window; they carried no scale bar (panel b kept its 2-μm scale bar). The body carried two dark longitudinal stripes that rotated with it, so a roll about the long axis could be seen. The four panels were arranged automatically: all guillotine arrangements were tried and the one that allowed the largest common scale of the panels was chosen (logged as, for example, `Layout H(a,H(b,V(c,d)))`), so the canvas stayed compact for elongated, flat and compact bacteria alike. `info.window_half_um` of `render_motion_animation` became a `[hx hy hz]` vector and `info.window_offset_um` was added. **Simulation results were unchanged.**

**8.6.2 (2026-10-18)** was a publication-oriented animation-layout release. Its motion animation had four panels on a fixed 16:9 canvas (2560 x 1440 logical px): **a** 3-D overview, **b** XZ overview with a 2-μm scale bar, **c** a co-moving 3-D window and **d** a co-moving XZ window, both centred on the centre of mass at constant magnification; the YZ panel is gone. Windows c and d carry ticks, tick labels and a grid in lab-frame coordinates, so the grid slides past the bacterium. A legend (body, every pilus, centre of mass, trajectory), a global time stamp with the playback speed relative to real time, panel titles flush on the upper-left corner of each image, scale bars and larger fonts / line widths were added; the gravity arrow is kept in panel a and the trajectory is drawn in panels a and b only. The video size became the fixed 16:9 canvas instead of a data-dependent size. **Simulation results are unchanged.**

**8.6.0 (2026-10-17)** was an animation-rendering performance release. Its 3-panel motion animation (3D + XZ + YZ) was rendered by the shared [`src/visualization/render_motion_animation.m`](src/visualization/render_motion_animation.m), which creates all graphics objects once and only updates their data on each frame (previously every frame ran `cla` and re-created the whole scene): a 101-frame test render dropped from about 28 s to about 13 s with unchanged video content, resolution and quality. `run_simulation.m` and `run_simulation_HPC_Linux.m` now share this implementation (platform-aware `MPEG-4` / `Motion JPEG AVI` profile), `simulation.log` gains a per-stage timing table, and [`post_processing/render_animation_from_results.m`](post_processing/render_animation_from_results.m) re-renders the animation from a saved `simulation_results.mat` without re-running the ODE. A defect in which the XZ/YZ border lines accumulated hidden duplicates on every frame is fixed. **Simulation results are unchanged.**

**8.5.2 (2026-10-11)** was a hydrodynamic-calibration release. The pilus regularisation length used to be `eps = r_pili` without calibration: for the Cortez 7/2 blob a straight rotlet-segment line has the axial-twist drag `2*pi*mu*eps^2`, only half of the no-slip-cylinder value `4*pi*mu*r^2`, and the translational mobility was 8-9 % above the slender-body theory of a cylinder. The new factor `material.hydro_radius_factor` (`cfg.hydro_radius_factor`, default `sqrt(2)`) sets `eps_pili = factor * r_pili`: the default gives the exact cylinder twist drag and translation within about 3 %, `sqrt(e)` gives exact translation, and `1` reproduces 8.5.1. The bending stiffness and the root gap still use the physical radius. **Default results of pili-bearing cases change** (pilus hydrodynamics by a few per cent, twist drag doubles; the one-pilus sink distance changes by 0.13 %); pilus-free runs are unaffected. The release also removes an unused assignment in `init_coupled_params.m`. Release **8.5.1 (2026-10-10)** was a numerical-performance and maintenance release: Brownian runs re-use the last exact finite-difference Jacobian across `ode15s` restarts (`src/core/LaggedJacobian.m`, runtime option `cfg.jacobian_reuse`, default on), cutting the ODE cost of a Brownian node by about 2-3.7x for multi-pilus cases (about 1.4x for one pilus) while the trajectories stay within the integration tolerance of v8.5.0 (`cfg.jacobian_reuse = false` reproduces v8.5.0 bit-for-bit; deterministic runs are unchanged). The release also cleans up the MATLAB code-analyzer messages of `run_simulation.m` / `run_simulation_HPC_Linux.m`, fixes `verify_migration.m` (it now lists `LaggedJacobian.m` and uses `datetime`), adds Markdown lint configuration (`.markdownlint.json`, `.markdownlint-cli2.jsonc`) and repairs the CHANGELOG structure. Release **8.5.0 (2026-09-04)** changes how body-only Brownian motion is applied: each Brownian interval's thermal increment is realised as a constant generalised force `f_th = M_eff \ (dU/dt_B)` acting on the body through the coupled ODE (`sample_brownian_body_force.m`, `ode_coupled_rhs` Step 8a) instead of an instantaneous displacement kick followed by root re-projection. Pili now follow the body's thermal motion continuously, so the spurious curvature kink between segments 1 and 2 that every kick used to inject — and the stiff high-wavenumber relaxation it forced onto `ode15s` — is gone; the Ito drift `kBT div(M)` is generated implicitly by the mobility gradient along the trajectory, so the explicit finite-difference drift (12 extra coupled solves per node) was removed. Brownian trajectories differ from earlier versions even for the same seed; deterministic runs are unchanged. The type-2 parabolic pilus template is also anchored exactly at the origin. Release **8.4.0 (2026-09-04)** corrected the `free_torsion` root twist equation: the root $\psi_1$ row is the whole-filament axial moment balance and now carries only the boundary (motor) torque, whereas the interior Kirchhoff twist moment used to be added as well, making the motor balance the drag of the first segment only and inflating the steady rotor speed by the segment count $N$. Release **8.3.0 (2026-09-04)** was a discretisation-accuracy and geometry-consistency update: the Kirchhoff curvature stencil is now node-collocated and second-order accurate (shared operator `compute_pili_curvature.m`, also used for `kb0`), and the 14 attachment-site maps in `pili_building.m` are proper rotations (sites 1, 2, 5, 7, 10, 11, 14 previously reflected chiral pili) with a right-handed `deflection_x`; planar pili are placed identically to before. The preceding release, **8.2.0 (2026-09-04)**, was a physics-correctness update from the full mathematical review: the flagellar motor now follows a signed (asymmetric) linear torque-speed law in the rotor speed *relative to the cell body* (full axial rate `dpsi_1 + cos(theta_1) dphi_1 - Omega_body . d3`), with a documented sign convention (positive torque = right-handed about the outward root tangent); the flow-field post-processor solves force densities from the velocity relative to the ambient flow and adds the background flow to the plotted field; `validate_simulation_params` rejects missing/invalid `shape_params` for pilus types 2-5 with a clear message; and the user guides document the `1/L_pili` and `EI/L_pili` units of all curvature and moment inputs.

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
