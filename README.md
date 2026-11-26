# BacEHD — Bacteria Elastohydrodynamic Simulation

English | [中文](README.zh-CN.md)

> **Program title:** BacEHD (see [`cpc/PROGRAM_SUMMARY.md`](cpc/PROGRAM_SUMMARY.md) and [`CITATION.cff`](CITATION.cff))  
> **Authors:** Shangping Lu (Lucien, [lucien-6@qq.com](mailto:lucien-6@qq.com)), Kun Zhao ([kzhao@uestc.edu.cn](mailto:kzhao@uestc.edu.cn))  
> **Repository:** [GitHub](https://github.com/Lucien-6/Bacteria_Elastohydrodynamic_Simulation)  
> **Created:** 2026-02-06  
> **Updated:** 2026-11-25  
> **Version:** 10.2.2  
> **Runtime:** Base MATLAB R2022a+ (no Statistics Toolbox, MEX, or Python required for simulations; optional Python plotting scripts in `post_processing/`).

This project simulates the low-Reynolds-number dynamics of a rigid cell body coupled to one or more flexible pili or flagella: settling under gravity, swimming driven by a motor torque or by active moments, motion in a background flow, and Brownian motion. The body is represented by regularized Stokeslets on a surface point cloud, while each pilus is represented as an inextensible Kirchhoff elastic rod coupled to the surrounding Stokes flow through regularized Stokeslet segments.

## Key Capabilities

- Rigid capsule or imported custom body geometry.
- Parametric or imported custom pili centerlines.
- Coupled RSM/RSS hydrodynamics with non-local segment-integral twist mobility, reciprocal Schur-complement coupling, and physical standard-rotlet scaling, and a calibrated pilus hydrodynamic radius (`material.hydro_radius_factor`).
- Kirchhoff elastic rods with intrinsic curvature and active moment options.
- Optional body-only Brownian dynamics: the 6-DOF body covariance uses the full coupled body-pili effective mobility with dimensionless SPD repair and an isolated reproducible random stream; each interval's increment is realised as a constant thermal generalised force acting through the coupled ODE (implicit Ito drift, no root kink). Pili internal elastic modes receive no thermal noise of their own unless the next option is on.
- Optional thermal fluctuations of the pili (9.0.0+, hardened in 10.0.0, `brownian.pili_thermal`): the slow bending modes of the pili receive thermal forces sampled jointly with the body from the coupled mobility and realised with the constrained probe's centreline force density and balancing axial torque density (`brown_pili_tau`), so these modes reach equipartition with negligible twist leakage. Supports time-scale cutoff (`brownian.pili_tau_cut_s`), near-degenerate pair protection (`brownian.pili_protect_pairs`), and optional per-pilus RNG streams (`brownian.pili_rng_seeds`). Twist and higher bending modes stay deterministic. For passive, naturally straight pili with clamped roots. See [`docs/theory.md`](docs/theory.md) §10.4.
- Body-frame Block Schur hydrodynamic precomputation (`build_body_precomp`, `build_pili_mobility_blocks`) and fast frozen-hydrodynamics 3-colour Jacobian (`ode_coupled_jacobian`, paired with `LaggedJacobian`) for large speedups on multi-pilus and Brownian runs (10.0.0).
- Automated CPC reviewer verification suite ([`cpc/run_cpc_tests.m`](cpc/run_cpc_tests.m), [`cpc/PROGRAM_SUMMARY.md`](cpc/PROGRAM_SUMMARY.md)), external slender-body / cantilever checks ([`validation/`](validation/)), and performance/parallel/GPU scaling benchmarks ([`validation/cpc_performance_benchmark.m`](validation/cpc_performance_benchmark.m), [`validation/cpc_parallel_and_flowfield_benchmark.m`](validation/cpc_parallel_and_flowfield_benchmark.m)).
- Optional background flow: off, uniform, or finite-thickness shear.
- Free-torsion root coupling and constant motor torque for motor-driven flagella.
- Biflagellate *Chlamydomonas* 3-D swimmer: measured-gait dual-k 1ω+2ω in-plane drive combined with 17° chiral beat-plane tilt and zero-twist proximal 1ω+2ω out-of-plane bending, reproducing 1–2 Hz spin and left-handed superhelical trajectories. See [`chlamydomonas/README.md`](chlamydomonas/README.md) and [`chlamydomonas/cases/run_3d_helix_2s.m`](chlamydomonas/cases/run_3d_helix_2s.m).
- Inline MATLAB configuration, JSON configuration, and programmatic batch execution.
- Publication-oriented exports: MAT, CSV, log, PNG, vector PDF, optional interactive HTML, and MP4 animation.
- Accelerated 6-panel 3-D/2-D streamline flow-field reconstruction & animation (`render_flowfield_from_results`, 10.2.0): direct Block Schur force extraction, batched 3-plane FP64 GPU/BLAS3 velocity kernels (`use_gpu = 'auto'`), exact analytical capsule cross-sections, optimised thread-pool parallel streamlines (`use_parallel = 'auto'`), and a slimmed ~1.6–4.0 MiB/frame cache (domain-dependent).
- Optional **parallel** batch drivers ([`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m), [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m)): `parfor` with automatic pool sizing (Parallel Computing Toolbox); serial fallback when the toolbox or pool is unavailable.
- **Linux HPC cluster execution:** [`run_simulation.m`](run_simulation.m) natively supports headless/Linux operation (platform-aware MP4/AVI video export via `render_motion_animation.m`, `.avi` detection, `cfg.silent`, and `DataQueue` progress heartbeats); [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) is retained as a thin forwarding alias for backward compatibility with [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m) (Linux-native CPU/memory probing, node-local `parpool` job storage, and heavy-output sampling).

## Directory Structure

```text
<project_root>/
├── main_coupled_simulation.m    # Interactive inline/JSON entry script
├── run_simulation.m             # Unified simulation engine (Windows & Linux HPC)
├── src/                         # Species-agnostic solver modules
│   ├── config/                  #   JSON loader/saver, parameter validator, initialiser
│   ├── geometry/                #   Capsule/custom body & equal-chord pili builders + S2 cache
│   ├── core/                    #   Block Schur linear system, RHS, 3-colour Jacobian, diagnostics
│   ├── physics/                 #   Kirchhoff elasticity, active moments, two-tier Brownian dynamics
│   ├── math/                    #   Quaternions, directors, Euler angles, adaptive SO(3) basis
│   ├── visualization/           #   4-panel motion animation & GPU/ThreadPool 6-panel flow field
│   └── utils/                   #   Version, state packing, ETA/progress helpers
├── configs/                     # Example & species JSON configurations
├── geometry/                    # Example, E. coli & Chlamydomonas point-cloud CSVs + generators
├── chlamydomonas/               # Pluggable 3-D Chlamydomonas inverse/forward swimmer module
├── cpc/                         # CPC Program Summary & 18+2 automated reviewer test suite
├── validation/                  # External SBT/cantilever checks & operator/parfor/GPU benchmarks
├── batch_sweep/                 # Multi-worker parfor & Linux HPC ensemble/sweep drivers
├── post_processing/             # Offline flow-field/animation redraw, MSD/MSAD & optional plot_*.py
├── docs/                        # Theory & numerical model documentation (EN & zh-CN)
└── MATLAB Add-On/               # Bundled MIT-licensed S2 Sampling Toolbox (A. Semechko)
```

## Requirements

- MATLAB R2022a or later (not compatible with GNU Octave); developed and tested with R2026a on Windows 11 and R2023a on Linux HPC. Earlier releases are not supported because runs write `config.json` with `jsonencode(..., 'PrettyPrint', true)` (R2021a) and the default flow-field figure uses `clim` (R2022a). Flow-field colour scaling uses pure base-MATLAB percentile interpolation (no Statistics and Machine Learning Toolbox dependency).
- PNG and PDF figures are exported with `exportgraphics`; optional interactive HTML figure export requires MATLAB R2026a or later.
- Capsule body-surface generation uses Anton Semechko's MIT-licensed **S2 Sampling Toolbox**, which is **bundled in this repository** under `MATLAB Add-On/Suite of functions to perform uniform sampling of a sphere/` (added to the path automatically by `cpc/run_cpc_tests.m` and validation scripts).
- Core MATLAB simulations require no Python, Conda, or pip environment; optional publication re-plotting scripts (`post_processing/plot_*.py`) use standard Python plotting packages listed in [`requirements.txt`](requirements.txt).
- **Parallel Computing Toolbox** (and an optional CUDA GPU) are optional: used automatically for parallel batch sweeps (`batch_sweep/`) and accelerated multi-frame flow-field animations (`use_gpu = 'auto'`, `use_parallel = 'auto'`), with seamless serial CPU fallback when unavailable.

## Quick Start

Before running a new case, confirm that the project root is the current MATLAB folder or is on the MATLAB path.

### Reviewer Verification & Benchmarks (CPC Test Run)

```matlab
addpath('cpc'); addpath('validation');
run_cpc_tests('SkipExternal', true)      % 18 built-in checks (~22.5 s)
run_cpc_tests                            % 18 built-in + 2 external checks (~1-2 min)
cpc_performance_benchmark                % Block Schur & 3-colour Jacobian speedup (~17 s)
cpc_parallel_and_flowfield_benchmark     % parfor strong scaling & GPU flow-field benchmark
```

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

Use [`configs/config_custom_example.json`](configs/config_custom_example.json) for imported body and pili geometry, and [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json) for free-torsion motor examples.  A realistic single-flagellum *E. coli* swimmer (2.5 µm × 0.88 µm body, left-handed helical filament with pitch 2.22 µm, motor torque 1300 pN·nm) is provided in [`configs/config_ecoli_single_flagellum.json`](configs/config_ecoli_single_flagellum.json).  A realistic biflagellate *Chlamydomonas* 3-D helical swimmer is started from the MATLAB desktop: open [`chlamydomonas/cases/run_3d_helix_2s.m`](chlamydomonas/cases/run_3d_helix_2s.m) (3-D helical beat, 2.0 s, left-handed superhelix) or [`chlamydomonas/cases/run_planar_sym50.m`](chlamydomonas/cases/run_planar_sym50.m) (coplanar symmetric reference) and press Run; see [`chlamydomonas/README.md`](chlamydomonas/README.md) for full preset options.

For machine-specific cases, keep copied JSON files under a local `configs/` path (for example `configs/my_case.json`). The repository tracks example JSON files under `configs/`, and the folder is not ignored by git, so a private case file shows up as untracked. Run products go to `output/`, which is ignored.

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
| `flowfield_panels.png` / `.pdf` / optional `.html` | XZ/YZ/XY flow-field panels (2×3). |
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
| [`src/core/ode_coupled_jacobian.m`](src/core/ode_coupled_jacobian.m) | Fast frozen-hydrodynamics 3-colour Jacobian for `ode15s` (10.0.0). |
| [`src/core/LaggedJacobian.m`](src/core/LaggedJacobian.m) | Lagged Jacobian handle for Brownian `ode15s` restarts, backed by `ode_coupled_jacobian` with automatic FD fallback. |
| [`src/core/build_body_precomp.m`](src/core/build_body_precomp.m) | Body-frame rigid-body RSM Cholesky precomputation for Block Schur acceleration (10.0.0). |
| [`src/core/build_pili_mobility_blocks.m`](src/core/build_pili_mobility_blocks.m) | Body-frame non-BB hydrodynamic mobility blocks for the Block Schur solver (10.0.0). |
| [`src/core/build_coupled_lin_sys.m`](src/core/build_coupled_lin_sys.m) | Coupled linear system assembly (Block Schur fast path + dense LU path). |
| [`src/core/build_force_balance.m`](src/core/build_force_balance.m) | Force/moment-balance operator of the coupled system, split from the assembly (9.0.0). |
| [`src/physics/compute_pili_curvature.m`](src/physics/compute_pili_curvature.m) | Shared node-collocated, second-order curvature/twist operator (moments and `kb0`). |
| [`src/physics/compute_pili_elastic_moments.m`](src/physics/compute_pili_elastic_moments.m) | Kirchhoff restoring moments with root boundary rules. |
| [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) | Full-coupled effective body mobility, Brownian covariance factor, and dimensionless mobility. |
| [`src/physics/sample_brownian_body_force.m`](src/physics/sample_brownian_body_force.m) | Per-interval thermal generalised body force (force realisation of the Brownian increment). |
| [`src/physics/init_pili_brownian_modes.m`](src/physics/init_pili_brownian_modes.m) | Slow elastohydrodynamic bending modes of the pili (orders or `tau_cut_s` selection with optional pair protection, 9.0.0–10.0.0). |
| [`src/physics/compute_effective_coupled_brownian.m`](src/physics/compute_effective_coupled_brownian.m) | Joint body-mode mobility, centreline force & balancing axial torque patterns (`P_load`/`P_tau`), and matrix-free held-load drift (9.0.0–10.0.0). |
| [`src/physics/sample_brownian_coupled_force.m`](src/physics/sample_brownian_coupled_force.m) | Per-interval thermal forces and balancing axial torques on the body and slow pilus modes (9.0.0–10.0.0). |
| [`src/physics/pili_thermal_rhs.m`](src/physics/pili_thermal_rhs.m) | Right-hand side of the coupled system for thermal pilus force and axial torque loads (9.0.0–10.0.0). |
| [`src/physics/background_flow_velocity.m`](src/physics/background_flow_velocity.m) | Shared background-flow velocity evaluator (solver and flow-field plots). |
| [`src/visualization/compute_plot_flowfield.m`](src/visualization/compute_plot_flowfield.m) | End-state flow-field evaluation at the achieved final time and figure export. |
| [`post_processing/render_flowfield_from_results.m`](post_processing/render_flowfield_from_results.m) | Recompute final, selected-time, or animated flow fields from a v8.7.0 `simulation_results.mat`. |
| [`src/visualization/render_motion_animation.m`](src/visualization/render_motion_animation.m) | Shared 4-panel motion-animation renderer (3-D and XZ overviews plus co-moving 3-D and XZ windows, legend, scale bar, a short dark-green waist patch as the roll mark, time stamp with playback speed): graphics objects are created once and only their data are updated per frame; logs per-stage timing. Co-moving ticks are never finer than 0.5 μm. |
| [`src/visualization/add_gravity_indicator_outside.m`](src/visualization/add_gravity_indicator_outside.m) | Gravity-arrow annotation shared by the preview, summary and animation figures (optional `head_scale` argument). |
| [`chlamydomonas/chlamy_run_preset.m`](chlamydomonas/chlamy_run_preset.m) | Authoritative *Chlamydomonas* 3-D and planar preset entry point (see [`chlamydomonas/README.md`](chlamydomonas/README.md)). |
| [`chlamydomonas/core/chlamy_run_case.m`](chlamydomonas/core/chlamy_run_case.m) | General *Chlamydomonas* runner: measured-gait table plus 3-D out-of-plane drive and chiral tilt overrides. |
| [`chlamydomonas/core/chlamy_swim_analysis.m`](chlamydomonas/core/chlamy_swim_analysis.m) | Unwrapped SO(3) spin, multi-turn helix kinematics, body clearance, and 3-D orbit diagnostics. |

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
