# Coupled Body-Pili Elastohydrodynamic Simulation

English | [中文](README.zh-CN.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Created:** 2026-02-06  
> **Updated:** 2026-05-07  
> **Version:** 7.4.2  
> **Runtime:** MATLAB-only. No Python or pip dependencies are required.

This project simulates the low-Reynolds-number settling dynamics of a rigid bacterial body coupled to one or more flexible pili or flagella. The body is represented by regularized Stokeslets on a surface point cloud, while each pilus is represented as an inextensible Kirchhoff elastic rod coupled to the surrounding Stokes flow through regularized Stokeslet segments.

## Key Capabilities

- Rigid capsule or imported custom body geometry.
- Parametric or imported custom pili centerlines.
- Coupled RSM/RSS hydrodynamics with rotlet-segment twist coupling.
- Kirchhoff elastic rods with intrinsic curvature and active moment options.
- Optional body Brownian dynamics using the full coupled effective mobility tensor.
- Optional background flow: off, uniform, or finite-thickness shear.
- Free-torsion root coupling and constant motor torque for motor-driven flagella.
- Inline MATLAB configuration, JSON configuration, and programmatic batch execution.
- Publication-oriented exports: MAT, CSV, log, PNG, vector PDF, optional interactive HTML, and MP4 animation.
- Optional **parallel** Brownian repeat batches ([`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m)): `parfor` with automatic pool sizing (Parallel Computing Toolbox); serial fallback when the toolbox or pool is unavailable.

## Requirements

- MATLAB R2020b or later is recommended.
- PNG and PDF figure export use `exportgraphics` and work from MATLAB R2020a onward.
- Optional interactive HTML figure export requires MATLAB R2026a or later.
- Capsule body-surface generation requires the MATLAB Add-On / File Exchange **S2 Sampling Toolbox** by Anton Semechko. Install it so `ParticleSampleSphere`, `SpiralSampleSphere`, `IcosahedronMesh`, `QuadCubeMesh`, and `SubdivideSphericalMesh` are on the MATLAB path.
- Python / Conda / pip are not required. See [`requirements.txt`](requirements.txt).
- **Parallel Computing Toolbox** is optional: only needed for the parallel branch of [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m) (`use_parallel = true`, default).

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

1. Copy [`config_example.json`](config_example.json) to a local configuration path, for example `configs/my_case.json`.
2. Set `config_file` in [`main_coupled_simulation.m`](main_coupled_simulation.m):

```matlab
config_file = 'configs/my_case.json';
main_coupled_simulation
```

Use [`config_custom_example.json`](config_custom_example.json) for imported body and pili geometry, and [`config_motor_flagella_example.json`](config_motor_flagella_example.json) for free-torsion motor examples.

For machine-specific cases, keep copied JSON files under a local `configs/` folder. The repository ignores that folder by default, so local paths and exploratory parameter sets are not accidentally committed.

### Option C: Programmatic Use

```matlab
cfg = load_sim_config('config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

Batch drivers are available under [`batch_sweep/`](batch_sweep/), and multi-case post-processing is available in [`post_processing/postprocess_cases.m`](post_processing/postprocess_cases.m).

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
| [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) | Full-system effective Brownian mobility and diffusion. |
| [`src/visualization/compute_plot_flowfield.m`](src/visualization/compute_plot_flowfield.m) | End-state flow-field evaluation and figure export. |

## Documentation

- [`README.zh-CN.md`](README.zh-CN.md): Chinese version of this entry page.
- [`USER_GUIDE.md`](USER_GUIDE.md): practical run guide, output notes, custom geometry notes, and batch workflow notes.
- [`USER_GUIDE.zh-CN.md`](USER_GUIDE.zh-CN.md): Chinese version of the user guide.
- [`docs/theory.md`](docs/theory.md): theory, numerical model, nondimensionalization, linear system, Brownian dynamics, and validation notes.
- [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md): Chinese version of the theory document.
- [`CHANGELOG.md`](CHANGELOG.md): release history and compatibility notes.
- [`docs/index.html`](docs/index.html): project website.

## Version Notes

The current release is **7.4.2 (2026-05-07)**. It corrects and completes inline documentation across multiple source files — output unit convention, `Pili_Matrix` row semantics, matrix-layout block annotations (`BP = RSM`, `PB = RSS`), Euler angle convention, adaptive-basis reflection vs. rotation, step-numbering cross-references, and preferred MAT variable-name lists. All changes are comment-only; **simulation outputs are bitwise identical vs. v7.4.1**. **v7.4.1** fixed the Brownian `IllConditioned` diagnostic by evaluating `rcond` on the dimensionless mobility `M_eff_dimless`. **v7.4.0** added optional `parfor` parallel scheduling to **`batch_repeat_brownian.m`**. For the full history, see [`CHANGELOG.md`](CHANGELOG.md).

## License

MIT License. See [`LICENSE`](LICENSE).
Copyright (c) 2026 Lucien.
