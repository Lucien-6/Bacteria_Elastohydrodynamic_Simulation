# BacEHD — CPC Reviewer Package

Short, self-contained verification and benchmark suite for the *Computer Physics Communications* (CPC) deposit of **BacEHD** (`v10.2.2`). Long Brownian ensembles and publication sweeps are not repeated here; this package verifies every core physical law, v10.0.0+ numerical acceleration, thermal operator, and full end-to-end workflow in seconds to a couple of minutes.

See [`PROGRAM_SUMMARY.md`](PROGRAM_SUMMARY.md) for the formatted CPC Program Summary metadata.

## Quick Run (Test Run Input)

From MATLAB, with the project root as the current folder:

```matlab
addpath('cpc')
run_cpc_tests('SkipExternal', true)   % 18 built-in checks (~22 s)
```

To also run the two external literature checks (Johnson slender-body helix propulsion and clamped cantilever relaxation in `validation/`):

```matlab
run_cpc_tests                         % 18 built-in + 2 external checks (~1-2 min)
```

To run the v10.0.0+ Block Schur / 3-colour Jacobian speedup & equivalence benchmark and the multi-worker `parfor` / GPU flow-field scaling benchmark:

```matlab
addpath('validation')
cpc_performance_benchmark               % writes output/cpc_checks/performance_benchmark_*.csv (~17 s)
cpc_parallel_and_flowfield_benchmark    % writes output/cpc_checks/parallel_scaling_benchmark.csv & flowfield_gpu_benchmark.csv
```

`run_cpc_tests` automatically adds `src/`, `validation/`, `chlamydomonas/`, and the bundled MIT-licensed S2 Sampling Toolbox (`MATLAB Add-On/Suite of functions to perform uniform sampling of a sphere`, Copyright (c) 2019 Anton Semechko) to the MATLAB path.

## What Is Checked

| Test | Reference / Invariant | Pass band |
| --- | --- | --- |
| `sphere` | Stokes sedimentation speed of a 0.4 µm sphere (`N_body = 600`) | 1% |
| `mobility_translation` | Stokes–Einstein translational diffusivity from RSM mobility | 2% |
| `mobility_rotation` | Stokes–Einstein rotational diffusivity (`N_body = 600`) | 2% |
| `capsule_speed`, `capsule_drift` | Axisymmetric rigid-body Stokes relation at 45° | `1e-6` relative |
| `capsule_spin` | Fore-aft symmetric capsule does not rotate while settling | `1e-8` rad/s |
| `uniform_flow` (4 checks) | Cell + 3 µm filament translates rigidly with a 1 µm/s uniform stream | `1e-6` to `1e-5` |
| `motor` | Signed steady rotor frequency [Hz] over 0.02–0.05 s vs `τ / [2π (2π μ (β r)² L + τ / ω_max)]` | 2% |
| `filament_sink` | Instantaneous speed ratio at 60° (3 µm filament, 30 segments) | `0.931 ± 0.02` |
| `schur_equivalence` | Body-frame Block Schur vs dense LU single-step `ode_coupled_rhs` relative difference | `1e-10` |
| `pili_thermal_pattern`, `pili_thermal_twist`, `pili_thermal_sampler` | v10.0.0 slow-mode thermal operator + balancing axial torque density (`pattern_resid`, `twist_row_resid`, SPD & finite sampler) | `1e-10` / exact |
| `end_to_end_pipeline` | Full `load_sim_config` → `validate_simulation_params` → `run_simulation` → MAT/CSV/JSON export | pass flag = 1 |
| `chlamy_3d_dry_run` | Custom spheroid + measured-gait `fit_1p00A.mat` 3-D *Chlamydomonas* `dry_run` | pass flag = 1 |
| `helix` *(external)* | Rigid-helix thrust, torque, and drag vs Johnson slender-body theory (`validation/cpc_helix_propulsive_matrix.m`) | 5% vs SBT, 2% grid conv. |
| `cantilever` *(external)* | Clamped-filament slowest bending time vs local Euler–Bernoulli formula (`validation/cpc_cantilever_relaxation.m`) | 10% vs analytic, 2% grid conv. |

## Expected Test Run Output

Reference output on MATLAB R2026a (Windows 11, PCWIN64) for `run_cpc_tests('SkipExternal', true)` (report saved to `output/cpc_package/cpc_test_report.csv`):

```text
================================================================
  BacEHD v10.2.2 — CPC Reviewer Verification Suite
  Platform: MATLAB 26.1.0.3346908 (R2026a) Update 5 on PCWIN64
================================================================
  sphere                 V / V_Stokes                 +0.998413  (expected +1 +/- 0.01)  PASS
  mobility_translation   D_t / D_Stokes               +0.998427  (expected +1 +/- 0.02)  PASS
  mobility_rotation      D_r / D_Stokes               +0.982656  (expected +1 +/- 0.02)  PASS
  capsule_speed          V(45 deg) / analytic         +1  (expected +1 +/- 1e-06)  PASS
  capsule_drift          drift / analytic             +1  (expected +1 +/- 1e-06)  PASS
  capsule_spin           angular speed (rad/s)        +8.11154e-18  (expected +0 +/- 1e-08)  PASS
  uniform_flow           V_x / U                      +1  (expected +1 +/- 1e-06)  PASS
  uniform_flow_cross     cross-flow speed / U         +2.24949e-17  (expected +0 +/- 1e-06)  PASS
  uniform_flow_spin      angular speed * L / U        +8.67229e-16  (expected +0 +/- 1e-06)  PASS
  uniform_flow_bend      angle rate * L / U           +2.65983e-12  (expected +0 +/- 1e-05)  PASS
  motor                  omega / analytic             +0.999994  (expected +1 +/- 0.02)  PASS
  filament_sink          instantaneous V / V_bare     +0.929507  (expected +0.931 +/- 0.02)  PASS
  schur_equivalence      max |dZ_schur-dZ_dense|/|dZ| +8.52171e-13  (expected +0 +/- 1e-10)  PASS
  pili_thermal_pattern   modal load pattern_resid     +1.00944e-15  (expected +0 +/- 1e-10)  PASS
  pili_thermal_twist     twist_row_resid              +5.04661e-21  (expected +0 +/- 1e-10)  PASS
  pili_thermal_sampler   SPD & finite sampler flag    +1  (expected +1 +/- 0)  PASS
  end_to_end_pipeline    JSON->run_simulation flag    +1  (expected +1 +/- 0)  PASS
  chlamy_3d_dry_run      3-D Chlamydomonas dry_run    +1  (expected +1 +/- 0)  PASS
BacEHD v10.2.2 CPC tests PASS in ~22.5 s.
```

---

Author: Lucien (lucien-6@qq.com)
Created: 2026-11-17
Updated: 2026-11-24
