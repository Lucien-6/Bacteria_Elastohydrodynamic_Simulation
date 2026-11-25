# User Guide — BacEHD (Bacteria Elastohydrodynamic Simulation)

English | [中文](USER_GUIDE.zh-CN.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Version:** 10.2.1
> **Updated:** 2026-11-24

This guide is for running simulations, choosing the right configuration path, finding outputs, and post-processing completed cases. For model equations and derivations, see [`docs/theory.md`](docs/theory.md). For release history, see [`CHANGELOG.md`](CHANGELOG.md).

## 1. Before You Run

### Requirements

- MATLAB R2022a or later. The project is developed and tested with R2026a on Windows 11 and R2023a on Linux HPC. Earlier releases are not supported: every run writes `config.json` with `jsonencode(..., 'PrettyPrint', true)` (R2021a), and the flow-field figure, which is on by default, uses `clim` (R2022a). Flow-field percentile scaling uses pure base MATLAB (no Statistics Toolbox dependency).
- PNG and PDF exports use `exportgraphics`; optional interactive HTML exports require MATLAB R2026a or later.
- Capsule body surface generation requires Anton Semechko's **S2 Sampling Toolbox** from MATLAB Add-On Explorer / File Exchange.
- No Python, Conda, or pip environment is required. [`requirements.txt`](requirements.txt) documents this explicitly.
- **Parallel Computing Toolbox** is optional: needed for the default parallel branches of [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m) and [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m) (`use_parallel = true`).
- **Running on a Linux HPC cluster (since 7.10.0):** use [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) / [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m) instead of the originals — see §5.3.

### Recommended Working Folder

Open MATLAB from the repository root or run:

```matlab
cd('path/to/project_root')
```

This keeps relative paths for JSON configs, geometry files, and outputs predictable.

## 2. Choose a Run Mode

| Goal | Recommended path |
| --- | --- |
| Quickly test or modify the default demo | Edit [`main_coupled_simulation.m`](main_coupled_simulation.m). |
| Reproduce or share a case | Use a JSON file copied from [`configs/config_example.json`](configs/config_example.json). |
| Import custom body and pili geometry | Use [`configs/config_custom_example.json`](configs/config_custom_example.json). |
| Test a free-torsion motor-driven flagellum | Use [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json). |
| Simulate a realistic *E. coli* (one polar flagellum, motor-driven) | Use [`configs/config_ecoli_single_flagellum.json`](configs/config_ecoli_single_flagellum.json); regenerate its geometry with [`geometry/generate_ecoli_geometry.m`](geometry/generate_ecoli_geometry.m). |
| Simulate a biflagellate *Chlamydomonas* in 3-D | Open [`chlamydomonas/cases/run_3d_helix_2s.m`](chlamydomonas/cases/run_3d_helix_2s.m) or [`chlamydomonas/cases/run_planar_sym50.m`](chlamydomonas/cases/run_planar_sym50.m) and press Run (see [`chlamydomonas/README.md`](chlamydomonas/README.md)). |
| Run many cases without prompts | Call [`run_simulation.m`](run_simulation.m) with `cfg.auto_confirm = true`. |

### Chlamydomonas biflagellate

The geometry is [`configs/config_chlamy_biflagellate_gz.json`](configs/config_chlamy_biflagellate_gz.json): a 5 µm by 4 µm spheroid, two clamped flagella, and a 13 degree angle between each root tangent and the cell axis. The JSON traveling-wave moment is a placeholder replaced by `chlamydomonas/data/fit_1p00A.mat` (dual-k 1ω+2ω measured-gait fit at 57.75 Hz) plus optional 3-D out-of-plane bending (`chlamydomonas/core/chlamy_make_perp_envelope.m`). See [`chlamydomonas/README.md`](chlamydomonas/README.md) for the full physics, parameters, and benchmarks.

```matlab
run('chlamydomonas/cases/run_3d_helix_2s.m')    % 3-D helical swimming, 2.0 s, with animation
run('chlamydomonas/cases/run_planar_sym50.m')   % Coplanar symmetric reference, 50 beats, with animation
```

`tilt_deg`, `drive_amp`, `drive_phase`, `drive_freq_ratio`, `drive_amp_perp`, and `perp_sign` are the per-flagellum overrides of `chlamy_run_case`. Each preset writes `swim_analysis.mat` beside the run.

## 3. Run with Inline Parameters

1. Open [`main_coupled_simulation.m`](main_coupled_simulation.m).
2. Edit the `cfg` block and `pili_config` entries.
3. Run:

```matlab
main_coupled_simulation
```

Use this mode when you are developing a new case interactively and want the configuration close to the entry script.

## 4. Run with JSON Configuration

1. Copy [`configs/config_example.json`](configs/config_example.json) to a local folder, for example `configs/my_case.json`.
2. Set `config_file` in [`main_coupled_simulation.m`](main_coupled_simulation.m):

```matlab
config_file = 'configs/my_case.json';
main_coupled_simulation
```

The repository tracks example JSON files under `configs/`. Copy them to paths such as `configs/my_case.json` for local edits. Files under `configs/` are not ignored by git, so a private case file shows up as untracked; run products go to `output/`, which is ignored.

### JSON Sections by Geometry Mode

**Parametric mode** (`geometry_mode: "parametric"`): use `body`, `pili`, `material`, `fluid`, `simulation`, `brownian`, `output`, `flowfield`, and `background_flow`. Do **not** set a non-empty `custom_geometry` block (since v7.4.4). Authoritative example: [`configs/config_example.json`](configs/config_example.json).

**Custom mode** (`geometry_mode: "custom"`): provide `custom_geometry` with imported body and pili point clouds. Do not set parametric body fields (`major_axis`, `N_body`, …) or a non-empty `Pili_Matrix`. Authoritative example: [`configs/config_custom_example.json`](configs/config_custom_example.json).

Since v7.4.4, `validate_simulation_params.m` rejects **both directions** of custom/parametric mixing with explicit errors when mixed fields are present on `cfg`.

**JSON note:** `load_sim_config.m` does not map a top-level `custom_geometry` block into `cfg` when `geometry_mode` is `"parametric"`. Stray `custom_geometry` keys in a parametric JSON file are therefore **silently ignored at load time** and never reach validation—delete the block from the JSON or set `"geometry_mode": "custom"`.

A typical top-level layout includes:

- `case_name`
- `geometry_mode`
- `body` and `pili` (parametric) **or** `custom_geometry` (custom)
- `material`, `fluid`, `simulation`, `brownian`, `output`, `flowfield`, `background_flow`

The `material` section also accepts the optional key `hydro_radius_factor` (default `sqrt(2)`, also `cfg.hydro_radius_factor`): the pilus hydrodynamic regularisation length is `eps_pili = hydro_radius_factor * r_pili`, while the bending stiffness always uses the physical `r_pili`. `sqrt(2)` reproduces the no-slip-cylinder twist drag exactly and the translational mobility within about 3 %; `sqrt(e)` makes the translation exact; `1` reproduces versions up to 8.5.1. Values outside [1, 3] raise a validation warning.

### Configuration Enumerations (from code)

| Field | Allowed values | Source |
| --- | --- | --- |
| `intrinsic_curvature` | bare string (`frozen`, `zero`, `constant`, `helix`) **or** struct with `.type` | `validate_simulation_params.m` |
| `intrinsic_curvature` (`constant`) | struct with 3-element `kappa` | `validate_simulation_params.m` |
| `intrinsic_curvature` (`helix`) | struct with finite scalar `kappa_d2`, `kappa_d3` | `validate_simulation_params.m` |
| `active_moment` | bare string (`off`, `traveling_wave`, `custom`) **or** struct with `.type` | `validate_simulation_params.m` |
| `motor_torque` | bare string (`off`, `constant`) **or** struct with `.type` | `validate_simulation_params.m` |
| `background_flow.type` | `off`, `uniform`, `shear_z` | `coupled_background_flow.m` |
| `root_coupling` | `clamped`, `free_torsion` | `build_coupled_lin_sys.m` |

### Units of Pilus Curvature and Moment Inputs

The elastic-rod equations of each pilus are nondimensionalised with the pilus
length $L_{\mathrm{pili}}$ (not the global body scale $L_0$): arc length runs
over $s \in [0, 1]$ with node spacing $1/N_{\mathrm{seg}}$, so every curvature
is a dimensionless quantity in units of $1/L_{\mathrm{pili}}$ and every
distributed moment is in units of $EI/L_{\mathrm{pili}}$.

| Input | Meaning | Unit / conversion |
| --- | --- | --- |
| `intrinsic_curvature.kappa = [k1, k2, k3]` (`constant`) | reference $\kappa_1 = \mathbf{d}_3\cdot\mathbf{d}_2'$ (bending about $\mathbf{d}_1$), $\kappa_2 = \mathbf{d}_1\cdot\mathbf{d}_3'$ (bending about $\mathbf{d}_2$), $\kappa_3 = \mathbf{d}_2\cdot\mathbf{d}_1'$ (twist) | dimensionless; physical value $= k_i / L_{\mathrm{pili}}$ [1/m] |
| `intrinsic_curvature.kappa_d2`, `kappa_d3` (`helix`) | uniform reference bending about $\mathbf{d}_2$ and reference twist | dimensionless; for a helix of radius $R$ and reduced pitch $p$ (pitch $=2\pi p$): `kappa_d2` $= L_{\mathrm{pili}}\,R/(R^2+p^2)$, `kappa_d3` $= L_{\mathrm{pili}}\,p/(R^2+p^2)$ |
| `active_moment.amplitude_Nm` | internal active moment amplitude | physical [N m]; `amplitude_dimless = amplitude_Nm * L_pili / EI` |
| `motor_torque.torque_pN_nm` | motor stall torque | physical [pN nm]; `torque_dimless = torque / (EI/L_pili)` |
| `motor_torque.max_omega_Hz` | zero-torque (run-away) speed | physical [Hz]; `max_omega_dimless = 2*pi*max_omega_Hz*T_0` |

Because curvature inputs are scaled by the pilus length, the same `kappa`
value on pili of different `length_um` describes different physical
curvatures.  With `intrinsic_curvature = "frozen"` (default) the reference
curvature is taken from the initial geometry and no unit conversion is needed.

### Pre-flight Parameter Checks

[`validate_simulation_params.m`](src/config/validate_simulation_params.m) runs before initialisation and aggregates all errors into a single report. The checks that most often stop a run:

| Check | Condition | Severity |
| --- | --- | --- |
| `Pili_Matrix` shape | not exactly eight rows (one column per pilus); a ninth, per-pilus stiffness row is rejected because all pili share one `E` and `r_pili` | error |
| `pili_ext` count consistency | `numel(pili_ext)` != number of pili (parametric: `size(Pili_Matrix,2)`; custom: `numel(custom_geometry.pili)`) | error |
| `delta` range | `delta <= 0` | error |
| `delta` range | `delta >= pi/4` (premature singularity event) | warning |
| `N_body` upper bound | `N_body > 2000` (`O(N^2)` memory) | warning |
| Background-flow orthogonality | `shear_z` requires finite nonzero `direction` and `shear_axis` with `abs(dot(direction_hat,shear_axis_hat)) <= 1e-12` | error |
| `alpha_0` range | outside `[-pi/2, pi/2]`; enforced in parametric mode and custom mode with empty `initial_orientation` | error |
| `AbsTol` / `RelTol` range | outside `(0, 1]`; recommended `1e-6 ~ 1e-3` | warning |

Brownian JSON fields mapped by [`load_sim_config.m`](src/config/load_sim_config.m): `brownian.enable`, `brownian.temperature` -> `Temper`, `brownian.rng_seed`, `brownian.dt_B_ratio`, `brownian.pili_thermal` -> `enable_pili_brown`, `brownian.pili_mode_orders` -> `pili_brown_orders`.

The authoritative schemas are the tracked example files and the loader implementation in [`src/config/load_sim_config.m`](src/config/load_sim_config.m).

## 5. Programmatic and Batch Runs

For scripts or sweeps, call the simulation engine directly:

```matlab
cfg = load_sim_config('configs/config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

Batch sweep **driver scripts** live in [`batch_sweep/`](batch_sweep/). The folder includes parameter sweeps for body resolution, body regularization, pilus length, pilus stiffness, segment density, initial angle, and Brownian repeats.

When those sweep drivers run, CSV/MAT summary products are written under **`batch sweep/`** (with a space) at the repository root. This output folder name is intentional in the sweep scripts and is separate from the source folder `batch_sweep/`.

### 5.1 Brownian repeat batch (`batch_repeat_brownian.m`)

From the repository root, open [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m). It runs independent repeats of the **same** Brownian-on configuration into `output/<case_name_base>/001`, `002`, … (zero-padded).

- **Parallel (default):** `use_parallel = true` runs a `parfor` loop when **Parallel Computing Toolbox** is licensed and a `parpool` can be created. Worker count is limited by `max_resource_frac` (default `0.90`) applied to probed physical cores and to available RAM via `mem_per_worker_gb` (default `2.5`). If PCT is missing or pool creation fails, the driver **warns and falls back** to the legacy serial `for` loop.
- **Serial:** `use_parallel = false` keeps the pre-v7.4 single-threaded loop (use this for strict reproducibility checks against older behaviour).
- **RNG / seeds:** In parallel mode, `brown_rng_seed_mode = 'fixed'` sets `cfg.brown_rng_seed = rng_seed_base + k` per repeat. For `'timestamp'`, the script assigns a clock-derived base plus `k` so workers cannot share the same wall-clock seed. In **serial** `'timestamp'` mode, `cfg.brown_rng_seed` is still omitted each run (same policy as a single interactive `run_simulation` call).
- **Worker console silence (parallel only):** the parfor branch sets `cfg_k.silent = true`, so each worker's `run_simulation` writes its full log to `output/<case>/<NNN>/simulation.log` only and does **not** echo to the shared stdout. This prevents ~15000 interleaved `[ODE] Segment N completed` lines from burying the useful `[k/100] OK` progress lines. Custom-geometry runs are fully silent on workers as of v7.9.5. The serial branch is unchanged (verbose) so interactive debugging remains informative.

A timestamped `batch_repeat_summary_<batch_id>.txt` under the case base folder records batch metadata, per-run status, timing, and parallel diagnostics. **Header (v7.9.5+):** `batch_id`, human-readable `started`, `config_source_type` / `config_json` or `config_source_note`, path to `batch_config_snapshot_<batch_id>.json` (full inline `cfg` copy for reproducibility), `project_root`, `script_path`, `brown_rng_policy`, `parallel_requested`, `execution_mode`, `summary_record_order`, optional `[POOL] Resource probe` fields (including on serial fallback), and `seeds_precomputed_unused` when parallel seeds were computed but serial ran. **Per-run records** are written **incrementally** (one `fopen(...,'a')` / `fclose` per run): OK lines include `out:` (relative path); FAILED lines include `Identifier` / `Message` / `Location`, and when available `outdir` / `log_hint` pointing at partial output. **Footer:** `finished`, `elapsed_total_s`, `success_count`, parallel speedup fields when applicable, `failed_count` / `failed_indices` / `failed_runs`, and `pending_count` / `pending_indices` for runs still `PENDING` after an abnormal exit. A mid-batch crash leaves the header plus all completed-run records on disk; only the footer may be missing. In parallel mode the per-run lines appear in worker-completion order (the leading `NNN` index preserves `k` for sorting); serial mode writes in `k` order.

### 5.2 JSON folder batch (`batch_run_from_json.m`)

From the repository root, run [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m). Select a folder of `*.json` files (or set `BATCH_JSON_DIR` at the top of the script). Each queued file is loaded with `load_sim_config`, forced headless (`auto_confirm = true`, `preview_initial = false`), and executed via `run_simulation`. Brownian RNG, animation, and flow-field flags follow each JSON.

- **Pre-scan:** files are sorted alphabetically. Duplicate **non-empty** `case_name` values are marked `SKIPPED` (first file in sort order runs; later duplicates are not queued). Empty `case_name` entries are never deduplicated (`run_simulation` assigns unique timestamped output folders).
- **Parallel (default):** `use_parallel = true` with `mem_per_worker_gb = 3.0`, `max_resource_frac = 0.75`, and `worker_blas_threads = 1`. Uses `parfor` when **Parallel Computing Toolbox** is available; otherwise warns and falls back to serial.
- **Real-time log:** `batch_summary_<timestamp>.log` in the JSON folder receives progress lines via `parallel.pool.DataQueue` (parallel) or per-iteration writes (serial), with `fflush` after each event.
- **Summary:** console and log report Total JSON files, Queued to run, Succeeded / Failed / Skipped, steady `V_sink` per file, and (when parallel) a speedup estimate.

```matlab
run('batch_sweep/batch_run_from_json.m');
```

### 5.3 Linux HPC cluster batch (`batch_repeat_brownian_HPC_Linux.m`, since 7.10.0)

For headless, single-node-exclusive Linux HPC cluster execution, use [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m). Since v10.0.0, all Linux/headless engine features (platform-aware AVI/MP4 video export, `.avi` detection, `cfg.silent`, and `DataQueue` progress heartbeats) live directly in `run_simulation.m`, and [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) is a thin backward-compatible forwarding wrapper to `run_simulation.m`. Usage is otherwise identical to §5.1 (same `case_name_base` / `n_repeats` / `brown_rng_seed_mode` / inline-or-JSON config settings).

Differences from `batch_repeat_brownian.m`:

- **Platform-aware video export:** `VideoWriter` uses `'MPEG-4'` (`.mp4`) on Windows/Mac as before, but `'Motion JPEG AVI'` (`.avi`) on Linux — MathWorks documents `'MPEG-4'` as Windows/Mac-only; it throws on any Linux system regardless of MATLAB release or display availability.
- **Linux-native resource probing:** physical-core count falls back to `lscpu -p=Core,Socket` and available memory falls back to `/proc/meminfo`'s `MemAvailable` field when the Windows-only `wmic`/`memory()` probes are unavailable, replacing a previous silent "assume 4 GB available" fallback that badly under-sized the parallel pool on real Linux nodes.
- **Independent CPU/memory headroom:** `cpu_resource_frac` (default `0.90`) and `mem_resource_frac` (default `0.65`, kept more conservative) replace the single `max_resource_frac`.
- **Heavy-output sampling:** `generate_animation` / `compute_flowfield` are forced off for every repeat except `heavy_output_sample_idx` (first, middle, and last by default), regardless of the inline config or JSON, to bound per-worker peak memory. Set `heavy_output_sample_idx = 1:n_repeats` to restore full output for every repeat once memory headroom is confirmed safe for your configuration.
- **Node-local `parpool` job storage:** on Linux, the local pool's `JobStorageLocation` is relocated to `tempdir()` (honouring `$TMPDIR`, typically node-local scratch set by the job scheduler) in a job-unique subfolder named from `SLURM_JOB_ID` / `PBS_JOBID`, avoiding a documented pool-startup race when the MATLAB default (under the user's home directory) sits on distributed/network cluster storage. Falls back to the original behaviour if relocation itself fails. Windows is unaffected.
- **Live sub-case progress heartbeat:** while `cfg.silent = true` suppresses all normal console echo for a running repeat, a throttled (≥ 30 s per repeat) `[HB] k=NNN running  PP.P%  elapsed=...` line is still printed, so a long-running repeat's progress is visible without waiting for it to finish. Console-only; not written to the summary file.

**Before trusting the auto-sized pool on a real cluster allocation**, `mem_per_worker_gb` (default `2.5` GB) is **not** empirically calibrated for any specific configuration. Run one standalone repeat and measure its real peak RSS (e.g. `/usr/bin/time -v matlab -batch ...` or `top`/`ps` on Linux), then set `mem_per_worker_gb` to that measured value plus a 30–50% margin. Under-budgeting this is the most likely root cause of workers being killed by the Linux OOM-killer mid-batch — since `parfor` has no built-in timeout or dead-worker recovery, this can leave the batch hanging indefinitely rather than erroring out. Heavy-output sampling (above) reduces the likelihood of hitting this; it does not eliminate it.

## 6. Custom Geometry Mode

Use custom mode when both the body and pili are imported from point-cloud files.

1. Copy [`configs/config_custom_example.json`](configs/config_custom_example.json).
2. Set `"geometry_mode": "custom"`.
3. Provide body and pili files under `"custom_geometry"`.
4. Use the example files in [`geometry/`](geometry/) as format references.

### Custom Geometry Rules

- Custom body and custom pili must be used together.
- Body point clouds support `.csv`, `.txt`, and `.mat` files with three spatial coordinates.
- Imported pili centerlines may have nonuniform source spacing.  Consecutive
  duplicates are removed, then a PCHIP curve is converted to one canonical
  equal-chord polygon while preserving both endpoints.
- Per-pilus `bending_stiffness` is not supported.  Every pilus uses the global
  `material.E` and `r_pili`, so a non-empty legacy field is rejected.
- For arbitrary body orientation, prefer `initial_orientation` over `alpha_0`.
- Custom-mode `deflection_z` is a twist about the aligned tangent, while custom-mode `deflection_x` is a tilt. In parametric mode these geometric roles are swapped.

## 7. Motor and Free-Torsion Runs

Use [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json) as the starting template.

Key fields under each pilus entry:

- `root_coupling`: `"clamped"` or `"free_torsion"`.
- `motor_torque.type`: `"off"` or `"constant"`.
- `motor_torque.torque_pN_nm`: motor stall torque in pN nm (signed, see below).
- `motor_torque.max_omega_Hz`: zero-torque speed of the linear torque-speed law in Hz.

For `free_torsion`, the root twist angle $\psi_1$ evolves dynamically. Its equation is the axial ($d_3$) moment balance of the **whole filament about the root**: the hydrodynamic torque on the filament is balanced by the torque the body transmits through the root, which is the motor torque (or zero for a frictionless bearing when `motor_torque.type = "off"`). The internal Kirchhoff twist moment does not appear in this row — it is carried by the interior nodes (v8.4.0; before that the interior twist moment was also added to the root row, which made the motor torque balance the drag of the first segment only and inflated the steady rotor speed by a factor of the segment count $N$). No separate reaction torque is added to `R_global(4:6)` because those rows represent the angular-momentum balance of the full coupled body-pilus system.

In steady state the rotor speed of a straight filament is set by the whole-filament torsional drag, $\omega \approx \tau_{\mathrm{eff}} / (\xi_t L)$ with $\xi_t$ the axial rotational drag per unit length of the regularized-Stokeslet-segment model (for the default $\epsilon = \sqrt{2}\,r_{\mathrm{pili}}$, set by `material.hydro_radius_factor`, this evaluates to $\xi_t = 2\pi\mu\epsilon^2 = 4\pi\mu r^2$, the no-slip-cylinder value; up to 8.5.1 the factor was 1 and $\xi_t$ was half of it.  The pre-flight stall diagnostic in `validate_simulation_params.m` uses the same coefficient).

**Torque-speed law (v8.2.0).** The effective motor torque follows the
measured linear torque-speed relation of the bacterial flagellar motor:

$$\tau_{\mathrm{eff}} = \tau_{\mathrm{motor}}\left(1 - \frac{\omega_{\mathrm{rel}}}{\omega_{\max}}\right),\qquad
\omega_{\mathrm{rel}} = \operatorname{sign}(\tau_{\mathrm{motor}})\,\big[(\dot\psi_1 + \cos\theta_1\,\dot\phi_1) - \boldsymbol{\Omega}_{\mathrm{body}}\cdot\mathbf{d}_3^{(1)}\big].$$

- $\omega_{\mathrm{rel}}$ is the rotor (segment 1) spin **relative to the cell body**, measured positive in the motor's driving direction. It uses the full axial angular velocity of segment 1, not $\dot\psi_1$ alone.
- The law is asymmetric: $\omega_{\mathrm{rel}} = 0$ gives the stall torque, $\omega_{\mathrm{rel}} = \omega_{\max}$ gives zero torque, forcing the rotor backwards ($\omega_{\mathrm{rel}} < 0$) raises the torque above stall, and spinning it faster than $\omega_{\max}$ turns the motor into a brake ($\tau_{\mathrm{eff}} < 0$). The former $|\omega|$ form treated both rotation senses identically.
- Because the law is linear, the self-consistent torque is obtained exactly from one small linear solve per RHS call (any number of motors).

**Sign convention.** A positive `torque_pN_nm` drives the rotor in the right-handed sense about the root tangent $\mathbf{d}_3$ pointing away from the body, i.e. counter-clockwise when viewed from the filament tip looking toward the cell (the *E. coli* run sense for a left-handed filament). A negative value drives it clockwise. Older configurations (before v8.2.0) had the opposite, undocumented sense.

## 8. Brownian Dynamics

Set `brownian.enable` in JSON or `cfg.enable_brown` in MATLAB.

By default Brownian dynamics is body-only: the 6×6 effective mobility `M_eff`
is obtained from the complete coupled body-pili deterministic system, but
thermal forcing is applied only to the body's three translations and three
rotations.  Pili internal Euler-angle modes receive no thermal noise of their
own.  This is not a full coupled thermal-equilibrium model and should not be
used to study pilus-mode equipartition.  The option `brownian.pili_thermal`
(v9.0.0, below) adds thermal forcing of the slow pilus bending modes.

**Force realisation (v8.5.0).**  Each Brownian interval `dt_B` starts by
probing `M_eff`, drawing the Euler–Maruyama increment
`dU = L_eff * eta * sqrt(dt_B)`, and converting it into the constant
generalised thermal force `f_th = M_eff \ (dU / dt_B)` that is added to the
body rows of the ODE right-hand side for the whole interval.  The body is
therefore moved by the coupled solver, and the pili follow it continuously
through their hydrodynamic and elastic coupling exactly as they follow
gravity or a motor torque.  Two practical consequences:

- No spurious root kink.  Earlier versions displaced the body instantaneously
  and re-projected only the first pilus segment onto the moved body; every
  kick injected a curvature jump between segments 1 and 2 that excited stiff
  high-wavenumber elastic modes and forced `ode15s` into very small steps
  after each node.  With force realisation the root constraint is never
  violated, so Brownian runs no longer need a shorter `dt_B` or finer
  segmentation to stay smooth.
- Implicit Ito drift.  The mobility gradient sampled along the trajectory
  reproduces the thermal drift `kBT*div(M_eff)` to first order in `dt_B`
  (see `docs/theory.md` §10.2), so the explicit finite-difference drift
  (12 extra coupled solves per node in v7.5.0–v8.4.0) is no longer computed.
  The `[BROWN]` log reports the sampled `|F_th|`, `|tau_th|` and the
  frozen-mobility target increment at the start of the first five intervals
  and every 100th interval thereafter.

SPD repair and Cholesky factorisation are performed in dimensionless
generalized coordinates before conversion to SI units.  The final Brownian
interval uses the actual remaining time when `T_obs` is not divisible by
`dt_B`.  There is no user-facing midpoint drift option.  Runs with a fixed
`brownian.rng_seed` remain reproducible, but trajectories differ from
v8.4.0 and earlier because the stochastic increment is realised differently.

Important options:

- `brownian.temperature`: positive fluid temperature in Kelvin; required when Brownian dynamics are enabled.
- `brownian.rng_seed`: fixed positive seed for reproducible runs.
- `brownian.dt_B_ratio`: Brownian time step relative to output step.
- `brownian.pili_thermal` (`cfg.enable_pili_brown`, default `false`): thermal fluctuations of the pili, see below.
- `brownian.pili_mode_orders` (`cfg.pili_brown_orders`, default `2`): number of bending orders per pilus that receive thermal noise.
- `brownian.pili_tau_cut_s` (`cfg.pili_brown_tau_cut_s`, optional, default unset): time-scale cutoff in seconds — keeps every bending mode with relaxation time `>= tau_cut_s` (recommended for unequal-length pili; enables pair protection by default).
- `brownian.pili_protect_pairs` (`cfg.pili_brown_protect_pairs`, optional): expand the mode cut across near-degenerate clusters (`gap_ratio < 1.15`) so a polarisation pair is never split.
- `brownian.pili_rng_seeds` (`cfg.brown_pili_rng_seeds`, optional): positive integer seeds (one per pilus, distinct from `brownian.rng_seed`) for independent per-pilus thermal streams; automatically assigned non-overlapping blocks across repeats in `batch_repeat_brownian*.m`.

**Thermal fluctuations of the pili (v9.0.0, upgraded in v10.0.0).**  With
`brownian.pili_thermal = true` the slow bending modes of the pili (by default
the first and second bending orders of every pilus, two polarisations each)
receive thermal forces of their own.  The modes are computed once from the
coupled operator at the initial state.  In every interval the thermal forces
on the body and on these modes are drawn together from the joint mobility and
realised through the constrained probe's centreline force density together
with its balancing scalar axial torque density (`brown_pili_tau`), so body–
filament correlations are kept, non-thermal modes are not kicked, and each
thermal mode relaxes to an elastic energy of `kT/2`.  Twist and higher
bending modes stay deterministic.  The option needs `brownian.enable = true`
and at least one pilus.

- **Choose `dt_B` inside the two-sided window printed in the log.**  The
  `[BROWN]` lines at start-up list the relaxation times of the thermal modes
  and the slowest non-thermal bending mode (`tau_det`).  Keep `dt_B` at or
  below a quarter of the fastest thermal relaxation time (0.002 s for 3 µm
  filaments with `configs/config_example.json` material; `run_simulation`
  warns for `x_fast > 0.5` and rejects `x_fast > 2.0` unless
  `cfg.allow_large_pili_dt_B = true`), and do not shrink `dt_B` well below
  `tau_det` without raising `pili_mode_orders` or lowering `pili_tau_cut_s`
  (`run_simulation` warns when `dt_B / tau_det < 0.5`).  Tiny closing
  remainders (`< 5 %` of `dt_B`) are absorbed into the final interval.
- **Outputs.**  `results.pili_modes` holds the mode amplitudes `a`
  (`energy_kT = a.^2 / (2 kT)`, equilibrium mean 1/2), total bending and
  twist energies (`elastic_energy_kT`), relaxation times, orders,
  `dom_pilus`, `pilus_order`, `selection_mode`, `pilus_share`, and run-time
  maximum diagnostics (`max_pattern_resid`, `max_hold_pattern_resid`,
  `max_twist_row_resid`, `max_hold_force`).
  `results.flow_state.brown_mode_force`, `brown_pili_force`, and
  `brown_pili_tau` store the modal forces, centreline force density, and
  balancing axial torque density at every output for exact flow-field redraws.
- **Limits.**  The scheme is built and verified for naturally straight,
  clamped, passive filaments without gravity or background flow.
  `root_coupling = free_torsion` (a rotating flagellum, with or without
  motor) is rejected with an error; use the body-only scheme for such
  cells.  Intrinsic-curvature overrides and active moments run with a
  warning: the thermal modes are those of the passive initial state, a
  force on the centreline leaves the twist untouched only for a naturally
  straight filament, and equipartition has not been verified there; check
  `elastic_energy_kT(2, :)` (twist) before relying on such a run.
  Modes faster than the thermal ones are treated as already relaxed when
  the noise is drawn, so the body diffusion below about 2 ms (the
  relaxation time of the third bending order of a 3 µm filament) is not
  resolved.  From the lag `dt_B` on it follows the linear theory, and at
  long times it equals the value for a cell of frozen shape
  (`docs/theory.md` §10.4 gives the numbers).  Raising `pili_mode_orders`
  to 3 needs `dt_B` of about 0.0005 s.

**Lagged Jacobian (v8.5.1).**  Brownian runs restart `ode15s` at every
`dt_B` node, and every restart used to pay a dense finite-difference
Jacobian (`n_ode` extra right-hand-side evaluations, 70–85 % of the run time
of a four-pilus case).  Since v8.5.1 the last exact Jacobian is re-used
([`src/core/LaggedJacobian.m`](src/core/LaggedJacobian.m)): it enters only the
Newton iteration matrix of the implicit step, so the solution stays within the
integration tolerance of the previous result while a Brownian node costs about
2–3.7× less (1 pilus ≈ 1.4×).  A fresh exact Jacobian is computed when
`ode15s` asks for one again within the same call (Newton failure), after an
Euler-basis rotation, and when the cache is 50 nodes old.  If a lagged
segment errors or stops prematurely it is re-run with the default Jacobian.
The `[JAC]` lines of `simulation.log` and `results.jacobian_reuse` report the
counters (`n_fresh`, `n_lagged`, `n_on_demand`, `n_invalidated`,
`n_fallback`).  Runtime-only options (not in JSON): `cfg.jacobian_reuse`
(default `true`; `false` reproduces v8.5.0 bit-for-bit) and
`cfg.jacobian_max_age` (default 50 nodes).  Deterministic runs are unaffected.

When Brownian dynamics are enabled, the realized seed is saved in `config.json`, `brownian_rng_seed.txt`, `simulation.log`, and `simulation_results.mat`.  It seeds a run-local stream, so cache creation and unrelated calls to MATLAB's global random stream cannot change the trajectory.

**Mobility conditioning warning (since v7.4.1):** Each Brownian step evaluates an `rcond`-based diagnostic inside [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) on the **dimensionless** effective mobility after symmetrisation (`M_eff_dimless`). Routine runs should therefore **not** spam `compute_effective_body_brownian:IllConditioned` from dimensional TT/TR/RR scaling alone; if MATLAB still emits this warning with a very small **`M_eff_dimless` rcond**, investigate degenerate kinematics or geometry (e.g. Euler roots near \(\theta = \pi\) for clamped pili rows, unusually extreme body aspect ratios) with `simulation.log` alongside the kinematic trajectory.

## 9. Main Outputs

Outputs are written under `output/<case_name>/` or a timestamped run directory.

| File | Purpose |
| --- | --- |
| `simulation_results.mat` | Full `results` and `params` structures. |
| `timeseries.csv` | Scalar time series for external analysis. |
| `config.json` | Saved configuration snapshot for reproduction. |
| `summary_plot.png` / `.pdf` / optional `.html` | Main summary figure. |
| `flowfield_panels.png` / `.pdf` / optional `.html` | Final XZ, YZ and XY flow-field panels (2×3) from the run. |
| `flowfield_panels_rerender_final.png` | Optional final-state panel from `render_flowfield_from_results`. |
| `flowfield_panels_tNNNNNN.png` | Optional selected-time panels from `render_flowfield_from_results`. |
| `flowfield_animation.mp4` | Optional flow animation (`.avi` on Linux). |
| `initial_pose.png` / `.pdf` / optional `.html` | Initial-pose preview when enabled. |
| `bacteria_3d_motion.mp4` | Optional 3D animation (`bacteria_3d_motion.avi` on Linux). |
| `bacteria_3d_motion_rerender.mp4` | Optional; written by `render_animation_from_results` (see below). |
| `simulation.log` | Complete run log. |
| `brownian_rng_seed.txt` | Brownian seed and reproduce hint when Brownian dynamics are enabled. |

The animation is rendered by `src/visualization/render_motion_animation.m`: all graphics objects are created once and only their data are updated per frame.  It has four panels on a fixed 16:9 canvas (2560 x 1440 logical px; the encoded size follows the display scaling, e.g. 3840 x 2160 at 150 %): **a** 3-D overview (orthographic, fixed limits, gravity arrow), **b** XZ overview (fixed limits, 2-μm scale bar), **c** co-moving 3-D window and **d** co-moving XZ window.  Windows c and d follow the centre of mass at constant magnification; they are not cubes or squares: along each lab axis the box extent is the reach of the body and all pili over all frames (relative to the centre of mass) plus an 8 % margin (at least 0.1 μm), so the bacterium fills the window, and the scale is equal on X, Y and Z (nothing is stretched).  Each window has ticks, tick labels and a grid in lab-frame coordinates (so the grid slides past the bacterium); on panels c and d the spacing is never finer than 0.5 μm.  The four panels are arranged automatically: their aspect ratios and relative sizes are fixed by the data, and the arrangement (nested side-by-side / stacked splits, reading orders a-b-c-d, a-c-b-d, a-c-d-b) that allows the largest common scale factor is chosen and written to `simulation.log` (for example `Layout H(a,H(b,V(c,d)))`), so the canvas stays compact for elongated, flat or compact bacteria alike.  The centre-of-mass trajectory is drawn in a and b only.  The body is uniform light green, with one short dark-green patch on one side of the waist (about 18 % of the body length and ±32° of azimuth). The patch rotates with the body, so a roll about the long axis can be seen in every panel.  The canvas carries a legend, a global time stamp *t* and the playback speed relative to real time (simulated seconds per video second = frame interval x `animation_fps`; labelled "slow motion" or "accelerated"), and panel titles a–d flush on the upper-left corner of each image.  `render_motion_animation` also accepts `opts.scale_bar_um` (default 2) and returns `window_half_um` (box half-extents `[hx hy hz]`), `window_offset_um` and `playback_speed` in its `info` structure.  `simulation.log` ends the animation phase with a per-stage timing table (axis-limit pre-scan, figure setup, data update, `drawnow`, `getframe`, frame resize, `writeVideo`), which shows where the render time goes.

To regenerate the animation from a completed case without re-running the ODE (for example with another frame stride or frame rate), run `addpath('post_processing')` and then `render_animation_from_results()` (file dialog) or `render_animation_from_results(mat_file, fps, skip)`.  Without `fps` / `skip`, the values are taken from the case folder's `config.json` (`output.animation_fps` / `output.animation_skip`), else 30 and 1.  The result is saved next to the `.mat` file as `bacteria_3d_motion_rerender.mp4` (`.avi` on Linux); the original video is never overwritten.

To recompute flow fields from a run saved by v8.7.0 or later, use `render_flowfield_from_results()` or `render_flowfield_from_results(mat_file, mode, requested_times_s, fps, skip)`.  `mode` is `final` (default), `times`, `animation`, or `all`.  Times are in seconds and match the nearest saved output; duplicate matches are drawn once.  `fps` and `skip` follow the same `config.json` fallback as the motion animation.  The colour scale is the 99th percentile of that final frame, of the requested times, or of every saved output for an animation, taken over the XZ, YZ and XY planes together.  From v8.7.3 the figure has six panels on a 2500 x 1400 canvas: 3-D views of XZ, YZ and XY on the top row and the matching 2-D views below; the colour bar sits beside panel f.  Pass `opts.drop_first_frame = true` in the optional sixth argument to omit saved output 1 from `times`, `animation` and `all` (`final` is unchanged).  Since v10.2.0, flow-field computation automatically leverages Block Schur force extraction, batched 3-plane FP64 GPU velocity evaluation (`opts.use_gpu = 'auto'`; set `false` to force CPU), thread-pool parallel batching for forces and 2-D streamlines (`opts.use_parallel = 'auto'`; set `false` for serial), and a slim ~1.6 MiB/frame grid cache (~1,300 frames in 2 GiB) with post-encode frame release; pass `opts.scale_from_video_frames = true` (or set it in `config.json`'s `flowfield` block) to compute colour-scale percentiles from rendered video frames only when `skip > 1`.  An animation prints one in-place Compute bar while the fields are evaluated and one Video bar while frames are written.  Output is written beside the `.mat` file and does not overwrite `flowfield_panels.png`.  The animation file is `flowfield_animation.mp4` on Windows and macOS, and `flowfield_animation.avi` on Linux.  On a headless machine, pass the MAT path; the file dialog needs a display.

`results.completed=false` with `results.termination_reason='timeout'`
identifies a segment-timeout return.  The available trajectory is saved and
returned normally, but animation/flow-field generation and steady-state
batch statistics are skipped.

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

The post-processor reads `timeseries.csv` and `config.json` (and optionally `results.q` from `simulation_results.mat` for orientation). For each selected case folder (one trajectory) it computes:

| Metric | Method |
| --- | --- |
| Sinking velocity | Mean / last-5% steady / std of \(-\texttt{Vsink\_um\_s}\) |
| Max horizontal displacement, 3-D tortuosity | Path geometry from \(X,Y,Z\) |
| Translational \(D_\mathrm{trans}\), \(V_\mathrm{msd}\) | 3-D TAMSD; fit \(\mathrm{MSD}=6D\tau+V^2\tau^2\) on \(0<\tau\le T/2\) (lags span full \(T\)) |
| Rotational \(D_\mathrm{rot}\) | XZ body-axis MSAD; fit \(\mathrm{MSAD}=2D_\mathrm{rot}\tau\) on \(0<\tau\le T/2\) |

Results go to one Excel workbook plus a multi-panel summary figure (PNG, PDF, optional HTML on R2026a+). Partial timeout cases (`results.completed=false`) are excluded from aggregate diffusion metrics.

## 11. Common Checks

| Symptom | Check |
| --- | --- |
| MATLAB cannot find a function | Confirm the repository root is on the MATLAB path. |
| Capsule generation fails | Confirm the S2 Sampling Toolbox is installed and visible on the MATLAB path. |
| JSON case runs without pili | Check that `pili` is a valid JSON array and that each pilus entry has the required fields. |
| Custom geometry is rejected | Check point-cloud orientation, duplicate/degenerate intervals, and canonical equal-chord construction convergence. |
| `geometry_mode` mixing error | Use parametric **or** custom fields only; do not combine `custom_geometry` with `body`/`pili` matrix fields (v7.4.4+). |
| `helix` intrinsic curvature error | Provide a struct: `{"type":"helix","kappa_d2":…,"kappa_d3":…}` with finite scalars; bare `"helix"` is invalid (v7.4.4+). |
| Output folder already exists | Use `auto_confirm = true` for unattended rename behavior, or choose interactively. |
| Brownian repeat results differ | Confirm whether the seed mode is fixed or timestamp-based. |
| Parallel Brownian batch runs serially | Install **Parallel Computing Toolbox**, fix `parpool` cluster settings, or set `use_parallel = false` in `batch_repeat_brownian.m`. |
| Parallel JSON folder batch runs serially | Same as above for `batch_run_from_json.m`, or set `use_parallel = false`. |
| JSON batch skipped as duplicate `case_name` | Only the alphabetically first JSON with that non-empty `case_name` runs; assign unique `case_name` values or leave empty for auto timestamp folders. |
| `IllConditioned` Brownian mobility warning (`M_eff_dimless rcond`) | v7.4.1+ diagnoses the dimensionless mobility; residual warnings mean genuinely ill-conditioned kinematics/hydrodynamics — cross-check Euler roots (\(\theta\approx\pi\)), body geometry, and `simulation.log` (see §8). |
| `NearSingularMatrix` warning (`rcond < 1e-14`) | v8.0.1+: recorded only; does **not** abort by itself. Abort requires joint large backward residual (`> 1e-6`) → `SingularMatrix`. Check `results.diagnostics` / `[DIAG]` log lines. |
| Warning `numerical_diagnostics:DissipationEigensolveFallback` | Harmless. The dissipation check uses the dense `eig`; the iterative `eigs` is tried only for a sparse operator with more than 600 rows, and if it does not converge the check falls back to `eig` with this warning. The run is not affected. |
| `run_simulation:SilentPreviewDeadlock` error (v7.9.6+) | `cfg.silent = true` with `preview_initial = true` requires `auto_confirm = true` (a headless/parallel worker cannot answer the interactive `questdlg`). Set `auto_confirm = true` or `preview_initial = false`. |
| `VideoWriter` `"The specified profile is not valid."` on Linux | Use [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) (§5.3) instead of `run_simulation.m` — `'MPEG-4'` is Windows/Mac-only; the HPC_Linux fork switches to `'Motion JPEG AVI'` on Linux automatically. |
| Worker process count drops mid-batch (`top`), batch hangs indefinitely (Linux) | Likely the Linux OOM-killer; `parfor` has no built-in timeout / dead-worker recovery. On `batch_repeat_brownian_HPC_Linux.m` (§5.3), calibrate `mem_per_worker_gb` from a real measured peak RSS (not the default) and/or reduce `heavy_output_sample_idx` coverage. Confirm with `dmesg \| grep -i "killed process"` or `journalctl -k \| grep -i "out of memory"`. |
| Slow `parpool` startup or pool-creation failure on a shared-home-directory cluster | On `batch_repeat_brownian_HPC_Linux.m` (§5.3), `ensure_parpool` already relocates `JobStorageLocation` to node-local `tempdir()` on Linux; confirm the console shows `[INFO] parpool JobStorageLocation (node-local): ...`. If it instead shows the relocation-failed warning, the fallback uses MATLAB's default location, which is a known slow-startup/race risk when the home directory is distributed/network storage. |

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
