# User Guide — Bacteria Elastohydrodynamic Simulation

English | [中文](USER_GUIDE.zh-CN.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Version:** 8.0.0
> **Updated:** 2026-07-13

This guide is for running simulations, choosing the right configuration path, finding outputs, and post-processing completed cases. For model equations and derivations, see [`docs/theory.md`](docs/theory.md). For release history, see [`CHANGELOG.md`](CHANGELOG.md).

Release note: v8.0.0 is a major mathematical-physics correctness and runtime-contract release. It introduces weighted-reciprocal translation/twist hydrodynamics for unequal segment lengths, equal-chord canonical custom-pilus geometry, objective full-director initialisation, dimensionless Brownian SPD repair, true-duration remainder kicks, isolated Brownian random streams, strict incompressible-shear validation, saved partial timeout returns, and run-scoped numerical diagnostics.

Breaking configuration changes: `Pili_Matrix` must contain exactly eight rows; per-pilus `bending_stiffness` is no longer supported; and `shear_z.direction` must be orthogonal to `shear_z.shear_axis`. A segment timeout now saves the achieved trajectory and returns normally with `results.completed=false` and `termination_reason='timeout'`; callers must not use that partial trajectory for steady-state statistics.

## 1. Before You Run

### Requirements

- MATLAB R2020b or later is recommended.
- PNG and PDF exports use `exportgraphics`; optional interactive HTML exports require MATLAB R2026a or later.
- Capsule body surface generation requires Anton Semechko's **S2 Sampling Toolbox** from MATLAB Add-On Explorer / File Exchange.
- No Python, Conda, or pip environment is required. [`requirements.txt`](requirements.txt) documents this explicitly.
- **Parallel Computing Toolbox** is optional: needed for the default parallel branches of [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m) and [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m) (`use_parallel = true`).
- **Running on a Linux HPC cluster (v7.10.0+):** use [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) / [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m) instead of the originals — see §5.3.

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
| Reproduce or share a case | Use a JSON file copied from [`configs/config_example.json`](configs/config_example.json). |
| Import custom body and pili geometry | Use [`configs/config_custom_example.json`](configs/config_custom_example.json). |
| Test a free-torsion motor-driven flagellum | Use [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json). |
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

1. Copy [`configs/config_example.json`](configs/config_example.json) to a local folder, for example `configs/my_case.json`.
2. Set `config_file` in [`main_coupled_simulation.m`](main_coupled_simulation.m):

```matlab
config_file = 'configs/my_case.json';
main_coupled_simulation
```

The repository tracks example JSON files under `configs/`. Copy them to paths such as `configs/my_case.json` for local edits. Only `output/` is listed in `.gitignore`.

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

### Configuration Enumerations (from code)

| Field | Allowed values | Source |
| --- | --- | --- |
| `intrinsic_curvature` | bare string (`frozen`, `zero`, `constant`, `helix`) **or** struct with `.type` (v7.9.0 P7) | `validate_simulation_params.m` |
| `intrinsic_curvature` (`constant`) | struct with 3-element `kappa` | `validate_simulation_params.m` |
| `intrinsic_curvature` (`helix`) | struct with finite scalar `kappa_d2`, `kappa_d3` | `validate_simulation_params.m` |
| `active_moment` | bare string (`off`, `traveling_wave`, `custom`) **or** struct with `.type` (v7.9.0 P7) | `validate_simulation_params.m` |
| `motor_torque` | bare string (`off`, `constant`) **or** struct with `.type` (v7.9.0 P7/P8) | `validate_simulation_params.m` |
| `background_flow.type` | `off`, `uniform`, `shear_z` | `coupled_background_flow.m` |
| `root_coupling` | `clamped`, `free_torsion` | `build_coupled_lin_sys.m` |

### Pre-flight Parameter Checks (v7.9.0)

[`validate_simulation_params.m`](src/config/validate_simulation_params.m) runs before initialisation and aggregates all errors into a single report. Notable checks added or tightened in v7.9.0:

| Check | Condition | Severity |
| --- | --- | --- |
| `pili_ext` count consistency (P1) | `numel(pili_ext)` != number of pili (parametric: `size(Pili_Matrix,2)`; custom: `numel(custom_geometry.pili)`) | error |
| `delta` range (P2) | `delta <= 0` | error |
| `delta` range (P2) | `delta >= pi/4` (premature singularity event) | warning |
| `N_body` upper bound (P4) | `N_body > 2000` (`O(N^2)` memory) | warning |
| Background-flow orthogonality (P5) | `shear_z` requires finite nonzero vectors and `abs(dot(direction_hat,shear_axis_hat)) <= 1e-12` | error |
| `alpha_0` range (P6/A1) | outside `[-pi/2, pi/2]`; enforced in parametric mode and custom mode with empty `initial_orientation` | error |
| `AbsTol` / `RelTol` wording (P9) | outside `(0, 1]`; recommended `1e-6 ~ 1e-3` | warning |

Brownian JSON fields mapped by [`load_sim_config.m`](src/config/load_sim_config.m): `brownian.enable`, `brownian.temperature` -> `Temper`, `brownian.rng_seed`, `brownian.dt_B_ratio`.

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

- **Parallel (default):** `use_parallel = true` runs a `parfor` loop when **Parallel Computing Toolbox** is licensed and a `parpool` can be created. Worker count is limited by `max_resource_frac` (default `0.75`) applied to probed physical cores and to available RAM via `mem_per_worker_gb`. If PCT is missing or pool creation fails, the driver **warns and falls back** to the legacy serial `for` loop.
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

### 5.3 Linux HPC cluster batch (`batch_repeat_brownian_HPC_Linux.m`, v7.10.0+)

For headless, single-node-exclusive Linux HPC cluster execution, use [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) / [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m) instead of the originals. This is a maintained fork, not auto-synced with `run_simulation.m` / `batch_repeat_brownian.m` — any future physics/algorithm fix to the originals must be ported here manually. Usage is otherwise identical to §5.1 (same `case_name_base` / `n_repeats` / `brown_rng_seed_mode` / inline-or-JSON config settings).

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
- `motor_torque.torque_pN_nm`: motor torque in pN nm.
- `motor_torque.max_omega_Hz`: optional stall-speed limit in Hz.

For `free_torsion`, the root twist angle $\psi_1$ evolves dynamically. The motor torque is applied to the pilus root $d_3$ moment-balance row only; no separate reaction torque is added to `R_global(4:6)` because those rows represent the angular-momentum balance of the full coupled body-pilus system.

## 8. Brownian Dynamics

Set `brownian.enable` in JSON or `cfg.enable_brown` in MATLAB.

Brownian updates use a body-only Ito approximation with explicit
finite-difference thermal drift `kBT*div(M_eff)`.  `M_eff` is obtained from
the complete coupled body-pili deterministic system, but stochastic
increments are applied only to the body's three translations and three
rotations.  Pili internal Euler-angle modes remain deterministic, and the
root constraints are re-enforced after each body kick.  This is not a full
coupled thermal-equilibrium model and should not be used to study pilus-mode
equipartition.

SPD repair and Cholesky factorisation are performed in dimensionless
generalized coordinates before conversion to SI units.  A final Brownian
step uses the actual remaining time when `T_obs` is not divisible by `dt_B`.
There is no user-facing midpoint drift option.

Important options:

- `brownian.temperature`: positive fluid temperature in Kelvin; required when Brownian dynamics are enabled.
- `brownian.rng_seed`: fixed positive seed for reproducible runs.
- `brownian.dt_B_ratio`: Brownian time step relative to output step.

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
| `flowfield_panels.png` / `.pdf` / optional `.html` | XZ and YZ flow-field panels. |
| `initial_pose.png` / `.pdf` / optional `.html` | Initial-pose preview when enabled. |
| `bacteria_3d_motion.mp4` | Optional 3D animation. |
| `simulation.log` | Complete run log. |
| `brownian_rng_seed.txt` | Brownian seed and reproduce hint when Brownian dynamics are enabled. |

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

The post-processor reads `timeseries.csv` and `config.json`, computes summary kinematic metrics, writes an Excel workbook, and exports a summary figure as PNG, PDF, and optional HTML.

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
