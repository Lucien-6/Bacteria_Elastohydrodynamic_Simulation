# User Guide — Bacteria Elastohydrodynamic Simulation

English | [中文](USER_GUIDE.zh-CN.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Version:** 8.6.6
> **Updated:** 2026-10-21
This guide is for running simulations, choosing the right configuration path, finding outputs, and post-processing completed cases. For model equations and derivations, see [`docs/theory.md`](docs/theory.md). For release history, see [`CHANGELOG.md`](CHANGELOG.md).

Release note: v8.6.6 corrects custom-body surface reconstruction, so a summary figure or animation shows the full body instead of a single facet, sets a 0.5 μm minimum tick spacing on the co-moving windows, and replaces the two body stripes with one short dark-green patch on one side of the waist. The ODE is unchanged and parametric results are unchanged. A custom body whose reconstructed surface changes also changes its enclosed volume, the volume-based length scale, and the sinking force when the cell is not neutrally buoyant. v8.6.5 adds a realistic single-polar-flagellum *E. coli* example (custom geometry, neutrally buoyant, `T_obs = 1 s`); the solver is unchanged and the results of existing cases are unchanged. v8.6.3 refined the animation: the co-moving windows are sized per axis from the reach of the bacterium (no longer cubes) and carry no scale bar, two dark stripes on the body make a roll about the long axis visible, and the four panels are arranged automatically on the canvas; simulation results are unchanged. v8.6.2 gives the motion animation a publication-oriented four-panel layout (3-D and XZ overviews plus co-moving 3-D and XZ windows with ticks and grid, legend, 2-μm scale bars, time stamp with playback speed; the YZ panel is removed); simulation results are unchanged. v8.6.0 speeds up the motion-animation rendering (graphics objects are created once and only updated per frame) and adds `render_animation_from_results` to re-render the animation from saved results; simulation results and video content are unchanged. v8.5.2 calibrates the pilus hydrodynamic radius (`material.hydro_radius_factor`, default `sqrt(2)`; see §4, §7 and `docs/theory.md` §5.2): the axial-twist drag of a straight pilus is now the no-slip-cylinder value instead of half of it and the translational pilus mobility is within about 3 % of slender-body theory, so the default results of pili-bearing cases change slightly (set `hydro_radius_factor = 1` to reproduce 8.5.1); pilus-free runs are unaffected. v8.5.1 re-uses the last exact finite-difference Jacobian across the `ode15s` restarts of Brownian runs, making a Brownian node about 2–3.7× cheaper for multi-pilus cases (about 1.4× for one pilus) within the integration tolerance (see §8; `cfg.jacobian_reuse = false` restores v8.5.0 bit-for-bit); it also cleans up the code-analyzer messages of the run scripts, fixes `verify_migration.m` and adds Markdown lint configuration. v8.5.0 changes how body-only Brownian motion is applied — each interval's thermal increment is realised as a constant generalised force acting through the coupled ODE instead of a displacement kick followed by root re-projection, which removes the spurious pilus root kink at every Brownian node and makes the Ito drift implicit (see §8). Brownian trajectories therefore differ from earlier versions even for the same seed; deterministic runs are unchanged. The type-2 parabolic pilus template is now anchored exactly at the origin (its first node was previously offset by ~4 nm/K). v8.4.0 corrects the `free_torsion` root twist equation — the motor torque is now balanced by the torsional drag of the whole filament instead of the first segment only, so steady rotor speeds of motor-driven runs drop by roughly the segment count $N$ relative to earlier versions (see §7); the pre-flight stall diagnostic uses the matching drag coefficient. v8.3.0 makes the pilus curvature stencil node-collocated and second-order accurate, and turns all 14 parametric attachment-site maps into proper rotations with a right-handed `deflection_x`. Planar pili (types 1-4, `deflection_x = 0`) are placed exactly as before; conical-spiral pili at sites 1, 2, 5, 7, 10, 11, 14 and any pilus with non-zero `deflection_x` now have the intended (previously mirrored) handedness — flip the sign of spiral parameter `c` or of `deflection_x` if a configuration was tuned against the old geometry. v8.2.0 was a physics-correctness update. The flagellar motor uses a signed linear torque-speed law in the rotor speed relative to the body with a documented sign convention (see §7); the flow-field post-processor accounts for the background flow; `shape_params` of pilus types 2-5 are validated up front; and the units of curvature/moment inputs are documented (§4). v8.1.0 extended batch post-processing (`postprocess_cases`): 3-D translational MSD fit on \(0<\tau\le T/2\), plus XZ-plane rotational diffusion \(D_\mathrm{rot}\) from MSAD.

v8.0.1 hardened runtime numerical diagnostics: dense dissipation checks use full `eig`; isolated near-singularity (`rcond < 1e-14`) warns as `NearSingularMatrix`; `SingularMatrix` aborts only under the joint gate of a large backward residual (`> 1e-6`) and near-singular `rcond`. Diagnostic snapshots store `rcond_near_singular_threshold`.

v8.0.0 introduced weighted-reciprocal translation/twist hydrodynamics for unequal segment lengths, equal-chord canonical custom-pilus geometry, objective full-director initialisation, dimensionless Brownian SPD repair, true-duration remainder kicks, isolated Brownian random streams, strict incompressible-shear validation, saved partial timeout returns, and run-scoped numerical diagnostics.

Breaking configuration changes (since 8.0.0): `Pili_Matrix` must contain exactly eight rows; per-pilus `bending_stiffness` is no longer supported; and `shear_z.direction` must be orthogonal to `shear_z.shear_axis`. A segment timeout now saves the achieved trajectory and returns normally with `results.completed=false` and `termination_reason='timeout'`; callers must not use that partial trajectory for steady-state statistics.

## 1. Before You Run

### Requirements

- MATLAB R2020b or later is recommended.
- PNG and PDF exports use `exportgraphics`; optional interactive HTML exports require MATLAB R2026a or later.
- Capsule body surface generation requires Anton Semechko's **S2 Sampling Toolbox** from MATLAB Add-On Explorer / File Exchange.
- No Python, Conda, or pip environment is required. [`requirements.txt`](requirements.txt) documents this explicitly.
- **Parallel Computing Toolbox** is optional: needed for the default parallel branches of [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m) and [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m) (`use_parallel = true`).
- **Running on a Linux HPC cluster (v8.2.0+):** use [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) / [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m) instead of the originals — see §5.3.

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
| Simulate a realistic *E. coli* (one polar flagellum, motor-driven) | Use [`configs/config_ecoli_single_flagellum.json`](configs/config_ecoli_single_flagellum.json); regenerate its geometry with [`geometry/generate_ecoli_geometry.m`](geometry/generate_ecoli_geometry.m). |
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

The `material` section also accepts the optional key `hydro_radius_factor` (default `sqrt(2)`, also `cfg.hydro_radius_factor`): the pilus hydrodynamic regularisation length is `eps_pili = hydro_radius_factor * r_pili`, while the bending stiffness always uses the physical `r_pili`. `sqrt(2)` reproduces the no-slip-cylinder twist drag exactly and the translational mobility within about 3 %; `sqrt(e)` makes the translation exact; `1` reproduces versions up to 8.5.1. Values outside [1, 3] raise a validation warning.

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

### 5.3 Linux HPC cluster batch (`batch_repeat_brownian_HPC_Linux.m`, v8.2.0+)

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

Brownian dynamics is body-only: the 6×6 effective mobility `M_eff` is
obtained from the complete coupled body-pili deterministic system, but
thermal forcing is applied only to the body's three translations and three
rotations.  Pili internal Euler-angle modes receive no thermal noise of their
own.  This is not a full coupled thermal-equilibrium model and should not be
used to study pilus-mode equipartition.

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
| `flowfield_panels.png` / `.pdf` / optional `.html` | XZ and YZ flow-field panels. |
| `initial_pose.png` / `.pdf` / optional `.html` | Initial-pose preview when enabled. |
| `bacteria_3d_motion.mp4` | Optional 3D animation (`bacteria_3d_motion.avi` on Linux). |
| `bacteria_3d_motion_rerender.mp4` | Optional; written by `render_animation_from_results` (see below). |
| `simulation.log` | Complete run log. |
| `brownian_rng_seed.txt` | Brownian seed and reproduce hint when Brownian dynamics are enabled. |

The animation is rendered by `src/visualization/render_motion_animation.m`: all graphics objects are created once and only their data are updated per frame.  It has four panels on a fixed 16:9 canvas (2560 x 1440 logical px; the encoded size follows the display scaling, e.g. 3840 x 2160 at 150 %): **a** 3-D overview (orthographic, fixed limits, gravity arrow), **b** XZ overview (fixed limits, 2-μm scale bar), **c** co-moving 3-D window and **d** co-moving XZ window.  Windows c and d follow the centre of mass at constant magnification; they are not cubes or squares: along each lab axis the box extent is the reach of the body and all pili over all frames (relative to the centre of mass) plus an 8 % margin (at least 0.1 μm), so the bacterium fills the window, and the scale is equal on X, Y and Z (nothing is stretched).  Each window has ticks, tick labels and a grid in lab-frame coordinates (so the grid slides past the bacterium); on panels c and d the spacing is never finer than 0.5 μm.  The four panels are arranged automatically: their aspect ratios and relative sizes are fixed by the data, and the arrangement (nested side-by-side / stacked splits, reading orders a-b-c-d, a-c-b-d, a-c-d-b) that allows the largest common scale factor is chosen and written to `simulation.log` (for example `Layout H(a,H(b,V(c,d)))`), so the canvas stays compact for elongated, flat or compact bacteria alike.  The centre-of-mass trajectory is drawn in a and b only.  The body is uniform light green, with one short dark-green patch on one side of the waist (about 18 % of the body length and ±32° of azimuth). The patch rotates with the body, so a roll about the long axis can be seen in every panel.  The canvas carries a legend, a global time stamp *t* and the playback speed relative to real time (simulated seconds per video second = frame interval x `animation_fps`; labelled "slow motion" or "accelerated"), and panel titles a–d flush on the upper-left corner of each image.  `render_motion_animation` also accepts `opts.scale_bar_um` (default 2) and returns `window_half_um` (box half-extents `[hx hy hz]`), `window_offset_um` and `playback_speed` in its `info` structure.  `simulation.log` ends the animation phase with a per-stage timing table (axis-limit pre-scan, figure setup, data update, `drawnow`, `getframe`, frame resize, `writeVideo`), which shows where the render time goes.

To regenerate the animation from a completed case without re-running the ODE (for example with another frame stride or frame rate), run `addpath('post_processing')` and then `render_animation_from_results()` (file dialog) or `render_animation_from_results(mat_file, fps, skip)`.  Without `fps` / `skip`, the values are taken from the case folder's `config.json` (`output.animation_fps` / `output.animation_skip`), else 30 and 1.  The result is saved next to the `.mat` file as `bacteria_3d_motion_rerender.mp4` (`.avi` on Linux); the original video is never overwritten.

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
| `DissipationEigensolveFailure` on dense mobility | Should not appear after v8.0.1 (dense path uses `eig`). If seen, confirm you are on ≥8.0.1 and that mobility is not unexpectedly sparse. |
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
