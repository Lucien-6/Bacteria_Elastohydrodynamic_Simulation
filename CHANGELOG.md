# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))

---

## [7.4.2] — 2026-05-07

Comment audit and documentation accuracy corrections — **no numerical change** to simulation outputs vs. v7.4.1.

### Fixed

- **`main_coupled_simulation.m`**: output unit convention header now includes `Angular velocity [rad/s]` alongside the existing five units (was previously omitted; `run_simulation.m` already listed it).
- **`main_coupled_simulation.m`**: Legacy `Pili_Matrix` raw-matrix alternative comment now lists per-row semantics (`Row 1: type, Row 2: length, Row 3: deflection Z, Row 4: deflection X, Rows 5–7: shape params, Row 8: site`) consistent with `parse_pili_config.m` and `pili_building.m`.
- **`parse_pili_config.m`**: `Pili_Matrix` output description now mentions the optional row 9 (custom bending stiffness) handled directly in `pili_building.m`, resolving the `[8 x M]` vs `[8-9 x M]` cross-file inconsistency.
- **`build_coupled_lin_sys.m`**: Step 3–4 section comment now explicitly explains why Step 5 (elastic restoring moments) is omitted from this subroutine (it is a RHS term, not part of the LHS; remains in `ode_coupled_rhs.m`).
- **`build_global_resistance.m`**: matrix-layout block annotation corrected — `BP = body-pili (RSM)` and `PB = pili-body (RSS)` (was previously reversed; the column method is determined by the force point per `eps_ij = eps_j`). Added explanatory note on column-method determination.
- **`init_coupled_params.m`**: `alpha_0` default-value comment now notes that the inline demo in `main_coupled_simulation.m` intentionally sets `alpha_0 = 0` (horizontal, non-default) to isolate motor-driven rotation, while the README sinking workflow uses the default `pi/2` (vertical).
- **`read_geometry_csv_mat.m`**: preferred MAT variable-name list now includes `'coordinates'` (was omitted in the comment; code checked five names) and notes case-insensitive matching.
- **`ode_coupled_rhs.m`**: evaluation-step header now explicitly maps Steps 1–4, 6–7 to `build_coupled_lin_sys.m` and Step 5 as not applicable (RHS term), with cross-file reference.
- **`adaptive_basis_check.m`**: the "Householder-like rotation" comment (line 91) replaced with an accurate description: `I + 2*K² = 2*v*v' − I` is a Householder **reflection** (det = −1), not a proper rotation; Walker (2020) uses it as a basis change because director d3 sign flips do not affect the physical curvature binormal.
- **`pili_euler_angles.m`**: added Euler angle convention specification (ZXZ intrinsic rotation: R_z(phi) → R_x(theta) → R_z(psi)) with reference to Walker (2020).

### Added (documentation)

- **`capsule_sphere_mesh_cache.m`**, **`capsule_sphere_cap_points.m`**, **`sphere_particle_sample_cached.m`**, **`snap_sphere_cap_count.m`**: added formal `Syntax / Input / Output` header sections (previously had only brief one-line descriptions).
- **`export_figure_png_pdf_html.m`**: added formal `Description / Syntax / Inputs` header (previously had only a one-line summary and inline syntax call).

### Compatibility

- **Bitwise identical numerical outputs vs. v7.4.1** for identical inputs. All changes are comment/documentation only; no executable code was modified.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**: version **7.4.2**, date **2026-05-07**; comment audit summary.

---

## [7.4.1] — 2026-05-06

Brownian mobility conditioning diagnostic fix — **no numerical change** to `L_eff` / `D_eff` / `M_eff_dim` outputs vs. v7.4.0.

### Fixed

- **`compute_effective_body_brownian.m`**: the `compute_effective_body_brownian:IllConditioned` warning evaluates `rcond(M_eff_dimless)` (after Frobenius symmetrisation of the probe-extracted mobility) instead of `rcond(M_eff_dim)`. The dimensional 6×6 SI mobility mixes TT, TR, and RR blocks scaled by \(\sim L_0^{-2}\)-separated factors from `S * S'`; MATLAB's **1-norm** `rcond` therefore reported spurious “\(\sim 10^{-13}\) near singular” alarms on routine Brownian steps even when the underlying dimensionless mobility was well conditioned. The dimensionless matrix is the intended diagnostic (healthy runs: \(\mathrm{cond} \sim 1..10^{3}\)). Log / `warning` message strings now cite `M_eff_dimless`.

### Compatibility

- **Bitwise identical numerical outputs vs. v7.4.0** for identical inputs; only spurious-console-noise suppression moved. Existing code that parses `warning` IDs still receives `IllConditioned` from this function.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`compute_effective_body_brownian.m`**, **`batch sweep/batch_repeat_brownian.m`**: version **7.4.1**, date **2026-05-06**; user guide §8 / §11 note on the Brownian `rcond` diagnostic.

---

## [7.4.0] — 2026-05-06

Parallel batch execution for independent Brownian repeat runs, plus documentation and site sync.

### Added

- **`batch sweep/batch_repeat_brownian.m` — optional `parfor` parallel batch**:
  - User switches: `use_parallel` (default `true`), `mem_per_worker_gb`, `max_resource_frac` (default `0.75`), `worker_blas_threads`.
  - **Pool sizing:** `n_pool = min(floor(n_phys × max_frac), floor(available_GB × max_frac / mem_per_worker_gb), n_repeats)` with probe source logged as `core_probe` (`wmic`, PowerShell **CIM** `NumberOfCores`, `NUMBER_OF_PROCESSORS/2`, `feature('numCores')/2`, etc.).
  - **Parallel Computing Toolbox** required for the parallel branch; missing license or `parpool` failure **falls back to the legacy serial `for`**.
  - **Brownian seeds in parallel:** `fixed` → `rng_seed_base + k`; `timestamp` → one `brown_rng_seed_from_clock` base + `k` so workers cannot share the same MMddHHmmss seed. Serial `timestamp` path still omits `cfg.brown_rng_seed` per run (unchanged physics intent).
  - **`parallel.pool.DataQueue`** progress lines on the client; summary rows written after the loop in ascending `k` order.
  - Worker init: `DefaultFigureVisible` off + `maxNumCompThreads` cap via `parfevalOnAll` with **`wait`**; defaults persist until `delete(gcp('nocreate'))` (no automatic restore on script exit).
- **`requirements.txt`**: documents **Parallel Computing Toolbox** as **optional** (parallel branch only).

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**: version **7.4.0**, date **2026-05-06**; batch / Brownian-repeat parallel notes.
- **`docs/theory.md`**, **`docs/theory.zh-CN.md`**: header version/date only (no theory change).
- **`docs/index.html`**: hero and release section updated for **v7.4.0**.
- **`run_simulation.m`**, **`main_coupled_simulation.m`**: log / header banner **7.4.0 (2026-05-06)**.
- **`batch sweep/batch_repeat_brownian.m`**: header `Updated` line (2026-05-06) notes parallel scheduling.

### Compatibility

- Single-run **`run_simulation` physics path is unchanged.** Only the batch driver `batch_repeat_brownian.m` gains scheduling code; `use_parallel = false` preserves the pre-v7.4 serial loop behaviour.

---

## [7.3.1] — 2026-05-04

Documentation, help-file, requirements, run-banner, and project website synchronisation release. No numerical algorithm or MATLAB execution path changes relative to v7.3.0.

### Documentation

- **`README.md`**: bumped the project banner to **7.3.1**, added a current-release changelog row, and marked the old v7.2 motor reaction-torque wording as superseded by the v7.3.0 system-angular-momentum audit.
- **`USER_GUIDE.md`**: corrected the documentation-sync release date to **2026-05-04** and kept the motor-torque narrative aligned with the v7.3.0 rule: no explicit reaction term belongs in `R_global(4:6)`.
- **`requirements.txt`**: updated the release banner to **7.3.1 (2026-05-04)** while preserving the MATLAB-only dependency statement.

### Website

- **`docs/index.html`**: refreshed the public project page for the current v7.3.x state, including the hero/version badges, release cards, footer, and bilingual copy for the `free_torsion + motor` physics correction.
- Improved the website wording around v7.3.0/v7.3.1 so the page clearly separates the physics/code correction (v7.3.0) from this documentation/site synchronisation (v7.3.1).

### Compatibility

- Documentation-only release. Simulation results, APIs, saved `.mat` structure, JSON schema, and numerical behaviour are unchanged from v7.3.0.

---

## [7.3.0] — 2026-05-04

Audit-driven physics-correctness fix and a batch of documentation / robustness clean-ups uncovered by the comprehensive review of v7.2.0.  The dominant change reverses the v7.0/v7.2 motor reaction-torque contribution to the system angular-momentum balance — that term was an internal moment that should NOT enter the system equation and was over-prescribing body rotation for `free_torsion + motor` configurations.  Additional fixes cover documentation, validator coverage, code clarity, and a small hot-path vectorisation.

### Fixed

- **`ode_coupled_rhs.m` — remove motor reaction torque from `R_global(4:6)` (Issue #1)**: rows 4-6 of `R_global` are the FULL SYSTEM angular-momentum balance (`B_global` rows 4-6 sum $\mathbf{r}\times\mathbf{f}$ over body + every pilus node with quadrature weights, see `build_coupled_lin_sys.m`). Motor torque between body and pilus is INTERNAL: by Newton's third law the action $+\tau$ on the pilus and the reaction $-\tau$ on the body cancel exactly when summed over the whole body+pilus system, so no contribution belongs in the system moment equation. v7.0.x added $-\tau_{\mathrm{eff}}\mathbf{J}_{\mathrm{euler},3}^{\mathsf{T}}$, v7.2.x added $-\tau_{\mathrm{eff}}\mathbf{d}_{3,\mathrm{lab}}^{(1)}$; both were physically incorrect — they injected an internal torque into the system equation, equivalent to claiming an external moment on the body+pilus that does not exist, and over-prescribed the body's rotation. The body's reaction motion is now propagated implicitly through the global mobility solve: the motor enters the pilus moment-balance row (still present at `R_global(r3) += tau_eff`), drives $d\psi_1/dt$, and through the root kinematic constraints ($\theta_1$, $\phi_1$ locked to body $\boldsymbol{\Omega}$) the body's angular velocity adjusts self-consistently. `geom.J_euler3{k}` is still populated by `build_coupled_lin_sys.m` for diagnostic comparisons and by the updated `tests/test_motor_reaction_torque.m`, but is no longer consumed by the production RHS.
- **`README.md` §11 — corrected $\kappa_{\alpha}$ discretisation formula (Issue #3)**: the displayed equation used $\mathbf{d}_{\alpha}\cdot\mathbf{d}_3'$ for both $\alpha=1,2$, which contradicts the code (`compute_pili_elastic_moments.m` computes $\kappa_1=\mathbf{d}_3\cdot\mathbf{d}_2'$, $\kappa_2=\mathbf{d}_1\cdot\mathbf{d}_3'$) and the underlying Walker (2020) Darboux-vector convention $\mathbf{d}_{\alpha}'=\mathbf{U}\times\mathbf{d}_{\alpha}$. v7.3.0 prints both rows correctly, prefixed with the explicit Darboux-vector definition, and adds a historical note flagging the v7.2 typo.
- **`pili_building.m` — accurate site map for sites 2-5 / 7-10 (Issue #4)**: the per-site comments labelled sites 2-5 and 7-10 as "Left/Right hemisphere side" but the SPoint X-coordinate is `±(major − minor)`, i.e. the cap/cylinder junction circle, not the curved hemisphere surface. Comment now states this explicitly and contrasts sites 11-14 (cylinder midpoint) with the same outward direction but `x = 0`. No code/numerical change.
- **`custom_pili_building.m` — explicit warning on non-zero `deflection_z` / `deflection_x` (Issue #5)**: in custom mode the deflection axis-bindings are SWAPPED relative to parametric `pili_type_building.m` (custom: `deflection_z = TWIST`, `deflection_x = TILT`; parametric: opposite). The semantics were already documented in the file header, but a user copy-pasting parametric values into a custom config would silently get rotated nonsense. v7.3.0 fires a runtime `warning('custom_pili_building:DeflectionAxisSwap', …)` per pilus when either deflection is non-zero, summarising both conventions.
- **`compute_plot_flowfield.m` — sanity check for `Z_final`/`params` shape consistency (Issue #8)**: the function reconstructs pili node velocities using `params.rot{k}`, which may have been rotated mid-simulation by `adaptive_basis_check.m`. The current single-call-at-end usage is safe, but a future intermediate-time call with stale `params.rot` would silently produce incorrect velocities. v7.3.0 asserts `numel(params.rot) == M_pili` and `numel(Z_final) == 7 + 3·sum(N_segs)` upfront with a clear error message.

### Added

- **`validate_simulation_params.m` — stall-curve diagnostic for `free_torsion + constant motor` (Issue #7+#9)**: under the linear stall model $\tau_{\mathrm{eff}}=\tau_{\mathrm{motor}}\left(1-|\omega|/\omega_{\max}\right)$, the steady-state ratio under RFT torsional drag alone is $\omega_{\mathrm{steady}}/\omega_{\max}\approx\tau_{\mathrm{motor}}/(4\pi\mu r_{\mathrm{pili}}^2L_{\mathrm{pili}}\omega_{\max}+\tau_{\mathrm{motor}})$. When this exceeds 0.4 the validator now warns that the explicit-lag `DPSI1_PREV` cache in `ode_coupled_rhs.m` may slow `ode15s` convergence; above 0.7 it strongly recommends raising `max_omega_Hz`, lowering `torque_pN_nm`, or extending the pilus. Available in parametric mode (pili lengths come from `Pili_Matrix(2,:)`); custom mode skips the check because per-pilus length requires reading the user centerline files.

### Changed

- **`coupled_background_flow.m` — vectorised translational and twist-rate blocks (Issue #10)**: replaced the `for i = 1:N_total` triple-write with a single `reshape(u_trans.', [], 1)`; replaced the per-segment `dot(d3_k(j,:), Ω_half)` loop with the matrix product `d3_k * Ω_half'`; in shear_z the inside-band mask is now applied via element-wise multiplication. Numerically equivalent; reduces the inner work in a hot path called once per ODE RHS evaluation. `s_axis_row` and `H_shear` are hoisted out of the pili loop.
- **`ode_coupled_rhs.m` — comment expansion around the removed reaction term**: the v7.2 explanatory comment was retained and extended to spell out (a) why no term should appear in `R_global(4:6)` (system vs internal torque), (b) how the body's reaction motion is propagated implicitly, and (c) the precise legacy expressions and their respective version tags, so the reasoning is self-contained for future readers.
- **`build_coupled_lin_sys.m` — comment update around `geom.J_euler3{k}`**: clarifies that v7.3.0 production code adds NO motor reaction term to `R_global(4:6)`; the field is retained purely for offline diagnostics.
- **`tests/test_motor_reaction_torque.m`**: rewritten as a forensic print-only probe documenting the two legacy candidate vectors (`J_euler(3,:)'` and `d3_lab_seg1`) and the v7.3.0 rule that NEITHER is added to `R_global(4:6)`. The script no longer claims a particular alignment is "correct"; instead it records the body Ω observed under the new rule and notes that anisotropic body drag plus pilus hydrodynamic feedback prevent any simple closed-form alignment prediction.

### Documentation

- **`README.md` §A.3 — re-derived $\ell_k$ factor in $\widetilde{\mathbf{A}}_{\mathrm{twist}}$**: the v7.2 README justified the $\ell_k$ factor only as "ensures dimensional consistency". v7.3.0 spells out the full derivation: per-segment torsional drag $c_t\Delta s_k$ non-dimensionalises to $\varepsilon^2/2$ in the global scale, and the conversion from global scale to the per-rod Walker (2020) scale (used for $\kappa_{\alpha}$, $\mathbf{m}_{\mathrm{active}}$, etc.) introduces an additional $\ell_k$ factor, giving $\ell_k\varepsilon_k^2/2=\varepsilon^2/(2\ell_k)$. The complementary $\ell_k^3$ factor that `build_coupled_lin_sys.m` applies to each pilus moment-balance row converts equations back from the per-rod scale to the global scale, making $\mathbf{L}\dot{\mathbf{z}}=\mathbf{R}$ dimensionally homogeneous (see A.4). The corresponding bullet in §13.4 ($\widetilde{\mathbf{A}}_{\mathrm{twist}}$ contains an extra $\ell_k$ factor) is also expanded.

### Compatibility

- **`free_torsion + motor` configurations**: the simulated body angular velocity, pilus twist rate $d\psi_1/dt$, and any quantities derived from them will differ from v7.0–v7.2 outputs because the spurious motor reaction term in `R_global(4:6)` has been removed. The clamped path is bit-exact unchanged. `free_torsion` without an active motor is also unchanged. Numerical results from prior `free_torsion + motor_torque = constant` simulations should be regenerated.
- All other changes are documentation, comments, validator warnings, or numerically-equivalent vectorisation; no API breaks.

---

## [7.2.0] — 2026-05-04

Audit-driven correctness fixes for the v7.0 free-torsion root coupling + motor torque feature. Three issues uncovered by the first review were addressed in v7.2.0: the then-assumed body-side reaction-torque direction, a silently-disabled motor speed-limiting curve caused by MATLAB's pass-by-value of structs, and an inconsistent ψ₁ initial condition for free-torsion pili at non-trivial body orientations. **Historical correction:** the v7.2.0 reaction-torque interpretation was superseded by v7.3.0, which removes any explicit motor reaction term from `R_global(4:6)` because those rows are the full coupled-system angular-momentum balance. The clamped path is unaffected.

### Fixed

- **`ode_coupled_rhs.m` — motor reaction torque on body (Issue 2, superseded by v7.3.0)**: v7.2.0 replaced the v7.0.x expression $-\tau_{\mathrm{eff}}\mathbf{J}_{\mathrm{euler},3}^{\mathsf{T}}$ with $-\tau_{\mathrm{eff}}\mathbf{d}_{3,\mathrm{lab}}^{(1)}$ on `R_global(4:6)`. That direction correction was useful diagnostically because `J_euler(3,:)` is a kinematic $\omega\to\dot{\psi}$ row, not the physical $\mathbf{d}_3$ vector. However, the later v7.3.0 audit showed that **neither** vector should be added to `R_global(4:6)`: those rows represent the full body+pili system angular-momentum balance, so motor action/reaction torques are internal and cancel in the system sum. Current production code therefore keeps only $+\tau_{\mathrm{eff}}$ on the pilus root $\mathbf{d}_3$ moment-balance row.
- **`ode_coupled_rhs.m` — motor speed-limiting cross-step state (Issue 1)**: The previous-step $d\psi_1/dt$ used by the linear $1-|\omega|/\omega_{\max}$ torque-reduction model was being assigned to `params.dpsi1_prev{k}` at the bottom of the RHS, but MATLAB structs are passed by value: every subsequent ODE call received the original `params` from the lambda capture, so the cache was effectively `[]` forever. As a result `omega_motor = 0` was always read back, `tau_eff = tau_motor` was constant, and the motor never decelerated as it approached `max_omega_Hz`. The cache now lives in a `persistent` cell array (`DPSI1_PREV`) inside `ode_coupled_rhs.m`, gated by a `params.run_token` introduced in `init_coupled_params.m` so cross-simulation contamination cannot occur. `run_simulation.m` additionally issues `clear ode_coupled_rhs` at the start of the integration loop as a belt-and-braces reset, mirroring the existing `clear odetpbar` pattern.
- **`enforce_clamped_geometry.m` / `init_coupled_state.m` — initial $\psi_1$ for free_torsion (Issue 3)**: When `alpha_0` is non-zero (parametric mode) or a custom `initial_orientation` is set, the BC requirement $\mathbf{D}_{\mathrm{lab}}^{(1)}=\mathbf{R}_{\mathrm{body}}\mathbf{D}_{\mathrm{local}}^{(k)}$ cannot be satisfied with $\psi_1=0$; v7.0 nevertheless left $\psi_1$ untouched at initialisation for free_torsion pili because the new conditional projection branch was unconditionally active. This produced a non-zero torsional elastic energy at t = 0 that "released" during the first ODE step, corrupting the early dynamics. `enforce_clamped_geometry` gains an optional third argument `force_full_projection` (default `false`); `init_coupled_state.m` now calls it as `enforce_clamped_geometry(Z0, params, true)` so all three first-segment Euler angles, including $\psi_1$, are projected at t = 0. All runtime callers (`run_simulation.m` post-ODE, post-quaternion-renorm, post-Brownian-kick) keep the default `false`, so free_torsion $\psi_1$ continues to evolve freely after the initial step.
- **`load_sim_config.m` — heterogeneous `pili[]` JSON in parametric mode**: When a `pili` JSON array contained objects with different field sets (e.g. `config_motor_flagella_example.json` where pilus 1 lacks `root_coupling`/`motor_torque` and pilus 2 has them), `jsondecode` returned a cell array of structs rather than a struct array.  The parametric-mode parser only checked `isstruct(jp)`, so the cell array fell through to the empty-pili fallback (`Pili_Matrix = zeros(8, 0)`) — every pilus was silently discarded, producing a body-only initial pose with no flagella.  v7.1 had already shipped `normalise_json_struct_array.m` and used it on the custom path; this fix routes the parametric path through the same helper before calling `parse_pili_config`.  Numeric `[8 x M]` matrices and homogeneous struct arrays continue to work unchanged.

### Added

- **`init_coupled_params.m` — `params.run_token` and `params.dpsi1_prev` placeholder**: a high-resolution monotonic token (`uint64(tic)`) lets `ode_coupled_rhs`'s persistent cache detect a fresh simulation in the same MATLAB session and self-flush. `params.dpsi1_prev = cell(M_pili, 1)` is also pre-allocated so downstream code can reference the field uniformly even though the authoritative cross-step state lives in the persistent variable.

### Changed

- **`ode_coupled_rhs.m` — local variable rename**: the BC-RHS-zeroing loop's offset variable was named `col_offset` but was actually used to index rows of `R_global`; it is now `row_offset_bc` to match the row-based bookkeeping in `build_coupled_lin_sys.m` Step 7. Numerical output is bit-exact unchanged.
- **`build_coupled_lin_sys.m` — comment update around `geom.J_euler3{k}`**: clarifies that the field is populated for diagnostic/tests use only and is no longer consumed by `ode_coupled_rhs.m` for the body reaction torque (see Fixed above).

### Documentation

- **`USER_GUIDE.md`** §7: documents the persistent-cache mechanism for speed-limiting, the then-current v7.2 $-\tau\mathbf{d}_{3,\mathrm{lab}}$ reaction-torque direction, the t = 0 $\psi_1$ projection, and a new **Interaction with `active_moment` and `intrinsic_curvature`** subsection clarifying that for `free_torsion + active_moment.d3_amp != 0 + motor_torque`, all three sources superpose to drive $\psi_1$ at the root. The reaction-torque wording is superseded by the v7.3.0 guide update.
- **`README.md`**: version bumped to **7.2.0** with a new bilingual changelog row summarising the fixes.

### Compatibility

- **Clamped pili path is bit-exact unchanged.** Only configurations containing at least one `root_coupling = 'free_torsion'` pilus are affected.
- **Free-torsion configurations**: v7.2.0 outputs differ from v7.0–v7.1 whenever segment 1's `d₃_lab` and `J_euler(3,:)'` disagree (which is essentially always). v7.3.0 further changes `free_torsion + motor` outputs by removing the explicit system-moment reaction term entirely; regenerate such simulations with the latest release.
- **`config_motor_flagella_example.json` demo timing retuned**: `T_obs` raised from 0.5 s to 1.0 s and `dt` reduced from 0.1 s to 0.001 s so that the motor stall transient and at least a handful of rotations per period are visible in the saved animation (rotation period approximately 3.33 ms at `max_omega_Hz = 300`, so `dt = 1 ms` yields about 3 frames per rotation). Inline `_comment_T_obs` and `_comment_dt` fields were added to the JSON to record the rationale; `load_sim_config.m`/`save_sim_config.m` ignore unknown sibling fields, so the round-trip is unaffected. Existing saved configs that pinned `T_obs = 0.5` / `dt = 0.1` continue to load with their stored values.
- **JSON `pili[]` schema robustness**: configurations that mixed entries with and without `root_coupling`/`motor_torque` (the canonical case in `config_motor_flagella_example.json`) are now parsed correctly in parametric mode. v7.0–v7.1 silently dropped all pili in such configs; users who worked around the issue by replicating the field set across all entries (or by feeding a numeric `[8 x M]` matrix) see no change in behaviour.

---

## [7.1.0] — 2026-05-04

Comprehensive hardening pass on the custom-geometry feature, covering correctness of pili kinematics, robustness of geometry I/O, richer JSON configurability, and the flow-field visualisation that previously assumed a capsule body.

### Added

- **`parse_initial_orientation.m`**: new helper that converts a user-supplied `initial_orientation` spec (`quaternion`, `axis_angle`, or `euler_zyx`) into a scalar-first unit quaternion. Custom mode now accepts arbitrary initial body orientation; when set it **overrides** the legacy `alpha_0` Y-axis-tilt path. `alpha_0` still works in both modes but a warning is printed in custom mode reminding that it assumes body-X is the long axis.
- **`normalise_json_struct_array.m`**: helper that reconciles MATLAB's `jsondecode` cell-vs-struct-array dichotomy. Callers now get a uniform struct array with the union of fields regardless of whether all JSON objects in an array share the same schema.
- **`geometry/` sample files**: pre-shipped `example_body_sphere.csv`, `example_body_capsule.csv`, `example_pilus_linear.csv`, `example_pilus_helix.csv` plus a MATLAB generator `generate_example_geometry.m`. The inline custom block in `main_coupled_simulation.m` and `config_custom_example.json` both reference these files so the example runs out of the box.
- **Pili segment-uniformity check (NEW REQUIREMENT)**: `custom_pili_building.m` and `validate_custom_geometry.m` now enforce that every user-supplied pilus centerline has equal-length segments (relative spread ≤ 1e-3). Non-uniform inputs are rejected with a clear error; this protects the PCHIP resampling and the inextensible Kirchhoff discretisation from ill-posed raw input.
- **Alignment-angle diagnostic (M6)**: `custom_pili_building.m` logs the rotation angle used to align the first segment tangent with `site.direction`, and warns if it exceeds 45°.
- **Attachment-site sanity checks (m7, m8)**: `validate_custom_geometry.m` now flags attachment sites whose position is far from the nearest body point (distance > 5 × mean spacing) and whose direction has a negative projection on the PCA-estimated local outward normal.
- **Per-pilus extended fields in custom mode (M1)**: `intrinsic_curvature`, `active_moment`, `root_coupling`, and `motor_torque` are now read from `custom_geometry.pili[*]`, plumbed through `cfg.pili_ext`, round-tripped by `save_sim_config.m`, and validated symmetrically for both modes via the new `validate_pili_ext` helper.

### Changed

- **`read_geometry_csv_mat.m` rewritten**: primary path now uses `readmatrix` + `detectImportOptions` with `CommentStyle = {'#', '%'}`, honours an explicit `ExpectedNumVariables = 3`, and pre-allocates rows in the manual fallback to eliminate the O(N²) `[all_data; vals']` growth pattern. MAT files prefer variables named `pts` / `points` / `XYZ` / `coords`, warning if a different variable is chosen. 3×3 inputs are now **rejected** (previously a silent warning); N×2 inputs are rejected by default but can be opted into via `'AllowPlanar2D', true`.
- **`custom_body_building.m` hardened (M4, m6, m9, m12)**: input point clouds are deduplicated via `uniquetol` (tol = 1e-6 of bounding-box extent) before the nearest-neighbour-spacing computation, eliminating the `h_min = 0` / RSM-singularity risk. The alpha-shape search was widened to a 15-point log scan over two orders of magnitude around `criticalAlpha`, with a basic non-manifold check rejecting bad alphas. An Euler-characteristic diagnostic warns when the reconstructed surface is not a closed genus-0 manifold. `L_0_suggested` now defaults to the equivalent-sphere radius from the triangulated volume (`L_0_suggest_mode = 'volume'`), with the legacy half-extent formula available via `'extent'`.
- **`custom_pili_building.m` rewritten (M2, M3, M5)**: (i) centerlines are deduplicated before PCHIP resampling; (ii) deflection rotations (`deflection_z`, `deflection_x`) are applied in the pilus **local basis** built via Gram-Schmidt around the aligned root tangent (e3 = aligned tangent, e1 = bending-plane normal), so they act independently of the lab frame instead of the lab frame after alignment; (iii) `Pili(j).root_dir` now stores the **actual** first-segment tangent after all transformations, eliminating the prior semantic mismatch; (iv) the uniform-segment-length check is enforced before resampling; (v) a clear error is raised when the first resampled segment would collapse to zero length. **Caveat on naming**: because the custom local-basis convention binds `deflection_z` to rotation about `e3` (= tangent, a *twist*) and `deflection_x` to rotation about `e1` (perpendicular, a *tilt*), whereas parametric `pili_type_building.m` builds the pilus along `+X` so that parametric `deflection_x` rotates about the *tangent* (twist) and parametric `deflection_z` tilts it, the two parameters' geometric roles are **swapped** between modes. Users porting a non-zero deflection between parametric and custom configurations must swap `deflection_x` ↔ `deflection_z`.
- **`compute_plot_flowfield.m` (C1)**: the previously capsule-hardcoded flow-field renderer now dispatches on `params.body_type` via a new `build_body_visualization` helper. For custom bodies, interior masking uses an `alphaShape`-based `inShape` test (with a bounding-sphere fallback), the 3D body mesh is drawn via `trisurf(params.body_faces, ...)` using the stored triangulation, and 2D silhouettes are extracted via a generic SDF-level contour. All `major_axis` / `minor_axis` references are now confined to the parametric branch.
- **`init_coupled_params.m`**: adds `params.q_0_init`, `params.orientation_source`, `params.body_verts_local`, and persists the body `alphaShape` when available. `init_coupled_state.m` now consumes `params.q_0_init` directly rather than recomputing the alpha_0 quaternion. A new warning fires when `alpha_0 ≠ 0` in custom mode, steering users towards `initial_orientation` for arbitrary bodies.
- **`validate_simulation_params.m`**: the per-pilus extended-field validation was factored into a shared `validate_pili_ext` helper, invoked for both `Pili_Matrix`-driven parametric pili and `custom_geometry.pili`-driven custom pili. The custom branch's Pili_Matrix rejection was tightened to `size(PM,2) > 0` (previously an all-zeros matrix slipped through).
- **`load_sim_config.m`** / **`save_sim_config.m`**: both now understand `custom_geometry.initial_orientation` and the per-pilus extended fields, providing full JSON round-tripping.
- **`main_coupled_simulation.m`**: the broken `F:\…\config_test_custom.json` reference was removed (default is empty = inline mode). The inline custom branch is populated with a working example that mirrors `config_custom_example.json`. The inline block now explicitly lifts `cfg.custom_geometry.initial_orientation` to the top-level `cfg.initial_orientation` so that the custom-mode C2 feature is equally effective in both the JSON path (via `load_sim_config.m`) and the inline path (previously only the JSON path performed this promotion, causing silent fallback to `alpha_0` when users set `initial_orientation` inline).

### Documentation

- **`USER_GUIDE.md`** §2 Option D rewritten to document: sample file locations, the uniform-segment-length requirement, the corrected `deflection_z/x` / `root_dir` semantics, the new `initial_orientation` field, per-pilus extended-field support, the `psi0 = 0` custom-pili assumption, and the flow-field renderer update.
- **`pili_type_building.m`** / **`custom_pili_building.m`**: cross-reference comments added documenting the left-handed sign convention of the `RX` matrix and — critically — the fact that the `deflection_z`/`deflection_x` parameter names bind to **different geometric axes** in the two modes (parametric: tangent = X; custom: tangent = e3). See the in-file comment blocks for the porting rule.

### Compatibility

- Existing parametric configurations are unchanged and bit-exact.
- Existing custom configurations may now **fail validation** if they ship non-uniform pili centerlines — this is intentional and reflects the new correctness requirement. Pre-discretise the centerline with equal-length segments to resolve.
- `deflection_x` / `deflection_z` behaviour change for custom mode: v6.x–v7.0.x applied the deflection rotations in the **lab frame** after site alignment; v7.1 moves them into the pilus **local basis** (e1, e2, e3). Users who tuned non-zero deflections against v7.0 custom configurations must re-tune after upgrading. Note: the v7.1 local-basis convention binds `deflection_x` to a tilt about `e1` (perpendicular to the tangent) and `deflection_z` to a twist about `e3` (the tangent). This is the OPPOSITE of the parametric-mode axis assignment; see the `pili_type_building.m` / `custom_pili_building.m` headers.

---

## [7.0.3] — 2026-04-30

### Fixed

- **`read_geometry_csv_mat.m` — CSV comma-delimiter parsing**: The local `read_csv_points` helper used `sscanf(line, '%f')` to parse numeric values from each line. MATLAB's `sscanf` with `%f` only treats whitespace (space, tab, newline) as field delimiters — commas are not recognised, causing comma-separated CSV data (e.g. `1.0e-6,0.0,3.98e-7`) to be silently truncated after the first number. Each subsequent data line yielded fewer than 3 values and was skipped, resulting in a "Point data must have 3 columns, got 0" error for all CSV geometry files. **Fix**: before header detection and numeric parsing, all data lines are normalised by replacing commas (`,`) and tabs (`\t`) with spaces, so `sscanf('%f')` correctly extracts every numerical field regardless of the original delimiter style. This makes the function accept comma-separated, space-separated, tab-separated, and mixed-delimiter CSV files without changing the output for any previously working format.

---

## [7.0.2] — 2026-04-30

### Fixed

- **`adaptive_basis_check.m` — `sprintf` string syntax error (line 114)**: In a prior edit, invisible Unicode characters (zero-width spaces or non-breaking spaces) were introduced around the line-continuation `...` operators in the `InsufficientRotation` warning block, causing MATLAB's parser to emit `'意外的 ''Post-rotation min(theta, pi-theta)...''。'` (unexpected string end). The affected `logf_h` and `warning` calls were rewritten using explicit `sprintf` string construction to eliminate multi-line string concatenation and the associated hidden-character issue.
- **Warnings restored to MATLAB default display format**: In v6.1.0, duplicate `warning()` calls were removed from `validate_simulation_params.m` and `init_coupled_state.m` to eliminate redundant console output. However, this caused all warning messages to appear in black (plain `fprintf`) rather than MATLAB's built-in orange warning color. The native `warning()` calls are now restored alongside the existing `logf_h` calls, preserving both the colored console output and the log-file record.

### Changed

- **`init_coupled_params.m` — log output cleanup**: Removed the redundant "Print summary" parameter block (~90 lines) and the "Computed adaptive gap distance" log lines. All derived parameter information is now centrally reported in `run_simulation.m`'s `[INIT]` section, eliminating duplicate log entries.
- **`main_coupled_simulation.m` — log header and `[PARAMS]` deduplication**: Trimmed duplicate log banners and redundant `[PARAMS]` blocks so the console log is clean and non-repetitive.

### Documentation

- **`README.md`**, **`CHANGELOG.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**: Version **7.0.2**, date **2026-04-30**.

---

## [7.0.1] — 2026-04-30

### Fixed

- **`pili_type_building.m` — Type 4 & 5 arc-length uniform resampling**: Types 4 (expansive sinusoidal) and 5 (conical spiral) previously sampled nodes at equal intervals of the curve parameter (x or t), producing **non-uniform arc lengths** between consecutive nodes. This violated the constant-segment-length assumption (`ds = L/N_seg`) of the inextensible Kirchhoff rod framework, causing cumulative reconstruction errors (up to 0.47 L₀ in Type 5). Both types now use **PCHIP-based arc-length uniform resampling**: a fine parameter sweep (Δ = 1e-4 over [0, 6]) computes cumulative arc length at sub-nanometer resolution; PCHIP interpolation then directly maps uniform arc-length targets to each coordinate component, bypassing the indirect arc-length → parameter → coordinate path. The sampling strategy ensures all `Ntail` nodes are equally spaced in arc length (`ds = L_helix / (Ntail-1)`), matching `pili_spatial_coords`'s reconstruction assumption (`ds = ell_k / N_k`). Types 1, 2, 3 remain unchanged (naturally arc-length uniform or explicitly arc-length-parameterised). The function header now documents the sampling strategy per type. An error is raised if the fine sweep's maximum arc length is shorter than the target `L_helix`. Reconstruction error for `[0.1, 0.5, 2π]` conical spiral reduced from ~0.47 to ~0.0135 L₀ (remaining error is due to PCHIP smoothing of the coordinate data; further reduction requires coordinate-level arc-length enforcement).

### Changed

- **`pili_type_building.m` — search range extended**: Fine parameter sweep range expanded from `[0, 4]` to `[0, 6]` (both Type 4 and 5), ensuring longer pili are fully covered.
- **`pili_type_building.m` — arc-length target**: Uniform sampling now uses the exact user-specified `L_helix_m` as the arc-length target (not the nearest discrete match `L_actual`), ensuring `ds = L/N_seg` is consistent with `pili_building.m`'s `ell = L/(N-1)`.

### Compatibility

- **Backward-compatible for Types 1, 2, 3**: Code paths unchanged; bit-exact results for existing configurations using these types.
- **Geometry change for Types 4, 5**: Node positions shift from parameter-uniform to arc-length-uniform distribution. The curve shape itself is unchanged — only the point spacing along the curve changes. Physical quantities (curvature, forces, moments) become more accurate because the constant-ds assumption is now satisfied.

---

## [7.0.0] — 2026-04-29

### Added

- **Free-torsion root coupling (`root_coupling = 'free_torsion'`)**: Per-pilus root boundary condition now supports two modes — `'clamped'` (default, all three Euler angles $\theta_1$, $\phi_1$, $\psi_1$ kinematically locked to the body) and `'free_torsion'` ($\theta_1$ and $\phi_1$ locked, but $\psi_1$ free to rotate under its own moment balance). The free-torsion mode models a "clamped-torsion-free" boundary condition suitable for motor-driven rotating flagella.
- **Motor torque (`motor_torque`)**: Per-pilus active torsional moment applied at the root segment. Currently supports `'off'` (default, no motor) and `'constant'` (user-specified torque magnitude in pN nm with optional stall-speed limit). For `'constant'` mode, the motor torque is nondimensionalised at initialisation; the effective torque is linearly reduced as the current $d\psi_1/dt$ approaches `max_omega_Hz` (speed-limiting model). Historical note: v7.0.0 attempted to apply a body-side reaction torque via `J_euler(3,:)`; this was first revised in v7.2.0 and then fully superseded by v7.3.0, where no explicit motor reaction term is added to the full-system angular-momentum rows.
- **`parse_pili_config.m` — `root_coupling` / `motor_torque` fields**: The `pili_ext` struct now includes `root_coupling` (default `'clamped'`) and `motor_torque` (default `struct('type', 'off')`) per pilus, parsed from user-friendly pili configuration structs.
- **`init_coupled_params.m` — `root_coupling` / `motor_torque` Name-Value parameters**: Two new `inputParser` parameters accept cell arrays per pilus. Nondimensionalisation converts `torque_pN_nm` → `torque_dimless` (scaled by EI/L₀) and `max_omega_Hz` → `max_omega_dimless` (scaled by 2π·T₀). Parameter summary now prints per-pilus `root_coupling` type and `motor_torque` details.
- **`validate_simulation_params.m` — root coupling & motor torque validation**: New rules ensure `root_coupling` ∈ {`'clamped'`, `'free_torsion'`}; validate `motor_torque.type` ∈ {`'off'`, `'constant'`}; for `'constant'` check `torque_pN_nm` finite non-zero and `max_omega_Hz` > 0; warn if motor torque is specified for a clamped pilus (no physical effect).
- **`build_coupled_lin_sys.m` — conditional root BC**: Root boundary condition assembly (Step 7) now branches on `params.root_coupling{k}`. `'clamped'` retains the original 3-row replacement ($\theta_1$, $\phi_1$, $\psi_1$). `'free_torsion'` replaces only 2 rows ($\theta_1$, $\phi_1$); the $\psi_1$ row is left as the natural moment-balance equation. `J_euler(3,:)` is stored in `geom.J_euler3{k}` for downstream reaction-torque computation.
- **`compute_pili_elastic_moments.m` — conditional $\kappa_3(1)$ zeroing**: For `'clamped'`, `kappa3_full(1)` is zeroed (overwritten by BC). For `'free_torsion'`, `kappa3_full(1)` is preserved — the torsional elastic moment contributes to the $\psi_1$ dynamics equation. New optional 7th argument `root_coupling_k` (default `'clamped'`).
- **`coupled_active_moment.m` — conditional M_active(:,1) zeroing**: For `'clamped'`, all three components ($\mathbf{d}_1$, $\mathbf{d}_2$, $\mathbf{d}_3$) of `M_active(:,1)` are zeroed. For `'free_torsion'`, only bending components ($\mathbf{d}_1$, $\mathbf{d}_2$) are zeroed; the torsional component $\mathbf{d}_3$ is preserved and drives free $\psi_1$ dynamics. New optional 11th argument `root_coupling_k` (default `'clamped'`).
- **`ode_coupled_rhs.m` — motor torque application**: New Step 8d applies motor torque for free-torsion pili. Calculates `tau_eff` from `params.motor_cfg{k}.torque_dimless` with speed-limiting (linear reduction if `abs(omega_motor)` approaches `max_omega_dimless`). Applies `tau_eff` to `R_global(r3)` (pilus torsional force). Historical note: v7.0.0 also applied `J_euler3{k}' * tau_eff` to `R_global(4:6)` as a body angular moment; v7.3.0 removes that pattern because `R_global(4:6)` is the full body+pili system moment balance and must not receive internal motor torques. Root BC RHS zeros: `R_global(r1)` and `R_global(r1+1)` always zeroed; `R_global(r1+2)` only zeroed for `'clamped'` pili. After the linear solve, `z_dot(psi1_col)` is stored in `params.dpsi1_prev{k}` for the next time step's speed-limiting estimate.
- **`enforce_clamped_geometry.m` — conditional $\psi_1$ projection**: For `'clamped'`, all three root-segment Euler angles ($\theta_c$, $\phi_c$, $\psi_c$) are projected to match the body orientation. For `'free_torsion'`, only $\theta_c$ and $\phi_c$ are projected; $\psi_c$ is left untouched (free to evolve dynamically).
- **`save_sim_config.m` — `root_coupling` / `motor_torque` persistence**: Per-pilus `root_coupling` and `motor_torque` fields are now written to the JSON output (stripping non-serialisable `fun` handles from motor configs).
- **`run_simulation.m` — root_coupling / motor_torque extraction**: Extracts `root_coupling_cfg` and `motor_torque_cfg` cell arrays from `cfg.pili_ext` and passes them as Name-Value pairs to `init_coupled_params`, following the same pattern as `intrinsic_curvature_cfg` / `active_moment_cfg`.
- **`main_coupled_simulation.m` — inline examples**: Commented-out inline configuration examples for `root_coupling` and `motor_torque` per pilus. Updated `pili_ext` default struct to include new fields.
- **`config_motor_flagella_example.json`** (NEW): Dedicated example configuration for motor-driven flagella with free-torsion coupling — one clamped pilus + one free-torsion pilus with constant motor torque (2000 pN·nm, stall speed 300 Hz).

### Changed

- **§13 "固支边界条件" → §13 "根部边界条件"**: The README section title and content are updated to reflect that the root boundary condition is no longer exclusively clamped — free-torsion is now a supported alternative per pilus.
- **§15 "时间积分"**: The algorithm pseudo-code Step 9 is updated from "施加固支边界条件" to "施加根部边界条件"; Steps 7b and 8 are extended with motor torque logic.

### Compatibility

- **Fully backward-compatible**: `root_coupling` defaults to `'clamped'` and `motor_torque` defaults to `struct('type', 'off')` for all pili when not specified. Existing configurations (inline or JSON) produce identical results without any modification. The `'clamped'` code path is bit-exact with the pre-v7.0 code at all checkpoints — only the branching logic was added, not the arithmetic.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**, **`CHANGELOG.md`**: Version **7.0.0**, date **2026-04-29**.

---

## [6.1.0] — 2026-04-29

### Added

- **Dual-channel logger (`logf` / `logf_file`)**: `run_simulation.m` now creates two logger closures — `logf` (console + file, `dual_log`) and `logf_file` (file only, `log_file_only`). Both write `[HH:mm:ss]` timestamps to the log file; console output remains timestamp-free for readability. The `logf` and `logf_file` handles are passed into `params` via `init_coupled_params` so sub-functions can write to `simulation.log` without needing direct access to the file handle.
- **`init_coupled_params.m` — `logf` / `logf_file` Name-Value parameters**: Two new `inputParser` parameters with `fprintf` / no-op fallbacks allow the logger handles to be injected from `run_simulation` while keeping standalone / batch calls working without modification. Handles are stored as `params.logf` / `params.logf_file` and a local shorthand `logf_h` is used throughout the function body.
- **`init_coupled_state.m` — `logf_h` parameter**: Optional second argument (function handle, `fprintf` fallback). All `fprintf` and `warning` calls replaced with `logf_h` for log-file capture while retaining `warning` for MATLAB's console alert system.
- **`validate_simulation_params.m` — `logf_h` parameter**: Optional second argument (no-op fallback). Warning loop now writes to both `logf_h` (log file) and MATLAB `warning` (console).
- **`build_coupled_lin_sys.m` — `get_params_logf` helper**: Euler-angle singularity `warning` now also writes to the log file via `logf_h = get_params_logf(params)` before calling MATLAB `warning`.
- **`compute_effective_body_brownian.m` — `get_brown_logf` helper**: Three diagnostic warnings (`AsymmetryHigh`, `IllConditioned`, `NotPosDef`) now also write to the log file via `logf_h` before calling MATLAB `warning`.
- **`adaptive_basis_check.m` — `get_basis_logf` helper**: `InsufficientRotation` warning now also writes to the log file via `logf_h` before calling MATLAB `warning`.
- **Animation progress logging**: `run_simulation.m` logs animation frame progress to `simulation.log` every 10% (file-only via `logf_file`, no console echo), complementing the existing console-only `render_progress_bar`.

### Changed

- **Default `E` (Young's modulus)**: Unified to **`3.279e8` Pa** across `init_coupled_params.m`, `load_sim_config.m`, and `main_coupled_simulation.m` (previously `3.67e8` in `init_coupled_params` and `1.768e7` in `load_sim_config`).
- **Default `N_body`**: Unified to **300** across `init_coupled_params.m`, `load_sim_config.m`, `main_coupled_simulation.m`, and `config_example.json` (previously 200 in `load_sim_config` / `config_example.json` and 300 in `init_coupled_params`).
- **`init_coupled_params.m`**: All `fprintf` calls replaced with `logf_h` so parameter summary and diagnostic messages are captured in `simulation.log`. Documentation comments for `E` and `N_body` defaults updated to match new values.
- **`run_simulation.m`**: `logf` / `logf_file` passed as Name-Value pairs to `init_coupled_params` and `logf` passed to `init_coupled_state` / `validate_simulation_params`. Before saving `.mat`, `params.logf` and `params.logf_file` are removed via `rmfield` (function handles are not serialisable).

### Fixed

- **`config_example.json` — `N_body`**: Updated from `200` to `300` to match the code default, so users who run directly from the example JSON get the standard body discretisation.
- **`config_example.json` / `config_custom_example.json` — `compute_flowfield`**: Both JSON templates now include `"compute_flowfield": true` in the `output` section, matching the code default and ensuring reproducibility when reloading a saved config.
- **Log-file completeness**: Previously, warnings and diagnostic messages from `init_coupled_params`, `init_coupled_state`, `build_coupled_lin_sys`, `compute_effective_body_brownian`, `adaptive_basis_check`, and `validate_simulation_params` were only visible on the MATLAB console (via `fprintf` / `warning`) and never written to `simulation.log`. All such messages are now captured in the log file with `[HH:mm:ss]` timestamps.

### Compatibility

- **Fully backward-compatible**: All logger handle parameters are optional with safe fallbacks (`fprintf` or no-op). Standalone calls to `init_coupled_params`, `init_coupled_state`, or `validate_simulation_params` without logger arguments work identically to before. Batch sweep scripts are unaffected.
- No physics, ODE right-hand side, linear-system assembly, Brownian kernel, or time-stepping path is modified. Trajectory outputs (`results.mat`, `timeseries.csv`) are bit-exact for identical configurations.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**, **`CHANGELOG.md`**: Version **6.1.0**, date **2026-04-29**.

---

## [6.0.0] — 2026-04-29

### Added

- **Custom geometry mode** (`geometry_mode = 'custom'`): Users can now import arbitrary 3D point clouds for both the bacterial body and its pili, moving beyond the existing parameterized capsule/analytic-pili models. Custom body and custom pili **must** be used together — mixing with parametric models is explicitly disallowed and validated at startup.
- **`read_geometry_csv_mat.m`** (NEW): Utility to read 3D point cloud data from `.csv`, `.txt`, or `.mat` files. Handles header rows, comment lines (`%`/`#`), and auto-detects matrix orientation (3×N vs N×3). Issues a warning for ambiguous 3×3 matrices (assumes N×3 = 3 points). Returns a validated N×3 matrix.
- **`rotation_between_vectors.m`** (NEW): Implements Rodrigues' rotation formula to compute the minimal 3×3 rotation matrix that aligns one 3D vector to another. Includes special handling for parallel and anti-parallel cases. Used by `custom_pili_building.m` to orient pili at attachment sites.
- **`custom_body_building.m`** (NEW): Processes user-provided custom body point clouds. Reads data via `read_geometry_csv_mat`, computes characteristic point spacing (`h_mean`), automatically triangulates the surface for visualization using `alphaShape` (with `convhulln` fallback) via a local `auto_triangulate_surface` helper, calculates the enclosed volume, and processes user-defined attachment sites. Outputs a `body_meta` struct with volume, point spacing, triangulation, and attachment site data.
- **`custom_pili_building.m`** (NEW): Constructs custom pili from user-provided centerline points. Reads data, calculates true arc length, resamples to the target segment density using PCHIP interpolation, aligns the pilus root tangent to the specified attachment site direction via `rotation_between_vectors`, applies optional `deflection_z`/`deflection_x` additional rotations, and translates to the body-local attachment point. Outputs a `Pili` struct formatted identically to `pili_building.m` for downstream compatibility.
- **`validate_custom_geometry.m`** (NEW): Dedicated validation for custom geometry configurations. Checks file existence and integrity, verifies attachment site consistency (non-zero directions, unique labels), and ensures pili site indices are valid attachment site references.
- **`config_custom_example.json`** (NEW): Full JSON configuration template for custom geometry mode, including `body_points_file`, `attachment_sites` (position + direction + label), and `pili` definitions (points_file + site_index + deflection angles).
- **`init_coupled_params.m`**: Added `geometry_mode` and `custom_geometry` parameters to `inputParser`. Introduced `switch` statements to dispatch body construction to `capsule_body_building` (parametric) or `custom_body_building` (custom), and pili construction to `pili_building` or `custom_pili_building`. Dynamically determines `L_0` and `V_body` based on the chosen mode. Stores geometry-specific metadata (`body_type`, `body_faces`, `body_verts`, `custom_attachment_sites`). Added a guard to error if `V_body` is `NaN` or non-positive (e.g., due to triangulation failure).
- **`load_sim_config.m`**: Reads the `geometry_mode` field. For `'custom'` mode, parses the `custom_geometry` node and explicitly sets parametric-specific fields (`minor_axis`, `major_axis`, `N_body`, `Pili_Matrix`) to `NaN` or empty to enforce the "no mixing" constraint. `alpha_0` and `eps_body` for custom mode are read from the `custom_geometry` node if present and promoted to the top-level `cfg` struct. `N_segments` derivation is now conditional on `geometry_mode`.
- **`save_sim_config.m`**: Writes the `geometry_mode` field and uses a `switch` statement to correctly save either the `body`/`pili` nodes (for parametric mode) or the `custom_geometry` node (for custom mode). For custom mode, `alpha_0` and `eps_body` are moved back into the `custom_geometry` sub-struct for consistency.
- **`validate_simulation_params.m`**: Added a top-level `geometry_mode` switch. For `'custom'` mode, calls `validate_custom_geometry.m` and includes explicit checks to error if any parametric-specific fields (`Pili_Matrix`, `minor_axis`, `major_axis`, `N_body`) are set. Common validations apply regardless of geometry mode. The Brownian Fixman midpoint warning now correctly checks for custom pili. The Stokes-Einstein `dt_B` warning applies only in parametric mode where `major_axis` is defined.
- **`run_simulation.m`**: Unpacks `geometry_mode` and `custom_geometry_cfg` from `cfg`. Parameter validation passes these new fields. Logging conditionally prints geometry details. Visualization routines (initial pose preview, summary plots, and animation) use `trisurf()` with the automatically generated mesh for custom bodies, replacing `surf()` for parametric capsules. Axis limit calculations adapt for both geometry types.
- **`main_coupled_simulation.m`**: Added comments explaining the `geometry_mode` options and a commented-out inline configuration example for custom mode. Inline parameter structure refined with an `if/else` block for parametric vs. custom-specific configuration, with common parameters defined outside the block. Explicit initialization of `cfg.alpha_0`, `cfg.Pili_Matrix`, and `cfg.pili_ext` in the custom branch ensures these top-level fields are always defined.

### Changed

- **"No mixing" enforcement**: When `geometry_mode = 'custom'`, any parametric-specific configuration fields (`Pili_Matrix`, `minor_axis`, `major_axis`, `N_body`) are explicitly rejected by `validate_simulation_params.m` with clear error messages.
- **`custom_pili_building.m`**: Removed redundant `kb`/`kb0` (curvature binormal) computation since `init_coupled_params` recomputes all curvature information from scratch. These fields are now initialized to `zeros(3,0)` as placeholders.
- **Custom body visualization**: The automatic surface triangulation (via `alphaShape` + `convhulln` fallback) generates `body_faces`/`body_verts` stored in `body_meta`, enabling `trisurf` rendering in `run_simulation.m` without requiring a separate mesh file.

### Fixed

- **Inline custom mode**: Missing `cfg.alpha_0`, `cfg.Pili_Matrix`, and `cfg.pili_ext` fields in the custom branch of `main_coupled_simulation.m` caused runtime "Reference to non-existent field" errors. All three fields are now explicitly initialized.
- **Custom body volume `NaN` propagation**: If `custom_body_building.m` failed to triangulate, `V_body = NaN` would silently propagate to physical parameters. `init_coupled_params.m` now errors explicitly on invalid volume.
- **Fixman midpoint warning in custom mode**: The warning for `enable_brown + pili + fixman_midpoint=false` was suppressed in custom mode because `cfg.Pili_Matrix` is `zeros(8,0)`. The check now also considers `cfg.custom_geometry.pili`.
- **Custom mode log output**: `run_simulation.m` no longer prints `NaN` for parametric-specific fields (`major_axis`, `minor_axis`, `N_body`) in custom mode; logging is conditional on `geometry_mode`.
- **`validate_custom_geometry.m` documentation**: Removed the misleading "volume plausibility" claim from the header (the function validates configuration only, not derived geometric properties).
- **`read_geometry_csv_mat.m` 3×3 ambiguity**: Refined the MAT file reading logic to specifically distinguish 3×N (N>3) matrices from ambiguous 3×3 matrices (issues a warning, assumes N×3 = 3 points).

### Compatibility

- **Parametric mode is fully backward-compatible**: All existing inline configurations, JSON configurations, and batch sweep scripts produce identical results without any modification. The `geometry_mode` defaults to `'parametric'` when not specified, and all new code paths are gated behind this switch.
- Custom geometry mode is a self-contained addition; no existing function signatures or output formats are changed.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**, **`CHANGELOG.md`**: Version **6.0.0**, date **2026-04-29**; README feature-history row **79**; file structure and quick-start sections updated with custom geometry files and configuration; root `.m` file count **46** (total **60**).

---

## [5.5.8] — 2026-04-22

### Added

- **`export_figure_png_pdf_html.m`**: Centralised figure export — **`exportgraphics`** saves **PNG** (600 dpi), **PDF** (`ContentType`, `vector`), and on **MATLAB R2026a or later** an **interactive HTML** web canvas (same base path, `.html` extension). Used by `run_simulation.m`, `compute_plot_flowfield.m`, and `post_processing/postprocess_cases.m` (`export_batch_summary_figure`).
- **Vector PDF alongside raster PNG** for `initial_pose`, `summary_plot`, `flowfield_panels`, and the post-process batch summary figure (same stem as the workbook).
- **Post-process batch summary figure** (`postprocess_cases.m`): 2×2 bar-chart panel (kinematic metrics per case) saved next to the chosen `.xlsx` as `.png` / `.pdf` / optional `.html` (R2026a+).

### Changed

- **`main_coupled_simulation.m` / `run_simulation.m`**: Documented path policy — `addpath(project root)` at the main script and defensively at the start of `run_simulation` so programmatic calls resolve shared helpers; `post_processing/` and `batch sweep/` keep local `addpath` for standalone runs.
- **Completion summaries list HTML only when present**: `run_simulation.m` *Files saved* block appends `*.html` only if `isfile`; `postprocess_cases.m` console and `questdlg` include the HTML path only when the file exists.
- **`run_simulation.m`**: Output-directory data whitelist extended with **`.pdf`** and **`.html`** (overwrite / rename probe).
- **`batch sweep/batch_sweep_eps_body_RSM.m`**, **`batch_sweep_N_body_eps_RSM.m`**: Header notes cross-reference the path policy in `main_coupled_simulation.m`.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**, **`CHANGELOG.md`**: Version **5.5.8**, date **2026-04-22**; README feature-history row **78**; §17 output table and quick-start tree updated for PDF/HTML and `export_figure_png_pdf_html.m`; file counts **55** `.m` files (root **41**).

### Numerical equivalence

- Export and logging/path changes only. No ODE, linear-system, or Brownian physics path is modified. Trajectory outputs (`results`, `timeseries.csv`) unchanged for identical configs.

---

## [5.5.7] — 2026-04-21

### Added

- **`run_simulation.m` — enriched `[PARAMS]` section**: Per-pilus logging now prints `intrinsic_curvature.type` and `active_moment.type` for every pilus; `traveling_wave` active moments additionally emit `amplitude_Nm`, `frequency_Hz`, and `wavelength`. A new `[PARAMS] Background flow:` block is appended at the end of the parameters section with three branches — `off`, `uniform` (logs `U_max`, `direction`), and `shear_z` (logs `U_max`, `direction`, `shear_axis`, `shear_height_m`). Two local helpers (`fmt_field_value`, `get_field_default`) format struct fields defensively (`'(n/a)'` fallback for missing/empty entries) so the log remains well-formed even if `cfg.pili_ext` or `cfg.background_flow` is only partially specified.
- **`run_simulation.m` — completed `Files saved` manifest**: The final "Files saved" block at the end of the run now lists `Config : config.json` and, only when `enable_brown = true`, `Seed : brownian_rng_seed.txt`, alongside the existing `Data / CSV / Plot / Pose / Video / Flow / Log` entries; column labels are re-aligned to a uniform width. This removes the mismatch between the human-readable summary and the files actually written into `output/<case>/`.
- **`save_sim_config.m` / `load_sim_config.m` — `compute_flowfield` persistence**: The flow-field on/off toggle is now saved into `config.json` under `output.compute_flowfield` (default `true`, matching `run_simulation.m`'s backward-compatible fallback) and read back by `load_sim_config` as `cfg.compute_flowfield`. Reloading a saved configuration no longer silently flips flow-field computation back on; reproducibility of batch reruns and JSON-driven sweeps is now complete for this flag. Header comment in `load_sim_config.m` is updated to list the new field.
- **`batch sweep/` — unified summary logs + per-run exception isolation (8 scripts)**: Every sweep driver (`batch_sweep_alpha`, `batch_sweep_N_body`, `batch_sweep_major_axis`, `batch_sweep_pili_E`, `batch_sweep_pili_length`, `batch_sweep_seg_per_um`, `batch_sweep_eps_body_RSM`, `batch_sweep_N_body_eps_RSM`) now wraps the `run_simulation(cfg_i)` call in a `try/catch` block that records the per-run status (`OK` / `FAILED` / `ABORTED`), wall-clock duration, and — on failure — the exception message, without aborting the enclosing sweep. A persistent timestamped summary log `batch sweep/<script>_summary_<ts>.log` is opened at the start of each sweep and receives (i) a header block (script name, sweep variable, sweep grid, start time), (ii) one status line per run, and (iii) a footer with total elapsed time, success/failure/abort counts, and paths to any CSV/MAT products. The two RSM scripts additionally persist `run_status` / `run_walltime` / `run_errmsg` (or `status_grid` / `walltime_grid` / `errmsg_grid` for the nested sweep) inside their `.mat` outputs. The stylistic convention and summary fields deliberately mirror the existing `batch_run_from_json.m` and `batch_repeat_brownian.m` logs for operator familiarity.

### Changed

- **`run_simulation.m` — overwrite detection whitelist**: The `data_exts` list inside the output-directory overwrite probe now includes `.json` and `.txt` (previously only `.mat/.csv/.png/.mp4/.log`). `config.json` and `brownian_rng_seed.txt` are therefore properly recognised as pre-existing run artefacts instead of being silently left behind when the user chooses `Overwrite`.
- **`run_simulation.m` — empty `case_name` log**: When `cfg.case_name` is blank, the `[PARAMS]` banner now prints the auto-generated `run-<MMddHHmmss>` identifier (which was already resolved by that point) instead of the placeholder `(auto — timestamp)`. The two-line banner behaviour is unchanged for non-empty `case_name`.
- **`post_processing/postprocess_cases.m` — `datestr` → `datetime`**: Deprecated `datestr(now, …)` calls replaced with `datetime('now', 'Format', …)` (two sites: the console banner timestamp and the default `.xlsx` filename date tag). Output identical; avoids upcoming `datestr` deprecation warnings in newer MATLAB releases.
- **`batch sweep/batch_sweep_N_body_eps_RSM.m` — syntax fix**: Stray leading character `d` on line 1 (`d%% batch_sweep_N_body_eps_RSM.m` → `%% batch_sweep_N_body_eps_RSM.m`) removed so the file parses and runs. The script had been unrunnable since it was introduced.

### Numerical equivalence

- Purely diagnostic, log-facing, and batch-plumbing changes. No physics, ODE right-hand side, linear-system assembly, Brownian kernel, or time-stepping path is touched. Re-running any pre-5.5.7 case (inline cfg or JSON) yields identical `results.mat` and `timeseries.csv` contents; the `.log`, `.json`, and newly-added `batch sweep/*_summary_*.log` files are richer but strictly additive.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**, **`CHANGELOG.md`**: Version **5.5.7**, date **2026-04-21**; README feature-history row **77** added covering the logging polish, sweep-summary unification, and `compute_flowfield` persistence; USER_GUIDE v5.5.7 paragraph added inline with the rest of the 5.5.x history.

---

## [5.5.6] — 2026-04-21

### Added

- **`Brownian_Implementation_Plan.md` (tracked companion)**: New repository document with **§1.5 — Scope & Simplifying Assumptions** (modelled body 6 DOF + pili-induced drag through the coupled mobility; deliberate omission of pili internal thermal modes; extension references) and **§2.4.1 — Stochastic Interpretation & Missing Drift** (Itô form, missing drift when `fixman_midpoint=false`, Fixman recovery, decision matrix). The file was previously listed in `.gitignore` and therefore missing from a clean clone; it is now versioned so cross-references in `validate_simulation_params.m` resolve.
- **`validate_simulation_params.m`**: Two new Brownian-related warnings. (a) When `enable_brown = true` and `Pili_Matrix` is non-empty but `fixman_midpoint = false`, a warning is emitted recommending Fixman to avoid the missing Itô drift. (b) When `dt_B` exceeds 10% of a Stokes-Einstein back-of-envelope estimate of `min(tau_trans, tau_rot)` (using the major semi-axis as a conservative radius), a warning suggests a smaller `dt_B` for accurate Brownian statistics. Both warnings are inert when `enable_brown = false`.
- **`compute_effective_body_brownian.m`**: Read-only numerical diagnostics injected into the existing computation path. (a) Before the mandatory symmetrisation of `M_eff_dimless`, the relative Frobenius asymmetry `||M-Mᵀ||_F / ||M||_F` is computed and compared against `1e-6`; exceeding that threshold emits `compute_effective_body_brownian:AsymmetryHigh` to prevent silent masking of `lin_sys` assembly bugs or catastrophic conditioning. (b) After dimensionalisation, `rcond(M_eff_dim)` is evaluated; values below `1e-10` emit `compute_effective_body_brownian:IllConditioned` to warn about near-singular effective mobility before the downstream Cholesky. Both hooks are diagnostic only — the numerical output (`L_eff`, `D_eff`, `M_eff_dim`) and function signature are unchanged.

### Changed

- **`.gitignore`**: Stopped ignoring `Brownian_Implementation_Plan.md` so the Brownian companion document ships with the repository (see **Added**).
- **`batch sweep/batch_repeat_brownian.m`**: Brownian seed mode `'default'` renamed to `'timestamp'` for self-documenting clarity (it is the wall-clock `MMddHHmmss` policy, never "the default"). The string `'default'` is retained as a silent backward-compatible alias: the script accepts it without warning and normalises it to `'timestamp'` internally before any downstream use. The default initial value of the `brown_rng_seed_mode` variable is updated from `'default'` to `'timestamp'`; existing user scripts or saved configs that still assign `'default'` continue to work identically.
- **`README.md`**, **`USER_GUIDE.md`**: Documentation of `brown_rng_seed_mode` updated to match the new canonical value; the alias is called out inline.
- **`Brownian_Implementation_Plan.md` §7**: File manifest entry for `validate_simulation_params.m` (lists the two new warnings) and `batch sweep/batch_repeat_brownian.m`.
- **`config_example.json`** (demo-only tweak): pili Young modulus `E = 1.768e7 → 3.279e8` (stiffer demo pili so the example reaches a visually interpretable regime faster); observation window `T_obs = 30.0 → 0.5` s (shorter default demo run). Purely illustrative — user-supplied JSON configs and inline cfg structs are unaffected; no numerical path or API change.

### Removed

- **`run_simulation.m` (post-5.5.6 cleanup)**: Experimental per-segment `results.seg_metrics` / "Phase 0 ODE profiler" instrumentation (briefly present in v5.5.6) has been removed. It only recorded `ode15s` bookkeeping statistics (`nsteps`, `nfevals`, `npds`, `first_dt`, `median_dt`) and did not affect physics; downstream code that never relied on `seg_metrics` is unchanged.

### Numerical equivalence

- All numerical paths are preserved bit-exactly. The `compute_effective_body_brownian.m` diagnostic hooks observe `M_eff_dimless` and `M_eff_dim` but do not mutate them; the existing symmetrisation, dimensionalisation, and Cholesky flow (Steps 1–5) is unchanged. `validate_simulation_params.m` adds only to the `warns{}` cell array; no error classification is altered. Re-running any pre-5.5.6 case with either Brownian off (`enable_brown=false`) or with `brown_rng_seed_mode = 'default'` / `'fixed'` yields identical trajectory output.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**, **`CHANGELOG.md`**: Version **5.5.6**, date **2026-04-21**; README feature-history row 76 added covering the Brownian polish and config-example tweaks; USER_GUIDE v5.5.6 paragraph added inline with the rest of the 5.5.x history. No other documentation sections required changes (file manifest and theory sections are unaffected by this release).

---

## [5.5.5] — 2026-04-21

### Performance

- **Direction 1 — `M_weighted` factorisation reuse (explicit LU via `decomposition`)**: The coupled mobility matrix `M_weighted` is now factorised **once per RHS call** via `decomposition(M_weighted, 'lu', 'CheckCondition', false)` and handed down to `ode_coupled_rhs.apply_Ainv_u` through `hydro.dM_weighted`. Rationale and effect:
    - `M_weighted` mixes **RSM** (body point-Stokeslets) and **RSS** (pili segment-integral Stokeslets) kernels, so the block matrix is **non-symmetric by construction** (empirically: $\lVert M-M^{\mathsf{T}}\rVert_F/\lVert M\rVert_F\approx0.56$). LU — not Cholesky — is the correct factorisation.
    - Before: every RHS call did **two independent** full-pivot LU factorisations of the same matrix (one inside `build_coupled_lin_sys` for `AQ_trans`, one inside `apply_Ainv_u` for the background-flow solve). Direction 1 reuses a single factorisation for both.
    - `'CheckCondition' = false` disables MATLAB's per-solve automatic `rcond` estimate (`O(N²)` saved per invocation); overall condition health is validated externally by the regression tests.
    - Even in the default **`bg = off`** workflow (no reuse), the explicit LU dispatch is marginally faster than generic backslash because it bypasses per-call symmetry / triangular heuristics.
    - **Measured on the reference benchmark (`benchmark_direction1.m`, DOFs = 1807)**:
        - `bg = off` (neutral path):  `ode_coupled_rhs` **767.5 ms → 731.6 ms (−4.7 %)** median.
        - `bg = uniform` (D1 hot path): `ode_coupled_rhs` **868.2 ms → 731.7 ms (−15.7 %)** median.

### Added

- **`hydro.dM_weighted`** (new field on the `hydro` struct returned by `build_coupled_lin_sys.m`): reusable LU decomposition object for `M_weighted \ rhs` solves. Consumed by `ode_coupled_rhs.apply_Ainv_u` (background-flow branch) to skip a second full LU factorisation per RHS call. The raw `hydro.M_weighted` is still provided for backward compatibility and diagnostics.
- **`_legacy/ode_coupled_rhs_v0.m`** (NEW, `_legacy/` folder): frozen pre-Direction-1 copy of `ode_coupled_rhs.m`, including a local `apply_Ainv_u_v0` helper that uses the raw `M_weighted \ u_t`. Used exclusively by the invariance-test and benchmark harnesses; never on the live simulation path.
- **`benchmark_direction1.m`** (NEW, root): direct wall-clock A/B benchmark comparing the legacy stack (`ode_coupled_rhs_v0` / `build_coupled_lin_sys_v0`) against the Direction-1 stack (`ode_coupled_rhs` / `build_coupled_lin_sys`) under both `bg = off` and `bg = uniform`. Reports median, min, max over `N_REP` interleaved iterations after `N_WARM` warm-up calls.

### Changed

- **`build_coupled_lin_sys.m`**: Replaced the single `AQ_trans = M_weighted \ Q_trans_corr` solve with an explicit two-step pattern — `dM_weighted = decomposition(M_weighted, 'lu', 'CheckCondition', false)` followed by `AQ_trans = dM_weighted \ Q_trans_corr`. The decomposition object is also exported as `hydro.dM_weighted`.
- **`ode_coupled_rhs.m::apply_Ainv_u`**: Now preferentially consumes `hydro.dM_weighted` (the reusable LU); falls back to the raw `hydro.M_weighted \` solve if the field is absent (backward compatibility).
- **`test_optimization_invariance.m`**: Extended with three new case groups:
    - **Case 1b / 2b (`AQ_trans`)**: bit-exactness of the optimised `dM \ Q` vs the legacy `M \ Q` at `Z0` and `Zp`.
    - **Case 3 (background-flow `apply_Ainv_u`)**: bit-exactness of `B_global · Ainv · u_bg` at `Z0` and `Zp` when `bg = uniform`.
    - **Case 4 (end-to-end `dZ`)**: `ode_coupled_rhs` vs `ode_coupled_rhs_v0` at `Z0` and `Zp`, for both `bg = off` and `bg = uniform`. `Zp` uses a condition-number-aware tolerance `DZ_TOL = 1e-9` (justified: `cond(lin_sys) ~ 1e9` at the perturbed state amplifies sub-ULP assembly noise inherited from the prior vectorisation round; the amplification is **not** introduced by Direction 1, which is bit-exact at all D1-owned checkpoints).

### Numerical equivalence

- Every D1-owned quantity (`M_weighted`, `AQ_trans`, `B_global · Ainv · u_bg`, `dZ` at `Z0`) passes at **`rel_err = 0`** (bit-exact). The `dZ` check at the perturbed state `Zp` passes at `4.2 × 10⁻¹¹` (well below the condition-number-amplified floor of `~1 × 10⁻⁹`); this residual is demonstrably carried over from the 5.5.4 `B_global` vectorisation round and is not introduced by Direction 1. No configuration, API, or physical behaviour is altered.

### Removed (post-optimisation cleanup, 2026-04-21)

- **One-shot developer harnesses** — purged from the live tree after the optimisation was validated and the results archived into this changelog:
    - `test_optimization_invariance.m` (equivalence regression; Cases 1/2/1b/2b/3/4).
    - `benchmark_direction1.m` (Direction-1 A/B wall-clock benchmark).
    - `profile_short_run.m` (short profiler-driver harness).
    - `compare_profiles.m` (profile A vs B comparator).
- **Retained**: `_legacy/ode_coupled_rhs_v0.m`, `_legacy/build_coupled_lin_sys_v0.m`, `_legacy/build_global_resistance_v0.m` — kept under `_legacy/` as frozen historical references; they are no longer referenced by any live or tool code and are **not** on the simulation run path.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**: Version **5.5.5**, date **2026-04-21**; feature-history row 75 added; file-structure section updated with the new `_legacy/ode_coupled_rhs_v0.m`, `benchmark_direction1.m`, and the `hydro.dM_weighted` reuse field.

---

## [5.5.4] — 2026-04-21

### Performance

- **`build_coupled_lin_sys.m` — hot-path vectorisation (7 sites)**: Replaced per-index `for` loops with strided indexing, triangular masks, and block assignments. End-to-end speed-up on the reference `ProfileShort` run: `build_coupled_lin_sys` **43.25 s → 35.88 s (−17.1 %)**, `ode_coupled_rhs` **43.22 s → 36.13 s (−16.4 %)**, `run_simulation` **95.97 s → 81.80 s (−14.8 %)**. Sites:
  1. **Body Q block** (identity + `-skew(r_i)` rows 1–3 × cols 1–6): strided direct writes `Q_global(idx_bx, 5) = rbz;` etc., no more per-body-point `skew_matrix` call.
  2. **Pili Q lower-triangular trig block**: `tril` mask + `cos/sin` vectorised element-wise fills (`Q12 = maskL .* (cT .* cP)` and the four siblings).
  3. **`row_order` interleave permutation**: `reshape([j_vec, (N_k+1)+j_vec, 2*(N_k+1)+j_vec]', [], 1)`.
  4. **`A_twist_total` assembly**: single `blkdiag(A_twist_cells{:})` call (O(M) vs the previous O(M²) iterative concatenation).
  5. **`B_global` body delta-force rows**: three strided writes `B_global(1, 1:3:3*N_body) = 1;` etc.
  6. **`B_global` angular-moment rows**: six strided writes using the precomputed quadrature weights (see "Added" below); replaces the per-point `get_trapezoidal_weight` lookup and per-iteration cross-product row.
  7. **`B_k` four sub-block assembly** (off-diagonal / diagonal / tip / tau): component-wise cross-product formulas with an upper-triangular mask for the off-diagonal block; MATLAB `cross(...,2)` for diagonal and tip; strided write for the tau block. Explicit write ordering (off-diag → diagonal → tip → tau) relies on IEEE-754 `x + 0 = x` so that mask-produced zeros in the diagonal column are correctly overwritten.
- **`build_global_resistance.m` — hot-path vectorisation (3 sites)**:
  1. **Stokeslet scatter** (`M_weighted` assembly): per-pilus block assignment via the shift-add pattern `[M1_k, 0] + [0, M2_k]` for each of the 9 `(a,b)` component matrices, replacing `N_seg_all × 18` cumulative element additions.
  2. **Rotlet scatter** (`M_rotlet` assembly): six strided writes `M_rotlet(row_x, 2:3:end) = Rxy;` etc., replacing the per-segment loop.
  3. **Torsional `Atilde_k` fill**: `sub2ind` + linear-index assignment `Atilde_k(lin_idx) = d3_k'` replaces the inner `for j = 1..N_k` loop.

### Added

- **`params.quad_weights`** (new field in `init_coupled_params.m`): Precomputed global trapezoidal quadrature weight vector of length `N_body + Σ(N_segments(k) + 1)`. Body points carry weight `1`; each pilus carries endpoint weight `ds/2` and interior-node weight `ds = ell(k)/N_segments(k)`. Values are static under the current (fixed-geometry) architecture and match the old `get_trapezoidal_weight` helper bit-exactly. Consumed by `build_coupled_lin_sys.m` angular-moment row assembly; eliminates ~`N_total × N_ODE_step` redundant weight lookups.
- **`test_optimization_invariance.m`** (NEW, root): Numerical equivalence test harness. Builds `params` / `Z0` with the same configuration as `profile_short_run.m`, then compares refactored `build_global_resistance` / `build_coupled_lin_sys` outputs against frozen reference copies `_legacy/build_global_resistance_v0.m` / `_legacy/build_coupled_lin_sys_v0.m` at two states (`Z0` and `Z0 + 1e-3*randn`). Reports relative Frobenius-norm error per output (`M_weighted`, `M_rotlet`, `A_twist_cells`, `lin_sys`, `hydro.*`, `geom.X_all`); **tolerance 1e-14**. Current pass: Case 1 (`Z0`) bit-exact (`rel_err = 0`); Case 2 (perturbed) `lin_sys` approximately `1.8e-20`, `hydro.B_global` approximately `3.4e-21`, all others `0`.
- **`compare_profiles.m`** (NEW, root): Side-by-side wall-time comparison of two `profile_info_*.mat` snapshots for hot-path functions (`build_coupled_lin_sys`, `build_global_resistance`, `ode_coupled_rhs`, `run_simulation`, `cross`, `blkdiag`).
- **`_legacy/build_coupled_lin_sys_v0.m`** and **`_legacy/build_global_resistance_v0.m`** (NEW, `_legacy/` folder): Frozen pre-optimisation copies with function names suffixed `_v0`. `build_coupled_lin_sys_v0` internally calls `build_global_resistance_v0`. Used exclusively by `test_optimization_invariance.m`; never on the live simulation path.

### Changed

- **`build_coupled_lin_sys.m`**: Removed the now-unused local helper `get_trapezoidal_weight` (fully superseded by `params.quad_weights`). Added structured code-block banners for the new vectorised sections with comments documenting the bit-exactness argument (IEEE-754 commutativity of two-operand addition; `0 + x = x`) and the overwrite ordering of the `B_k` sub-blocks.

### Numerical equivalence

- All changes are bit-exact (rel. Frobenius `= 0`) at the canonical initial state; under arbitrary perturbation the largest observed deviation is `1.8 × 10⁻²⁰`, well below the `1 × 10⁻¹⁴` equivalence threshold and several orders of magnitude below machine epsilon. No configuration, API, or physical behaviour is altered — all existing JSON configs, sweep drivers, and ODE tolerances continue to reproduce identical trajectories.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**: Version **5.5.4**, date **2026-04-21**; feature-history row 74 added; file-structure section updated with the new `_legacy/`, `test_optimization_invariance.m`, `compare_profiles.m`, and the `params.quad_weights` precomputation note.

---

## [5.5.3] — 2026-04-20

### Changed

- **Default ODE timeout `tlim`**: Project-wide default is now **3600 s** (was **360 s**) — `load_sim_config.m` JSON fallback, `main_coupled_simulation.m` inline config, `config_example.json`, `docs/index.html` JSON snippet, all `batch sweep/*` scripts that set `cfg_base.tlim`, and **`README.md`** §2 parameter table.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**: Version **5.5.3**, date **2026-04-20**; `tlim` default and feature-history row 73.

---

## [5.5.2] — 2026-04-20

### Changed

- **`capsule_body_building.m` — cylinder circumferential resolution**: After the area-based estimate `TN = ceil(CPN / Turns)`, the code now snaps **`TN`** to a multiple of 4 with **`TN = max(round(TN / 4) * 4, 4)`** before constructing the azimuthal grid (`Alpha`, `TempY`, `TempZ`, `CLoc`). Total cylinder points become **`TN * Turns`** with the snapped **`TN`**, so actual **`N_body`** may differ further from the requested input.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`docs/index.html`**: Version **5.5.2**, date **2026-04-20**; capsule cylinder grid note updated.

---

## [5.5.1] — 2026-04-18

### Added

- **`sphere_particle_sample_cached.m`** (NEW): Disk cache for `ParticleSampleSphere` results when `100 < SPN < 1000`. On first call for a given `SPN` the converged point set is saved to `cache/sphere_particle/spherical_points_NXXXX.mat`; subsequent calls load from disk and bypass the iterative Riesz-energy optimisation entirely. Each cache file stores `Vtool` (N×3, toolbox frame) plus a `meta` struct containing the toolbox-default tolerances (`Etol`, `Dtol_deg`, `Nitr`), the initial-position strategy (`Vo_method`), and a format version. Any mismatch between stored meta and current constants triggers recomputation and overwrites the stale file; the old file is also deleted on any structural or meta validation failure so corrupt entries do not persist. All I/O errors fall back silently to direct computation, leaving the simulation unaffected.

### Changed

- **`capsule_sphere_cap_points.m` — `ParticleSampleSphere` call**: Removed custom `Etol` / `Dtol` / `Nitr` overrides; the toolbox defaults (`Etol=1e-5`, `Dtol=1e-4 deg`, `Nitr=1e4`) now apply. Initial particle positions (`'Vo'`) are supplied as `SpiralSampleSphere(SPN, false)` — fully deterministic, no `RandSampleSphere` / RNG dependency. For `100 < SPN < 1000` the call is routed through `sphere_particle_sample_cached`.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**: Version **5.5.1**, date **2026-04-18**; sphere-cap discretisation description updated.

---

## [5.5.0] — 2026-04-18

### Added

- **Capsule sphere-cap discretisation (S2 Sampling Toolbox)** — `capsule_body_building.m` now builds hemispherical-cap Stokeslet points via Anton Semechko’s add-on (`ParticleSampleSphere`, `SpiralSampleSphere`, optionally Ico / quad meshes from `IcosahedronMesh` / `QuadCubeMesh` + `SubdivideSphericalMesh`). New helpers: **`capsule_sphere_mesh_cache.m`** (persistent cache of subdivided meshes and **measured** vertex counts for `k = 0…4`), **`snap_sphere_cap_count.m`** (snap area-weighted cap count to nearest native total if within ±5), **`capsule_sphere_cap_points.m`** (mesh / particle / spiral branches; documented `Etol`, `Dtol`, `Nitr` for `ParticleSampleSphere`). Native vertex totals match the **installed** toolbox version (no closed-form assumption).

### Changed

- **`capsule_body_building.m`** — Duplicates cap points on the **strict** equator **`x == 0`** onto **both** hemispherical caps after translation. Final `NBody` is `numel(rxB)` (may differ from requested `NBody` when snap or duplication applies). Cylindrical patch rule unchanged.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**: Version **5.5.0**, date **2026-04-18**; S2 toolbox prerequisite and cap-algorithm summary.

---

## [5.4.4] — 2026-04-15

### Changed

- **`main_coupled_simulation.m`**: Inline example case updated to **`OnePili-Active-HighFre-Test`** with a single pilus at site 1 and a higher-frequency `traveling_wave` active moment, giving a simpler active-motion demonstration template.

### Fixed

- **`run_simulation.m` — MP4 animation export (`bacteria_3d_motion.mp4`)**: Hardened Windows HiDPI / display-scaling handling. After the figure is resized and probed, the encoded frame size is now clamped to a codec-safe, even-valued MPEG-4 box (`<= 3840 x 2160`). If `getframe` still returns oversized physical frames, `frame.cdata` is resized before `writeVideo`, preventing failures where logical figure shrinking alone did not reduce the captured physical pixel dimensions enough.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**: Version **5.4.4**, date **2026-04-15**; animation export notes updated to match the current HiDPI-safe implementation.

---

## [5.4.3] — 2026-04-15

### Changed

- **Header comment cleanup (all 47 `.m` files)**: Removed all `%  Updated:` history lines from file headers. Header comments now end with `Author:` / `Created:` only, keeping the source tree concise and relying on `CHANGELOG.md` for change history.

### Fixed

- **`build_coupled_lin_sys.m`**: Syntax comment now lists both 2-output and 3-output calling forms (`[lin_sys, geom]` and `[lin_sys, geom, hydro]`), matching the actual function signature.
- **`brown_rng_seed_from_clock.m`**: Added standard `Description:` / `Syntax:` / `Inputs:` / `Outputs:` header sections (previously had only a brief one-liner).
- **`post_processing/postprocess_cases.m`**: Standardised header comment to project style (`%POSTPROCESS_CASES`, two-space indent, `Syntax:` / `Inputs:` / `Outputs:` sections, `Author:` / `Created:`; was using inconsistent `% POSTPROCESS_CASES` + `Date:` format).

---

## [5.4.2] — 2026-04-15

### Added

- **`traveling_wave` active moment**: optional **`phase_rad`** [rad] (default 0) — constant phase added to \(\Phi = k_w s - \omega t\) in `coupled_active_moment.m`; set in `init_coupled_params` per pilus. JSON / `validate_simulation_params` accept a finite scalar.

---

## [5.4.1] — 2026-04-15

### Removed

- **`breaststroke` active moment mode** from `coupled_active_moment.m`, `init_coupled_params.m`, and `validate_simulation_params.m`. Per-pili active moments remain: `off`, `traveling_wave`, and `custom`. Second-pilus example in `config_example.json` now uses `traveling_wave`.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`docs/index.html`**, **`run_simulation.m`**: Dropped references to `breaststroke`.

---

## [5.4.0] — 2026-04-14

### Added

- **`coupled_active_moment.m` — `breaststroke` mode** (NEW): Two-harmonic asymmetric beat pattern modelling Chlamydomonas-like breaststroke flagellar motion. Moment density: $m(s,t) = A \cdot E(s) \cdot [\cos\Phi + C_2\cos(2\Phi+\varphi_2) + C_0]$ where $\Phi = k\,s - \omega\,t$ and $E(s) = (1-\beta)+\beta\,s$ is a linear spatial envelope. The linear envelope admits **exact antiderivatives** (same precision as `traveling_wave`). Three mechanisms absent from the single-harmonic `traveling_wave`:
  - **Temporal asymmetry** (`C2_ratio`, `C2_phase`): second harmonic creates distinct power/recovery stroke durations.
  - **Spatial envelope** (`envelope_slope`): amplitude ramps from base to tip.
  - **Mean curvature offset** (`mean_curv`): beat oscillates around a curved rest shape.
  Setting `envelope_slope = 0`, `C2_ratio = 0`, `mean_curv = 0` reduces exactly to `traveling_wave`.
- **`init_coupled_params.m`**: `case 'breaststroke'` — same `amplitude_dimless` / `omega_dimless` conversion as `traveling_wave` plus defaults for `envelope_slope` (0.5), `C2_ratio` (0.25), `C2_phase` (0), `mean_curv` (0).
- **`validate_simulation_params.m`**: Validation rules for breaststroke — `envelope_slope` in [0, 1], `C2_ratio >= 0`; warnings for `C2_ratio > 1`, `|mean_curv| > 10`.
- **`main_coupled_simulation.m`**: Commented inline breaststroke configuration example.
- **`config_example.json`**: Second pili entry with full breaststroke JSON template.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**: Version **5.4.0**, date **2026-04-14**.

---

## [5.3.2] — 2026-04-14

### Added

- **`compute_node_directors.m`** (NEW): Standalone utility that interpolates per-segment director frames to per-node directors. Internal nodes use neighbour averaging; boundary nodes use linear extrapolation with `min(2, N_k)` safety for single-segment pili. All averaged frames are projected onto the nearest proper rotation via **SVD polar decomposition** with an explicit `det < 0` correction to prevent reflections. Shared by `ode_coupled_rhs.m` (active moment evaluation) and `build_coupled_lin_sys.m` (geometry reconstruction).

### Fixed

- **`build_coupled_lin_sys.m` — node director interpolation**: The inline implementation lacked the `det(R) < 0` SVD fix, so polar decomposition could silently return an improper rotation (reflection) instead of the nearest proper rotation. Also had a potential index-out-of-bounds for single-segment pili (`N_k = 1`) in the boundary extrapolation. Both issues fixed by replacing inline code with a call to the new shared `compute_node_directors.m`.

### Changed

- **`ode_coupled_rhs.m`**: (1) Elastic-moment assignment to `R_global` and active-moment subtraction **vectorised** (replaced per-element `for` loops with `reshape` + array indexing). (2) Local `compute_node_directors` helper removed; calls shared `compute_node_directors.m` instead — eliminates code duplication with `build_coupled_lin_sys.m`.
- **`validate_simulation_params.m`**: Added warnings for `traveling_wave` active moment parameters: negative `wavenumber` (wave propagates tip → root) and non-zero `d3_amp` (drives torsion rather than bending, unusual for standard pili).

### Documentation

- **`coupled_active_moment.m`**: Clarified `d3_amp` physical meaning — non-zero value drives twist (torsional component), not bending; typically zero for standard flagellar/pili beating.
- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**: Version **5.3.2**, date **2026-04-14**.

---

## [5.3.1] — 2026-04-14

### Fixed

- **`run_simulation.m` — MP4 animation (`bacteria_3d_motion.mp4`)**: Two failure modes fixed:
  1. **Large Brownian trajectories** push computed figure width beyond Windows MPEG-4 / H.264 encoder limits. Logical figure size is now **proportionally scaled** to fit within a **3840 x 2160** codec-safe box.
  2. **Hi-DPI / display scaling** (e.g. 125%): `getframe` returns **physical pixels** (e.g. 2250 x 1366) larger than the logical `Position` (1800 x 1200). The codec limit was checked against logical pixels only, which was ineffective. Fix: after creating the figure, a **probe `getframe`** measures the actual physical size; if it exceeds codec limits, the figure `Position` is **shrunk** and re-probed. `VideoWriter` is opened **after** the physical frame size is known. In the render loop, `getframe(fig)` is called **without a rect** (avoids logical/physical coordinate mismatch). A fallback **`normalize_anim_frame_cdata`** (`imresize` or `interp2`) catches any residual per-frame size drift.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**: Version **5.3.1**, date **2026-04-14**.

---

## [5.3.0] — 2026-04-13

### Added

- **`batch sweep/batch_sweep_N_body_eps_RSM.m`**: Nested sweep of **`N_body`** and **`eps_body = k · h_mean`** (RSM regularisation) on a sphere with no pili; compares steady **`V_sink`** to analytic Stokes settling speed; writes `batch_sweep_N_body_eps_RSM_grid.csv`, `batch_sweep_N_body_eps_RSM_summary.csv`, and `.mat` results.

### Changed

- **Default body RSM regularisation** (`init_coupled_params.m`): **`eps_body`** default is now **`h_mean / 2.6667`** (was **`h_mean / 3`**), i.e. \(\varepsilon/\Delta s_{\mathrm{avg}} \approx 0.375\), from the above sweep vs Stokes reference.
- **`build_global_resistance.m`**: Header comment aligned with the new default ratio.

### Documentation

- **`README.md`** (§8.5, JSON/`eps_body` notes, feature history), **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**: Version **5.3.0**, date **2026-04-13**.

---

## [5.2.1] — 2026-04-10

### Changed

- **Default `N_body`**: project-wide default for body surface Stokeslet points is now **200** (was **250**) — `load_sim_config.m`, `init_coupled_params.m`, `main_coupled_simulation.m`, `config_example.json`, `batch sweep/*.m`, `docs/index.html`, `README.md`.
- **Default `tlim`**: ODE segment timeout limit default is now **360** s (was **120** s) in the same locations.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**, **`run_simulation.m`**: Version **5.2.1** (log banner).

---

## [5.2.0] — 2026-04-10

### Added

- **Brownian RNG seed policy** (`run_simulation.m`, `brown_rng_seed_from_clock.m`): When `enable_brown` is true, MATLAB's global generator is initialised with **`rng(seed, 'twister')`**. Default seed is **local wall-clock time** encoded as integer **`MMddHHmmss`** (via `brown_rng_seed_from_clock`). Optional **`cfg.brown_rng_seed`** (positive scalar) or JSON **`brownian.rng_seed`** forces a user-defined reproducible seed.
- **Seed persistence**: The realised seed is written to **`output/<case>/config.json`** (`brownian.rng_seed`), **`brownian_rng_seed.txt`**, **`simulation.log`** (`[PARAMS] Brownian RNG`), **`results.brown_rng_seed`**, and **`params.brown_rng_seed`** in `simulation_results.mat`.
- **`load_sim_config.m` / `save_sim_config.m`**: Parse and emit optional **`brownian.rng_seed`** ↔ flat **`cfg.brown_rng_seed`**.

### Changed

- **`batch sweep/batch_repeat_brownian.m`**: Brownian seeds are applied only inside **`run_simulation`** (not via a separate `rng()` call). **`brown_rng_seed_mode`**: `'default'` omits **`cfg.brown_rng_seed`** so each repeat uses the same **MMddHHmmss** time default as a single run (`brown_rng_seed_from_clock`); `'fixed'` sets **`cfg.brown_rng_seed = rng_seed_base + k`** for reproducible batches.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**: Version **5.2.0**, date **2026-04-10**; Brownian seed and output files documented.

---

## [5.1.1] — 2026-04-08

### Changed

- **Default discretisation**: project-wide defaults for **`N_body`** (body surface Stokeslet points) and **`seg_per_um`** (RSS segments per micron along each pilus) are now **250** and **10** respectively (`load_sim_config.m`, `init_coupled_params.m`, `main_coupled_simulation.m`, JSON templates, `batch sweep/*.m`, `docs/index.html`, `README.md`).
- **`run_simulation.m`**: After `exportgraphics` saves `summary_plot.png`, the summary figure is **`close`d** so batch runs do not accumulate visible figure windows.
- **`compute_plot_flowfield.m`**: After `exportgraphics` saves `flowfield_panels.png`, the flow-field figure is **`close`d** for the same reason.
- **`run_simulation.m`**: Log banner reports **Version 5.1.1**.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**: Version **5.1.1**, date **2026-04-08**; §17 notes batch-friendly figure cleanup.
- **JSON / config docs**: **`config_example.json`** and **`configs/OnePili-P1-L3.json`** updated to include **`background_flow`** (v5.0) and per-pilus **`intrinsic_curvature`** / **`active_moment`** defaults; **`README.md`** §2.2 documents nested JSON keys vs flat `cfg`; **`USER_GUIDE.md`** §6 adds JSON fragments; **`docs/index.html`** quick-start JSON aligned with `save_sim_config`, **`traveling_wave`** field names corrected to **`d1_amp`** / **`d2_amp`** / **`d3_amp`**.

---

## [5.1.0] — 2026-04-07

### Added

- **`post_processing/postprocess_cases.m`**: Batch post-processing of **multiple** simulation result folders under a chosen parent directory (non-recursive: only first-level subfolders containing `timeseries.csv`). Reads each case's `timeseries.csv` and `config.json`, computes summary metrics (sinking velocity statistics, maximum horizontal displacement, 3-D trajectory tortuosity, translational diffusion coefficient via MSD fit $6D\tau + V^2\tau^2$), and exports a **Summary** sheet to Excel (`.xlsx`). On Windows, optional COM formatting (bold headers, borders, row shading) and optional user prompt before opening the file. English UI dialogs; progress via `waitbar` and console output.

### Changed

- **`run_simulation.m`**: Log banner reports **Version 5.1.0**.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**, **`docs/index.html`**: Version **5.1.0**, date **2026-04-07**; post-processing workflow documented.

---

## [5.0.0] — 2026-04-03

### Added

- **`coupled_background_flow.m`** (NEW): Evaluates a steady (time-invariant) background flow field at all hydrodynamic evaluation points (body surface + pili nodes). Supports three modes: `'off'` (default, zero flow), `'uniform'` (spatially constant velocity), and `'shear_z'` (piecewise-linear shear). Output vector matches the DOF layout of `build_coupled_lin_sys` (translation block + twist block). All values pre-converted to dimensionless during initialisation.
- **`coupled_active_moment.m`** (NEW): Computes per-pilus time-varying active internal moments at the current dimensionless time. Supports `'off'` (default), `'traveling_wave'` (sinusoidal traveling wave with configurable amplitude, frequency, wavenumber, and per-component weights), and `'custom'` (user-provided function handle). Algorithm follows `3D-Flexible-Filaments/integrated_internal_moments.m` (antiderivative → segment integrals → reverse cumsum → director-basis projection).
- **Intrinsic curvature configuration**: Per-pilus static intrinsic curvature (`params.kb0{k}`) can now be set to `'frozen'` (default, keep initial geometry), `'zero'`, `'constant'` (user-specified `[k1, k2, k3]`), or `'helix'` (helical curvature from `kappa_d2`/`kappa_d3`). Configured via `pili_config(k).intrinsic_curvature` in the inline config or per-pilus JSON fields.
- **Active moment configuration**: Per-pilus time-dependent active moments configured via `pili_config(k).active_moment` (struct with `.type`). Pre-converts physical parameters (amplitude [N·m], frequency [Hz]) to dimensionless at initialisation.
- **Background flow configuration**: Global steady background flow configured via `cfg.background_flow` (struct with `.type`, `.U_max`, `.direction`, `.shear_axis`, `.shear_height_m`). Pre-converts physical velocity to dimensionless at initialisation.

### Changed

- **`build_coupled_lin_sys.m`**: Added optional 3rd output `hydro` struct (contains `B_global`, `M_weighted`, `A_twist_total`, `M_rotlet`, `n_trans`) for background flow computation. Backward-compatible: existing 2-output calls are unaffected (`nargout >= 3` guard).
- **`ode_coupled_rhs.m`**: Integrated active moments (Step 8b: subtracted from `R_global` per-pilus before clamped BC) and background flow (Step 8c: `B*Ainv*u_bg` subtracted from `R_global`). Added local helpers `apply_Ainv_u` and `compute_node_directors`.
- **`init_coupled_params.m`**: Added three Name-Value pairs (`'intrinsic_curvature'`, `'active_moment'`, `'background_flow'`). Intrinsic curvature overrides `params.kb0{k}` per mode; active moment configs pre-convert physical units to dimensionless; background flow pre-computes `U_max_dimless` and `shear_height_dimless`. Summary printout extended with new feature status per pilus.
- **`parse_pili_config.m`**: Added optional 2nd output `pili_ext` (struct array with `.intrinsic_curvature` and `.active_moment` per pilus). Backward-compatible: 1-output calls still work; `Pili_Matrix` [8×M] unchanged.
- **`load_sim_config.m`**: Parses `background_flow` JSON section and collects `pili_ext` from `parse_pili_config`.
- **`save_sim_config.m`**: Writes `background_flow` section and per-pilus `intrinsic_curvature`/`active_moment` fields to JSON. Function handles (`'custom'` mode) are stripped before serialisation.
- **`validate_simulation_params.m`**: Added validation rules for background flow type/parameters, per-pilus intrinsic curvature type/kappa, per-pilus active moment type/amplitude/frequency; warnings for background flow + Brownian combination and active moments without pili.
- **`main_coupled_simulation.m`**: Added inline config examples for background flow, per-pilus intrinsic curvature, and per-pilus active moment (all default to `'off'`/`'frozen'`). Uses `[cfg.Pili_Matrix, cfg.pili_ext] = parse_pili_config(pili_config)`.
- **`run_simulation.m`**: Extracts new feature fields from `cfg`, passes `'intrinsic_curvature'`, `'active_moment'`, `'background_flow'` to `init_coupled_params`; adds new fields to `val_cfg` for validation. Log banner reports **Version 5.0.0**.

### Compatibility

- All three new features default to `'off'`/`'frozen'`: existing configurations and sweep scripts produce **identical results** without modification.
- Background flow is deterministic and does not affect the Brownian diffusion tensor (`compute_effective_body_brownian` probes `lin_sys` which excludes flow contributions).
- Active moments are applied before clamped BC rows are zeroed, avoiding any boundary condition conflict.
- `adaptive_basis_check.m` preserves `params.am_cfg{k}` and `params.bg_flow` across basis rotations (static config, not state-dependent).

---

## [4.1.0] — 2026-04-03

### Added

- **`batch sweep/batch_repeat_brownian.m`**: Batch independent repeats of the **same** Brownian-on configuration. User sets repeat count and `case_name_base`; each run writes to `output/<case_name_base>/001/`, `002`, … (zero-padded). Requires `enable_brown = true`; forces `preview_initial = false` and `auto_confirm = true`; sets `rng(rng_seed_base + k)` before each repeat for reproducible independent noise streams. Optional JSON config via `config_json`, else inline `cfg_base` in the script. Writes `batch_repeat_summary_*.txt` under `output/<case_name_base>/`.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**: Version **4.1.0**, date **2026-04-03**; file tree and feature bullets updated for the new batch script.
- **`docs/index.html`**: Version badge/footer aligned with release **4.1.0**.

### Changed

- **`run_simulation.m`**: Log banner reports **Version 4.1.0** (aligned with `README.md`).

---

## [4.0.0] — 2026-03-31

### Added

- **`build_coupled_lin_sys.m`** (NEW): Extracted Steps 1-7 from `ode_coupled_rhs.m` into a reusable subroutine that assembles the coupled body-pili linear system with clamped BC. Shared by the ODE RHS and the Brownian mobility computation.
- **`compute_effective_body_brownian.m`** (NEW): Computes the full-system effective 6×6 body mobility by probing the coupled linear system, dimensionalises it, and returns the Cholesky factor for Euler-Maruyama Brownian dynamics. Implicitly captures pili hydrodynamic drag via Schur complement.
- **`fixman_midpoint`** configuration option: Optional Fixman midpoint correction for configuration-dependent diffusion drift (2× mobility evaluation per Brownian step).

### Changed

- **`ode_coupled_rhs.m`**: Refactored to call `build_coupled_lin_sys` for Steps 1-7, retaining only elastic moment computation and state-rate packing.
- **`run_simulation.m`**: Brownian kick block replaced: uses per-step dynamically computed `L_eff` (6×6 Cholesky) from full-system effective friction; generates 6D correlated noise (translation-rotation coupling included); supports optional Fixman midpoint correction.
- **`init_coupled_params.m`**: Removed one-time body-only `compute_body_brownian_params` call at initialisation; Brownian diffusion is now computed dynamically per step; added `fixman_midpoint` parameter.
- **`load_sim_config.m`**: Parses new `fixman_midpoint` field from the `brownian` JSON block.
- **`config_example.json`**: Added `fixman_midpoint: false` to the `brownian` block.
- **`validate_simulation_params.m`**: Added validation for `fixman_midpoint` parameter.
- **`Brownian_Implementation_Plan.md`**: Added upgrade notice documenting the transition from body-only to full-system effective friction.
- **`README.md`**: Updated section 16, file listings, architecture diagram, and feature descriptions.

### Removed

- **`compute_body_brownian_params.m`**: Deleted. Old body-only RSM friction matrix implementation fully superseded by `compute_effective_body_brownian.m`.
- **`test_effective_brownian.m`**: Deleted. Verification tests passed (7/7); script no longer needed.

### Fixed

- **`main_coupled_simulation.m`**: Added missing `cfg.fixman_midpoint` to inline configuration block (runtime crash when using inline config path).
- **`configs/OnePili-P1-L3.json`**: Added missing `fixman_midpoint` field to the `brownian` block.

---

## [3.3.5] — 2026-03-30

### Changed

- **`batch sweep/batch_sweep_major_axis.m`**: semi-major axis sweep grid **`0.4:0.2:1.6` µm** (7 cases); **`minor_axis`** fixed at **0.4 µm** (first case is a sphere, `major_axis == minor_axis`).

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**: version **3.3.5**, date **2026-03-30**; feature table row 56.
- **`run_simulation.m`**: log banner **Version 3.3.5**.

---

## [3.3.4] — 2026-03-30

### Removed

- **`batch sweep/batch_sweep_pili_site.m`**: removed; attachment-site sweeps use **`batch_sweep_pili_length.m`** with user-defined **`site_list`**.

### Changed

- **`batch sweep/batch_sweep_major_axis.m`**: **body-only** major-axis sweep — **`Pili_Matrix` empty** (no pili).

### Documentation

- **`README.md`**: file counts (38 `.m` files; 6 sweep drivers + `batch_run_from_json`), feature table row 55; `batch_sweep_major_axis` described as no-pili body.

---

## [3.3.3] — 2026-03-30

### Changed

- **`batch sweep/batch_sweep_alpha.m`**: **Body-only** initial tilt sweep — **`Pili_Matrix` empty** (no pili); fixed ellipsoidal capsule geometry unchanged.

### Documentation

- **`README.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**: version **3.3.3**, date **2026-03-30**; feature table row 54; `batch_sweep_alpha` description clarified.
- **`run_simulation.m`**: log banner **Version 3.3.3**.

---

## [3.3.2] — 2026-03-30

### Added

- **`batch sweep/batch_sweep_N_body.m`**: sweep `N_body` for a **spherical** body (radius **0.4 µm**, `major_axis == minor_axis`) with **no pili** (`Pili_Matrix` empty).
- **`batch sweep/batch_sweep_seg_per_um.m`**: sweep `seg_per_um` (default grid 10:5:35).

### Removed

- **`batch sweep/batch_sweep_pili_count.m`**: removed (multi-pili / site control covered by `batch_sweep_pili_length.m` + JSON batching).

### Changed

- **All `batch_sweep_*.m` drivers** (excluding `batch_run_from_json.m`): `case_name` tokens joined with **`-`**; **`enable_brown = false`**; **`T_obs = 200`**; **`generate_animation = true`**, **`animation_skip = 5`**; default pili Young's modulus **`E = 3.279e8` Pa**.
- **`batch_sweep_pili_length.m`**: user-defined **`site_list`** for arbitrary per-pilus attachment sites (equal length per run).
- **`batch_sweep_pili_E.m`**: E sweep **`1e5`–`1e9` Pa**, tenfold per step (`10.^(5:9)`).
- **`validate_simulation_params.m`**: allow **`major_axis == minor_axis`** (sphere); error only when **`major_axis < minor_axis`**.
- **`README.md`**: batch sweep inventory (7 sweep scripts + JSON runner), file counts, `N_body` script description.

---

## [3.3.1] — 2026-03-30

### Added

- **`CHANGELOG.md`**, **`USER_GUIDE.md`**, **`requirements.txt`**: release notes, concise user guide, and MATLAB-only environment note (no pip dependencies).

### Changed

- **`README.md`**: Version/date bump; §18 flow-field description; §19 file tree; new §25 companion docs; feature table rows 51–52.
- **Flow-field 2D quiver** (`compute_plot_flowfield.m`): Removed the white underlay (double-draw halo) on XZ and YZ panels; arrows are drawn once in dark gray for clearer interpretation on the jet colormap. Arrow length, line width, and head size are scaled **+50%** relative to the previous dark layer for improved visibility.
- **`run_simulation.m`**: Log banner reports **Version 3.3.1** (aligned with `README.md`).

---

## [3.3.0] — 2026-03-18

### Added

- Parameter sweep scripts consolidated under `batch sweep/` (six sweep drivers + shared patterns).
- `batch_run_from_json.m`: batch runner over a folder of JSON configs with summary log.

### Fixed

- Minor batch-script hygiene (`batch_run_from_json.m` dead code / redundant assignments).

---

## Earlier releases

Feature and fix history prior to 3.3.0 is summarized in `README.md` (sections *22. 验证总结* — bugfix and feature tables).
