# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))

---

## [Unreleased]

---

## [8.7.1] — 2026-10-30

Flow-field figure patch.  **The ODE trajectory and the saved
`results.flow_state` are unchanged.**

### Changed

- **Body interior.**  The custom-body mask is a solid volume.  Velocity is
  not drawn inside the cell.  The white outline is the surface–plane
  section, resampled with a periodic spline, instead of the stair-stepped
  contour of the velocity grid.
- **Axis window.**  Each plane covers the configured range and is extended,
  with a margin, so the body and every flagellum node lie inside the axes.
- **Figure layout.**  Panel letters sit outside the top-left of each axes
  box.  The only colour bar is vertical, beside panel d.  The figure title
  is the time at 0.0001 s resolution.  Lighting on the 3-D panels comes
  from `[1, -1, 0]` toward the origin.  The footer note is removed.

### Compatibility

- No new dependency.  Results from 8.7.0 can be redrawn with the new
  figure.  Results from 8.6.6 or earlier still cannot reconstruct a flow
  field.

### Documentation

- `project_version.m` bumped to 8.7.1 (2026-10-30); README (EN/ZH), USER_GUIDE
  (EN/ZH), theory (EN/ZH), `requirements.txt`, the project website and
  presentation version strings, and the headers of `run_simulation*.m` and
  `main_coupled_simulation.m` updated.

---

## [8.7.0] — 2026-10-29

Flow-field post-processing release.  **The ODE trajectory is unchanged.**
Every new run stores the state needed to reconstruct the flow at each saved
output.  Results written before 8.7.0 do not contain `results.flow_state`
and cannot be used for this post-processing.

### Added

- **Saved flow-field state.**  `results.flow_state` (version 1) records, at
  every stored output, the full dimensionless state `Z`, the pilus basis
  `rot{k}`, and the single Brownian generalised force that produced that
  output.  The initial sample has zero thermal force.  The last sample keeps
  the force of the Brownian interval that produced it; the saved
  `params.brown_force_dimless` is still cleared after integration.
- **`post_processing/render_flowfield_from_results.m`.**  Independent entry,
  with the same file-dialog style as `render_animation_from_results`.
  Modes: `final` (default), `times` (nearest saved output in seconds;
  duplicate matches are computed once), `animation` (`fps` and `skip`, last
  frame always kept), and `all`.  Static multi-time output is PNG only.
  The animation is MPEG-4 `.mp4` on Windows and macOS, and Motion JPEG AVI
  `.avi` on Linux.  New files do not overwrite `flowfield_panels.png`.
- **Four-panel flow figure.**  3-D and 2-D XZ/YZ views on a fixed canvas,
  panels (a)–(d), one colour bar, and a note that the colour is the full
  speed while the arrows are normalised in-plane directions.  The window
  follows the centre of mass.  The final panel uses that frame's 99th
  percentile; selected times use only the requested frames; an animation
  uses every saved output.

### Changed

- **One hydrodynamic assembly per flow frame.**  `ode_coupled_rhs` can
  return the mobility matrices already built for `dZ`.  Callers that request
  only `dZ` do not copy them, and the discrete trajectory is unchanged.
  `flowfield_from_saved_state` solves the force and torque coefficients from
  those matrices with the existing backslash path.
- **One grid evaluation per required frame.**  Final and selected-time
  figures are drawn from the grids just computed.  An animation evaluates
  every saved output once for the global colour scale and keeps grids only
  for frames that will be shown, up to 512 MiB.  Frames beyond that limit
  are redrawn from the saved force coefficients; the velocity kernel is
  unchanged.
- **End-of-run flow figure.**  `run_simulation` and
  `run_simulation_HPC_Linux` still draw only the final state when
  `compute_flowfield` is true, and that figure now uses the thermal force
  recorded for the last saved output.

### Compatibility

- The ODE, Brownian sampling, and saved trajectories of quantities other
  than `results.flow_state` are unchanged.
- A `simulation_results.mat` from 8.6.6 or earlier can still re-render the
  motion animation.  It cannot reconstruct a flow field.
- Linux HPC batch sampling of `compute_flowfield` is unchanged.  Every new
  repeat still saves `flow_state`.
- No new dependency.

### Documentation

- `project_version.m` bumped to 8.7.0 (2026-10-29); README (EN/ZH), USER_GUIDE
  (EN/ZH), theory (EN/ZH), `requirements.txt`, the project website and
  presentation version strings, and the headers of `run_simulation*.m` and
  `main_coupled_simulation.m` updated.

---

## [8.6.6] — 2026-10-21

Visualisation fix.  **The hydrodynamic solver is unchanged.**  Parametric
results are unchanged.  Custom-body surface reconstruction can change, and
with it the enclosed volume, the volume-based length scale, and the sinking
force when the cell is not neutrally buoyant.  Motion-animation videos are
not pixel-comparable with 8.6.3–8.6.5.

### Fixed

- **Custom-body alpha shape.**  `try_alpha_shape` used to accept the first
  closed piece it found while scanning from a small alpha upward.  On the
  *E. coli* capsule that piece was a 4-vertex tetrahedron, so the summary
  figure and the animation drew a single facet and the enclosed volume was
  \(4.1\times10^{-24}\,\mathrm{m}^3\).  The scan now starts at the one-region
  critical alpha and keeps the tightest closed genus-0 surface that still
  uses at least 90 % of the input points.  The same capsule reconstructs as
  326 vertices and 648 faces, with volume \(1.31\times10^{-18}\,\mathrm{m}^3\)
  (about 97 % of the analytic capsule).  The regularised-Stokeslet point
  cloud is still the input cloud; only the display mesh and the
  divergence-theorem volume change.

### Changed

- **Co-moving windows (c, d).**  Tick spacing is never finer than 0.5 μm.
  The step is a 1-2-5 value at or above that floor, refreshed as the window
  follows the centre of mass, so the grid still slides past the bacterium.
- **E. coli example time limit.**  `configs/config_ecoli_single_flagellum.json` sets `simulation.tlim` to 72000 s.  The previous 3600 s limit can expire before `T_obs = 1 s` finishes.
- **Roll marker.**  The two dark longitudinal stripes are removed.  The body
  is uniform light green again (`FaceAlpha` 0.35).  One short dark-green
  patch lies on one side of the waist: the central 18 % of the body length
  and ±32° of azimuth, 2 % outside the local radius.  It rotates with the
  body.  The legend key is the light-green body with that small patch.

### Compatibility

- Parametric geometry, the ODE, and results of parametric cases are
  unchanged.
- A custom body whose previous alpha shape was already a proper closed
  surface can still receive a different mesh, volume, and volume-based
  `L_0`.  When `rho_bac` differs from `rho_fluid`, `Eg` changes with the
  volume.  The *E. coli* example is neutrally buoyant, so its sinking force
  stays zero; a re-run is not bit-identical because the length scale used
  by that case was taken from the old tetrahedron volume.
- Saved `simulation_results.mat` files keep the mesh they were written
  with.  Re-rendering an old file does not rebuild the body surface.
- No new dependency.

### Documentation

- `project_version.m` bumped to 8.6.6 (2026-10-21); README (EN/ZH), USER_GUIDE
  (EN/ZH), theory (EN/ZH), `requirements.txt`, the project website and
  presentation version strings, and the headers of `run_simulation*.m` and
  `main_coupled_simulation.m` updated.

---

## [8.6.5] — 2026-10-20

Example-case release.  **No change to the solver, and no change to the results of
existing cases.**

### Added

- **Realistic *E. coli* single-flagellum case** (`configs/config_ecoli_single_flagellum.json`, `geometry/generate_ecoli_geometry.m`, `geometry/ecoli_body_capsule.csv`, `geometry/ecoli_flagellum_helix.csv`).  Custom-geometry case with a 2.5 µm × 0.88 µm sphero-cylinder body (`N_body = 300` requested, laid out by `capsule_body_building`, 326 points generated; neutrally buoyant, `T_obs = 1 s`) and one left-handed helical filament at the tail pole (pitch 2.22 µm, helix diameter 0.41 µm, 3 turns, ≈ 7.6 µm contour, filament radius 10 nm, EI = 3.5 pN·µm²), driven by a `free_torsion` root with a 1300 pN·nm / 300 Hz linear torque-speed motor.  The helix amplitude is ramped over the first half turn so the root leaves the pole along the cell axis.

### Compatibility

- No source file under `src/` changed.  Configurations from 8.6.3 run unchanged.

### Documentation

- `project_version.m` bumped to 8.6.5 (2026-10-20); README (EN/ZH), USER_GUIDE
  (EN/ZH), theory (EN/ZH), `requirements.txt`, the project website and presentation
  version strings, and the headers of `run_simulation*.m` and
  `main_coupled_simulation.m` updated.  The guides now point at the E. coli case.

---

## [8.6.3] — 2026-10-19

Animation-layout refinement.  **No change to simulation results.**  The animation
video changes (window ranges, body colouring, arrangement of the panels), so videos
rendered by 8.6.3 are not pixel-comparable with those of 8.6.2.

### Changed

- **Co-moving windows (c, d) are no longer cubes / squares.**  Along each lab axis
  the box extent is the reach of the body and all pili over all frames (relative
  to the centre of mass) plus an 8 % margin (at least 0.1 μm), centred on the
  occupied range (constant offset from the centre of mass), so the bacterium fills
  the window and its shape is visible.  The scale is equal on X, Y and Z (nothing
  is stretched or compressed); the axes rectangle of each window has the aspect
  ratio of its box (panel c: orthographic projection with a manual camera whose
  field of view is the shorter side of the projected box; panel d: XZ).
  Previously panel d was a square of half-width 1.15 x the largest distance from
  the centre of mass and panel c a cube of half-width 1.05 x that distance.
- **No scale bar in panels c and d** (panel b keeps its 2-μm scale bar,
  `opts.scale_bar_um`).
- **Adaptive panel layout.**  The four axes rectangles keep their natural aspect
  ratios and relative sizes (panels a, b: 810 x 1076 px; panels c, d: 440 x 478 px
  at scale 1), but their arrangement on the canvas is chosen automatically: all
  guillotine arrangements (nested side-by-side / stacked splits) of the panels in
  the reading orders a-b-c-d, a-c-b-d and a-c-d-b are evaluated; the one that fits
  the largest common scale factor (at most 2) of the four rectangles into the free
  canvas region, including the space for tick labels, titles and the gravity arrow,
  wins.  The slack is spread over the gaps (up to 120 px) and the group is centred.
  This removes the large empty areas that the former fixed "a | b | c over d"
  template left for elongated, flat or very wide / tall trajectories.  The chosen
  arrangement and scale are written to the log, for example
  `Layout H(a,H(b,V(c,d))): scale 0.98` (H: side by side, V: stacked).

### Added

- **Roll marker**: two dark longitudinal stripes on the body (a wide one, ±30°,
  and a narrow one, ±15°, on opposite sides; smooth profile, faded towards the
  end caps), stored as per-vertex colour and opacity of the body surface, so a
  rotation about the long axis is visible in every panel at no extra cost per
  frame (custom bodies: long axis = principal axis of the vertex cloud).  The
  legend entry reads "Body (dark stripes mark roll)".
- `render_motion_animation` returns `info.window_offset_um` (centre of the
  co-moving box relative to the centre of mass, 1 x 3, μm); the log reports the box
  half-extents and the offset.

### Compatibility

- `info.window_half_um` is now a 1 x 3 vector `[hx hy hz]` (half-extents of the
  co-moving box, μm) instead of a scalar.
- Simulation results, `simulation_results.mat`, CSV and all figures other than the
  motion animation are unchanged.  The static preview and summary figures keep the
  uniform body colour.
- No new dependency; the animation still needs MATLAB R2018b or newer.
- The Linux `Motion JPEG AVI` branch could not be exercised on the Windows
  development machine.

### Documentation

- `project_version.m` bumped to 8.6.3 (2026-10-19); README (EN/ZH), USER_GUIDE
  (EN/ZH), theory (EN/ZH), `requirements.txt`, the project website and presentation
  version strings, and the headers of `run_simulation*.m`,
  `main_coupled_simulation.m`, `render_motion_animation.m` and
  `render_animation_from_results.m` updated.  The animation descriptions now
  cover the per-axis window ranges, the roll-marker stripes, the missing scale bar
  in c and d and the automatic arrangement.

---

## [8.6.2] — 2026-10-18

Publication-oriented animation layout.  **No change to simulation results.**  The
animation video itself changes (layout, size, content of the panels), so videos
rendered by 8.6.2 are not pixel-comparable with those of 8.6.0.

### Changed

- **Four-panel animation on a fixed 16:9 canvas** (2560 x 1440 logical px) in
  `src/visualization/render_motion_animation.m`: **a** 3-D overview (now
  orthographic; the axes rectangle is sized to the projected data box), **b** XZ
  overview, **c** co-moving 3-D window and **d** co-moving XZ window.  The YZ panel
  is removed.  Previously the video size depended on the data (for example
  3840 x 1524 px).
- **Co-moving windows (c, d)** follow the centre of mass with constant magnification.
  Panel d has a fixed half-width of 1.15 x the largest distance from the centre of
  mass to any body or pilus point over all frames (at least 1.1 x the scale-bar
  length); panel c is an orthographic cube of half-width 1.05 x that distance around
  the centre of mass, zoomed to fill its axes (manual camera).  Both have an axes
  box with ticks, tick labels, a grid and axis labels in lab-frame coordinates, so
  the grid slides past the bacterium.  The centre-of-mass trajectory is drawn in
  panels a and b only; the gravity arrow only in panel a.
- Panel titles **a–d** (bold letter plus short descriptor) sit flush on the
  upper-left corner of each image; the time stamp, playback note and legend are
  horizontally centred on the canvas.
- Typography and line widths for the larger canvas: ticks 16 pt, axis labels 18 pt,
  legend and notes 16 pt, time stamp and panel titles 24 pt; pilus lines 2.0 pt,
  trajectory 1.5 pt, scale bars 4 pt.

### Added

- **Legend** (body, every pilus, centre of mass, trajectory), drawn manually on the
  canvas (90 px between entries, rows centred, automatic wrapping).
- **2-μm scale bars** in panels b, c and d (`opts.scale_bar_um`, default 2).
- **Global time stamp** *t* = ... s and the **playback speed** relative to real time
  (simulated seconds per video second = frame interval x `animation_fps`; labelled
  "slow motion" or "accelerated").
- `render_motion_animation` returns `info.window_half_um` and `info.playback_speed`;
  the log reports the co-moving window size and the playback speed.
- `add_gravity_indicator_outside` accepts an optional `head_scale` argument
  (default 1, unchanged behaviour).

### Documentation

- `project_version.m` bumped to 8.6.2 (2026-10-18); README (EN/ZH), USER_GUIDE
  (EN/ZH), theory (EN/ZH), `requirements.txt`, the project website and presentation
  version strings, and the headers of `run_simulation*.m` and
  `main_coupled_simulation.m` updated.  The animation descriptions (three panels
  with a YZ view) in README, USER_GUIDE and the website are replaced by the
  four-panel layout; USER_GUIDE gains a description of the panels, the playback
  speed and `opts.scale_bar_um`.

### Compatibility

- Simulation results, `simulation_results.mat`, CSV and all figures other than the
  motion animation are unchanged.
- The animation code now needs MATLAB R2018b or newer (`vecnorm`, `max(..., 'all')`); no new toolbox.
- The encoded frame size follows the display scaling (for example 3840 x 2160 at
  150 % scaling for the 2560 x 1440 canvas); the codec-limit clamp is unchanged.
- The Linux `Motion JPEG AVI` branch could not be exercised on the Windows
  development machine.

---

## [8.6.0] — 2026-10-17

Animation-rendering performance release.  **No change to simulation results or
to the video content** (resolution, frame count, frame rate and quality are
unchanged).

### Changed (performance)

- **Motion-animation rendering no longer rebuilds the scene on every frame.**  The
  3-panel animation (3D + XZ + YZ) used to call `cla` and re-create the body
  surface, every pilus, the centre-of-mass marker, the trajectory, four lights,
  the labels, titles, limits and border lines on all three axes for each frame.
  The new `src/visualization/render_motion_animation.m` creates all graphics
  objects once and only updates their `XData`/`YData`/`ZData` (or `Vertices`) and
  the 3-D title string per frame.  Measured on a 101-frame test run (capsule, two
  pili, 3840 x 1524 px video): 28.3 s -> about 13 s total, i.e. about 280 ms ->
  about 100 ms per frame; the per-frame cost no longer grows with the frame index.
  The remaining cost is dominated by `getframe` and `writeVideo`.
- The axis-limit pre-scan keeps running minima/maxima instead of growing
  concatenated coordinate arrays (identical limits).

### Added

- `src/visualization/render_motion_animation.m`: single shared implementation of the
  animation, now used by both `run_simulation.m` and `run_simulation_HPC_Linux.m`
  (previously two copies).  Releases its figure and `VideoWriter` itself on error.
- Animation phase of `simulation.log` now ends with a per-stage timing table
  (axis-limit pre-scan, figure setup, data update, `drawnow`, `getframe`, frame
  resize, `writeVideo`).
- `post_processing/render_animation_from_results.m`: re-renders the animation from a
  saved `simulation_results.mat` without re-running the ODE (optional `fps` / `skip`,
  defaults from the case folder's `config.json`); writes
  `bacteria_3d_motion_rerender.mp4` (`.avi` on Linux) and never overwrites the
  original video.
- `src/visualization/add_gravity_indicator_outside.m`: the gravity-arrow helper, moved
  out of the two run scripts (the local copies were identical).  `verify_migration.m`
  lists both new visualization files.

### Fixed

- The XZ and YZ animation panels re-drew their two black border lines every frame with
  `HandleVisibility = 'off'`, so `cla` never removed them and N duplicate lines
  accumulated by frame N.  They are now drawn once.

### Fixed (documentation)

- **The Ito drift of the body-only Brownian scheme is not "identically zero in an
  unbounded fluid".**  The statement in `docs/theory.md` / `docs/theory.zh-CN.md`
  section 10.2, the project website, the presentation and the v8.5.0 notes
  (`CHANGELOG.md`, `requirements.txt`) only holds for the rotational rows.  The
  translational rows equal `kBT * sum_{j,n} eps_ijn (M_TR)_nj`, which is non-zero
  when the tracking point `X_cm` is not the hydrodynamic centre (about 0.02 um/s
  for the default capsule with one 3 um pilus; zero for the pilus-free capsule
  and for symmetric pili arrangements).  It averages to zero over orientations.
- The code is unchanged: the force-realisation scheme integrates the actual
  `M(theta(t))` and therefore reproduces this term (Monte-Carlo test of the
  scheme as coded on a 2-D rigid body with `M_TR != 0`: mean displacement
  rate within 1.3 standard errors of `kBT R J b`; zero for `M_TR = 0`).  Only the
  documentation is corrected.

### Documentation

- `project_version.m` bumped to 8.6.0 (2026-10-17); README (EN/ZH), USER_GUIDE
  (EN/ZH), theory (EN/ZH), `requirements.txt`, the project website and
  presentation version strings, and the headers of `run_simulation*.m` and
  `main_coupled_simulation.m` updated.  README and USER_GUIDE gain notes on
  `render_animation_from_results`, the shared animation renderer and the
  per-stage timing log; `docs/index.html` lists the new files.

### Compatibility

- Video content, resolution, frame count, frame rate and `Quality = 100` are unchanged.
  Decoded frames of the new and the old renderer were compared for a capsule case with
  two pili, a pili-free capsule and a custom-mesh case: identical or within 4 grey
  levels (PSNR >= 92 dB over all 101 frames; first and last frame pixel-identical).
- `run_simulation.m` now selects the video profile by platform like the HPC script
  (`MPEG-4` .mp4 on Windows/macOS, `Motion JPEG AVI` .avi on Linux); previously it
  failed on Linux.  The final `Video  :` log line reports the actual file name.
- The Linux AVI branch could not be exercised on the Windows development machine.

---

## [8.5.2] — 2026-10-11

Calibration of the pilus hydrodynamic radius.  **Changes default results**
(see Compatibility).

### Fixed

- **The pilus regularisation length was `eps = r_pili` with no calibration.**
  For the Cortez (2018) 7/2 blob a straight rotlet-segment line has the axial
  twist drag `xi_t = 2*pi*mu*eps^2` per unit length, i.e. **half** the no-slip
  cylinder value `4*pi*mu*r^2`; the translational mobility of the Stokeslet
  line was 7.9 % (perpendicular) and 9.4 % (parallel) above the slender-body
  theory of a cylinder of radius `r` (mid-filament, `r/L = 3e-3`, 100
  segments; the reference is the RPY line with `a_hat = r*e^(3/2)/4` of Maxian
  et al.).  Both numbers were measured with `build_global_resistance` in
  `output/benchmark_v850/radius_check.m`.
- The pilus regularisation length is now `eps_pili = beta*r_pili` with the new
  calibration factor `beta = hydro_radius_factor`, default `sqrt(2)`: the twist
  drag equals the cylinder value exactly and the translational mobility is
  2.4 % / 2.9 % above the reference.  `beta = sqrt(e)` reproduces the
  translation exactly (twist drag then 36 % above the cylinder), `beta = 1`
  reproduces 8.5.1.

### Changed

- `init_coupled_params.m`: new name-value pair `hydro_radius_factor` (default
  `sqrt(2)`); `params.epsilon` is now the pilus regularisation length
  `beta*r_pili/L_0`; new fields `params.slenderness = r_pili/L_0` and
  `params.hydro_radius_factor`.  The bending stiffness `EI` and the body-pilus
  root gap `dis` still use the physical `r_pili`, so the geometry is unchanged.
- New optional JSON key `material.hydro_radius_factor` (also written to
  `config.json`) and `cfg.hydro_radius_factor`; validated as a positive finite
  scalar (a warning outside [1, 3]).  `configs/config_example.json` lists it.
- The pre-flight motor stall diagnostic in `validate_simulation_params.m` uses
  the calibrated twist drag; the `[PARAMS]`/`[INIT]` log lines report the factor.
- `init_coupled_params.m`: removed an unused assignment (`ic_input = {}` in the
  no-pili branch; code-analyzer NASGU, no behavioural change).

### Validation

- With `hydro_radius_factor = 1` the one-pilus sinking case (`D1`) and the motor case
  (`D4`) are bit-for-bit identical to 8.5.1 (single-thread BLAS); the HPC fork
  `run_simulation_HPC_Linux.m` reproduces `run_simulation.m` bit-for-bit for `D1`
  with the default factor (run on Windows).
- With the default factor the one-pilus sink distance changes by 0.13 %
  (`-0.0210405` to `-0.0210139` in the case's dimensionless units); in the motor
  example case (`D4`) the body counter-rotation angle over the run changes by a
  factor of 1.95 (0.0976 to 0.190 rad).
- The numerical check reproduces `xi_t/(2*pi*mu*eps^2) = 1.0000` and the
  analytic result `eps = sqrt(e)*r` for matching the translational mobility
  of the RPY line (perpendicular and parallel).

### Documentation

- `project_version.m` bumped to 8.5.2 (2026-10-11); README, USER_GUIDE, theory
  (EN/ZH), `requirements.txt`, the project website and presentation version
  strings, and the headers of `run_simulation*.m` and
  `main_coupled_simulation.m` updated.  The USER_GUIDE gains a note on the
  optional `material.hydro_radius_factor` key.

### Compatibility

- **Default results change** for pili-bearing cases: pili hydrodynamics changes
  by a few per cent, and the axial twist drag doubles.  Set
  `hydro_radius_factor = 1` to reproduce 8.5.1.  Pilus-free runs are unaffected.
- Previous Brownian and sinking batches were computed with `beta = 1`.

---

## [8.5.1] — 2026-10-10

Numerical-performance and maintenance release.  **No physical or
discretisation change**: the discrete equations, error control, tolerances and
Brownian model are those of 8.5.0.  Deterministic runs are bit-for-bit
unchanged.

### Changed

- **Lagged Jacobian for Brownian runs (`src/core/LaggedJacobian.m`,
  `run_simulation.m`, `run_simulation_HPC_Linux.m`).**  A Brownian run restarts
  `ode15s` at every `dt_B` node, and every restart paid a dense
  finite-difference Jacobian (`n_ode` extra RHS evaluations — 70–85 % of the
  ODE time of a four-pilus case, measured).  The Jacobian only enters the
  modified-Newton iteration matrix `I - gamma*h*J`; a stale `J` changes the
  convergence rate of that iteration, not the implicit equations it solves
  nor the local error estimate.  The last exact finite-difference Jacobian is
  therefore handed back to `ode15s` at the start of the next segments.  A
  fresh exact Jacobian is computed (i) when `ode15s` requests `J` a second
  time inside the same call (Newton failure), (ii) after an Euler-basis
  rotation (`adaptive_basis_check`), (iii) when the cache is
  `cfg.jacobian_max_age` Brownian nodes old (default 50, so a cached Jacobian
  serves at most the 49 following nodes; invalid values raise
  `LaggedJacobian:InvalidMaxAge`).  A segment whose lagged
  integration throws an error, or stops prematurely without an event or
  timeout, is re-run from the same state with the default Jacobian.
- New runtime-only options (not stored in JSON, like `cfg.silent`):
  `cfg.jacobian_reuse` (default `true`; `false` reproduces 8.5.0
  bit-for-bit) and `cfg.jacobian_max_age`.  The counters are written to
  `simulation.log` (`[JAC]` lines) and `results.jacobian_reuse`
  (`n_fresh`, `n_lagged`, `n_on_demand`, `n_age_refresh`, `n_invalidated`,
  `n_fallback`, `n_fd_rhs`).  Only Brownian runs use the cache.
- `src/utils/verify_migration.m` lists `src/core/LaggedJacobian.m` in the
  expected core files (it reported the new file as `UNEXPECTED` before) and
  uses `datetime` instead of `datestr(now)`.
- MATLAB code-analyzer clean-up in `run_simulation.m` and
  `run_simulation_HPC_Linux.m` (no behavioural change): obsolete
  `%#ok<NASGU>` / `%#ok<AGROW>` suppressions removed, two dead assignments
  (`anim_fig_w/h = new_log_w/h` after the codec-limit resize, never read
  again) removed, and the documented `output_dir = ''` default marked with
  `%#ok<NASGU>`.  `checkcode` reports 0 messages for both files,
  `LaggedJacobian.m` and `verify_migration.m`.

### Fixed

- `CHANGELOG.md` structure: the missing `## [7.9.1] — 2026-06-22` heading
  (date taken from the 7.9.1 commit) was restored, which also removes a
  duplicated `### Documentation` heading inside the 7.9.2 entry; nested-list
  indentation of the 5.5.5 entry and two multiple-blank-line spots were
  corrected.

### Added

- Markdown lint configuration: `.markdownlint.json` (`MD024` checked among
  sibling headings only, as required by the Keep-a-Changelog layout; `MD013`
  line length off) and `.markdownlint-cli2.jsonc` (third-party
  `MATLAB Add-On/` documentation and generated `output*/` folders are not
  linted).  `markdownlint-cli2` reports 0 issues for the project documents.

### Documentation

- `project_version.m` bumped to 8.5.1 (2026-10-10); README, USER_GUIDE and
  theory (EN/ZH; new `docs/theory.md` §10.3), `requirements.txt`, the project
  website and presentation version strings, and the `run_simulation*.m` /
  `main_coupled_simulation.m` headers updated.

### Validation (Windows, MATLAB R2026a, i9-13900KF, 1 BLAS thread per worker)

- **Deterministic regression** (sinking 1 and 4 pili, traveling-wave active
  moment, free-torsion motor, shear flow): outputs of 8.5.0 (git `HEAD`),
  8.5.1 with `jacobian_reuse=false` and `true` are bit-for-bit identical.
- **Brownian, same seed, 100 nodes (5 s)**: `jacobian_reuse=false` is
  bit-for-bit identical to 8.5.0 (1 and 4 pili).  `true` vs `false`:
  max centre-of-mass difference 0.007–0.14 nm (relative to path length
  ≤ 2e-5), body rotation ≤ 2.1e-4 rad, pilus-node position ≤ 0.72 nm — far
  below the Brownian displacement and well inside `RelTol = 1e-4`.  Cases:
  1 pilus, 4 pili (N_body 200/5 seg·µm⁻¹ and 300/10), motor + Brownian,
  active moment + Brownian, shear flow + Brownian, soft pili (E = 1e5, 1e6 Pa).
- **ODE wall-time speed-up** (same workload, `false` / `true`): 4 pili
  3.16× (N200/S5) and 3.66× (N300/S10), motor 2.8×, active moment 1.95×,
  soft pili 1.9–2.0×, 1 pilus 1.4×, shear 1.4×.  Soft pili (E = 1e5 Pa)
  trigger 24 basis rotations in 100 nodes; the cache is invalidated each
  time and the run still completes with no fallback.
- **Long runs** with the lagged Jacobian: 1 pilus, 200 s (4000 nodes) and
  4 pili, 60 s (1200 nodes) completed with 4 and 3 basis rotations and
  0 fallbacks.
- **Statistical equivalence**: 20 seeds, 1 pilus, 10 s — paired differences of
  the MSD-fitted D and of V_z between the two modes are ≤ 4e-5 (relative) and
  1.5e-5 µm/s, i.e. ~5e-5 of the between-seed spread.  Bootstrap 95 % CIs
  (20000 resamples) of the group means coincide: D = 0.2655 [0.2265, 0.3095]
  vs 0.2655 [0.2265, 0.3092] µm²/s; the CI of the paired mean difference of D
  is [-2.1e-6, 2.3e-6] µm²/s.  Only the one-pilus configuration was tested
  statistically; the four-pilus case is covered by the same-seed comparisons.
- **Fault injection** (patched copy): an error thrown from every 7th Jacobian
  request triggers the fallback branch 16 times in 100 nodes and the run
  completes with deviation of the same size as above; replacing every lagged
  Jacobian by a random garbage matrix (entry-wise ×(1+2·randn)) still gives
  the same trajectory (on-demand refresh every node).
- After the maintenance edits, with single-thread BLAS the deterministic
  sinking (1 pilus), motor and shear cases and the 100-node Brownian
  one-pilus case reproduce the validated results bit-for-bit, and animated
  runs of `run_simulation.m` and `run_simulation_HPC_Linux.m` complete and
  write the video.  Results are bit-reproducible only for the same BLAS
  thread count.
- A non-silent run (progress bar and console output enabled) was smoke-tested
  with the lagged Jacobian; all other validation runs used `cfg.silent=true`.
  Failed-step rates were measured only in the stand-alone integrator study
  (0.06 failed steps per node, 0 for the default Jacobian), not in the
  `run_simulation` validation runs.
- Newton-contraction analysis of the drift of `J` (4 pili, 20 s): spectral
  radius of `I - W_old^{-1} W_true` ≤ 0.17 for cached Jacobians up to 64 nodes
  old and ≤ 0.34 for `h` up to 1e-2 (dimensionless), versus actual `ode15s`
  steps of median 1.7e-4 and p99 1.7e-3.

### Compatibility

- Brownian trajectories are not bit-for-bit identical to 8.5.0 unless
  `cfg.jacobian_reuse = false` (differences ≈ 1e-5 relative, see above).
- `run_simulation_HPC_Linux.m` carries the identical patch; its Linux
  execution was not re-run here (a Windows smoke test of the fork passed).
- No new MATLAB toolbox dependency.

---

## [8.5.0] — 2026-09-04

Physics/algorithm release for the second "major" finding of the full
mathematical/physical code review: the body-only Brownian displacement kick
followed by root re-projection injected a spurious elastic kink at the pilus
root on every Brownian node.  The stochastic increment is now realised as a
constant generalised thermal force acting through the coupled ODE.  Also
closes the remaining minor geometry finding (type-2 parabolic template not
anchored at the origin).

### Changed

- **Brownian dynamics — force realisation instead of displacement kicks
  (`run_simulation.m`, `run_simulation_HPC_Linux.m`, new
  `src/physics/sample_brownian_body_force.m`, `ode_coupled_rhs.m` Step 8a,
  `compute_effective_body_brownian.m`, `init_coupled_params.m`).**  At the
  start of every Brownian interval `[T_last, T_last + dt_B]` the 6×6
  effective body mobility `M_eff` is probed once, the frozen-mobility
  Euler–Maruyama increment `dU = L_eff * eta * sqrt(dt_B)` is drawn (six
  normal deviates from the isolated stream, as before), and the constant
  generalised force `f_th = M_eff \ (dU / dt_B)` (dimensionless units
  `F_0`, `F_0*L_0`) is stored in the new field `params.brown_force_dimless`.
  `ode_coupled_rhs` subtracts it from `R_global(1:6)` with the same sign
  convention as gravity for the whole interval; the field is reset to zero
  when the interval closes and after the main loop.  The body is therefore
  moved by the coupled solver and the pili follow it continuously through
  their hydrodynamic and elastic coupling — the root constraint is never
  violated, so the post-kick `enforce_clamped_geometry` re-projection (which
  rotated segment 1 onto the moved body while segments `2..N` kept their
  orientation, creating an `O(|dTheta|/ell)` curvature jump that relaxed on
  the stiff `ell^4` time scale and forced `ode15s` into very small steps
  after every node) is no longer needed.  A basis-rotation or timeout
  restart inside an interval keeps the already-sampled force.
  `compute_effective_body_brownian` gains a fourth output `M_eff_dimless`.
- **Ito drift is now implicit.**  Expanding the constant-force ODE over one
  interval gives `E[dU_real] = kBT div(M) dt_B + O(dt_B^{3/2})` because
  `E[f_th f_th'] dt_B = 2 kBT M^{-1}` (docs/theory.md §10.2), so the
  scheme is a consistent weak first-order integrator of the body-only Ito
  SDE.  The explicit centred finite-difference `kBT*div(M_eff)` correction
  introduced in 7.5.0 (`compute_body_ito_drift` / `perturb_body_pose`, 12
  extra coupled solves per node) is removed from both run entry points;
  adding it again would double count the drift.  Per-node cost drops from
  13 to 1 effective-mobility probe.
- The `brown_kick_applied` post-kick sample override in the output-time
  loop is removed (the ODE solution is already the realised state).  Log
  lines change from `[BROWN] Kick #n ...` to `[BROWN] Interval #n start ...`
  (sampled `|F_th|`, `|tau_th|`, frozen-mobility target `|dX|`,
  `|dTheta|`) and `[BROWN] Interval #n closed ...`; the `[PARAMS]` Brownian
  line reports the scheme.  `results.brownian_step_sizes_s` keeps its
  meaning (completed interval lengths).
- `docs/theory.md` / `theory.zh-CN.md` §10 restructured (10.1 effective
  mobility, 10.2 force realisation with the kink and implicit-drift
  derivations, plus the remark that in an unbounded fluid
  `M_eff = R M_shape R'` so the body-pose Ito drift `kBT div(M)` vanishes
  identically); `USER_GUIDE` §8 (EN/ZH), READMEs and `project_version.m`
  bumped to 8.5.0.
- Documentation catch-up for 8.2.0–8.5.0: `requirements.txt` release
  notes, the README core-file tables (`compute_pili_curvature.m`,
  `compute_pili_elastic_moments.m`, `sample_brownian_body_force.m`,
  `background_flow_velocity.m`), the version/documentation-release headers
  of `main_coupled_simulation.m` / `run_simulation*.m`, and the project
  website (`docs/index.html` release cards, Brownian section rewritten for
  the force realisation and implicit Ito drift; `docs/presentation.html`
  Brownian slides and module notes) were still at 8.1.0 and described the
  displacement-kick scheme; all now reflect 8.5.0.

### Fixed

- **`pili_type_building.m` — type-2 parabolic template anchored at the
  origin.**  The least-squares arc-length inversion polynomial carries a
  residual constant term (`-0.004127`), so `x(s = 0) != 0` and the first
  template node sat `~4.1 nm / K` off the origin, contradicting the
  origin-anchored template assumption of `pili_building`,
  `custom_pili_building` and `enforce_clamped_geometry`.  The inverted curve
  is now shifted by its own `s = 0` value (uniform in `x`, node spacing
  unchanged).

### Behaviour changes / migration notes

- Brownian runs change quantitatively even for a fixed `brownian.rng_seed`:
  the same six deviates per interval are consumed, but the increment is
  realised through the coupled trajectory (mobility evaluated along the
  path, pili deforming with the body) instead of on the frozen mobility.
  Statistics (MSD, MSAD, `D_rot`) agree with <= 8.4.0 to the `O(dt_B)`
  weak-order accuracy of both schemes; trajectories do not.  Deterministic
  runs (`brownian.enable = false`) are bit-for-bit unchanged.
- Brownian runs should be noticeably faster: no finite-difference drift
  probes and no post-node stiff transient.  Configurations that had
  shortened `dt_B` or refined `N_segments` solely to tame the root kink can
  revert to their physically motivated values.
- Parabolic (type-2) pili are translated by `0.004127 / K` um along the
  template tangent relative to <= 8.4.0.

## [8.4.0] — 2026-09-04

Physics-correctness release for the first "major" finding of the full
mathematical/physical code review: the `free_torsion` root twist equation
double-counted the interior Kirchhoff twist moment.  The remaining major
finding (body-only Brownian kick creating a spurious root kink) was fixed in
8.5.0.

### Fixed

- **`compute_pili_elastic_moments.m` / `coupled_active_moment.m` — root
  row zeroed for `free_torsion` as well as `clamped`.**  Row `i` of the
  pilus moment balance is the balance of the rod portion `[s_i, 1]` about
  node `i`; the internal moment `M(s_i)` that enters it is the moment
  exerted by `[0, s_i)` on `[s_i, 1]`.  At the root node `[0, s_1)` is the
  cell body, so `M(0)` is the torque the body applies through the
  bearing/motor — the motor torque added in `ode_coupled_rhs` Step 9 (zero
  for a free bearing).  The discrete `kappa_3` available at node 1 is
  evaluated between segments 1 and 2 and is already the internal moment
  used by the node-2 balance; adding it to the root row too (`kappa3(1)
  preserved`, since v7.x) meant that subtracting the node-2 row left
  `tau_motor = drag of segment 1 alone`, so the whole filament span up to
  the rate at which the first segment's drag balanced the motor.  Verified
  on a straight 2 um filament (`N = 10`, 1 pN nm, no speed limit): steady
  `omega_rel` was `8.47e4 rad/s` before and is `8.84e3 rad/s` after the
  fix (ratio exactly `N`); the steady twist profile is linear from node 2
  to the tip with `GJ*kappa_3(s_2) = tau*(1 - 1/N)`; a free bearing with a
  `helix` intrinsic twist and no motor now gives zero net axial
  hydrodynamic torque (`sum(dpsi) / max|dpsi| = 1e-5` at `t = 0`) and
  relaxes to rest at `kappa_3 = kappa_d3`.  The integrated active moment is
  likewise an internal moment, so its root `d3` component is dropped for
  the same reason (a purely internal twist wave can no longer spin a
  free-bearing filament).  Both functions keep their `root_coupling_k`
  argument for interface stability.
- **`validate_simulation_params.m` — stall diagnostic drag coefficient.**
  The steady-state estimate now uses `xi_t = 2*pi*mu*eps^2*L`, the axial
  rotational drag of a straight regularised-rotlet-segment line with the
  Cortez (2018) 7/2 blob (`eps = r_pili`), which the simulator reproduces to
  `< 1e-3`; the former `4*pi*mu*r^2*L` (no-slip cylinder) under-predicted
  the simulated speed by exactly 2.  The `free_torsion` + `motor off`
  warning now describes the frictionless-bearing behaviour.

### Changed

- Comments in `ode_coupled_rhs.m` (BC RHS block) and
  `build_coupled_lin_sys.m` (Step 7) describe the `psi_1` row as the
  whole-rod boundary-torque balance.
- `docs/theory.md` / `theory.zh-CN.md` §8 gain a "Root twist row for
  `free_torsion`" subsection with the derivation and the `xi_t` note;
  `USER_GUIDE` §7 (EN/ZH), READMEs and `project_version.m` bumped to 8.4.0.

### Behaviour changes / migration notes

- All `free_torsion` runs change quantitatively.  With a motor, steady
  rotor speeds (and the induced body counter-rotation, swimming speed and
  flow fields) drop by roughly the segment count `N` relative to <= 8.3.0;
  motor configurations tuned to reproduce a target rotation rate under the
  old code must be re-tuned (the pre-flight stall estimate is now a
  reliable guide).  Without a motor the root acts as a true frictionless
  bearing.  `clamped` runs are unaffected.

## [8.3.0] — 2026-09-04

Discretisation-accuracy and geometry-consistency release addressing the two
"medium" findings of the full mathematical/physical code review: the
half-segment offset of the Kirchhoff curvature stencil and the handedness
(chirality) flip introduced by some attachment-site transforms.  The
hydrodynamic kernels, the coupled linear system and the Brownian stepping
are untouched.  The two remaining major findings (free-torsion root moment
double counting, body-only Brownian kick kink) are still tracked separately.

### Added

- **`src/physics/compute_pili_curvature.m`** — single shared operator that
  returns the raw nodal curvature/twist `[kappa_1, kappa_2, kappa_3]`
  (`[(N_k+1) x 3]`, units `1/L_pili`) from the segment directors.  Nodal
  frames come from `compute_node_directors` (SVD polar mean); the arc-length
  derivative at interior node `i` is formed from the two **segment** frames
  meeting at that node, `(d_seg(i) - d_seg(i-1)) / ds`, which is centred at
  the node.  Both factors of the projection are therefore collocated and the
  stencil is second-order accurate.  End nodes replicate the nearest interior
  stencil (node 1 is overwritten by the root BC; node N+1 is unused).

### Changed

- **`compute_pili_elastic_moments.m` — node-collocated curvature stencil.**
  Previously the nodal directors were differenced backwards,
  `(dN(i) - dN(i-1)) / ds`, which is centred at `s_i - ds/2`, i.e. half a
  segment away from the nodal director it was projected onto.  On an
  analytic constant-Darboux-vector frame the observed convergence order was
  0.98; with the new operator it is 2.00 and the absolute error at N = 40
  drops from `3.4e-2` to `3.7e-4` (|Omega| ~ 2.6).  The function now calls
  `compute_pili_curvature` and only applies the intrinsic-curvature
  subtraction, root-row rules and the `1/(1+sigma)` twist factor.
- **`init_coupled_params.m` — `kb0` via the same operator.**  Both the
  parametric and custom branches replace their duplicated stencil code by a
  call to `compute_pili_curvature`, so `kappa_live - kb0` at `t = 0` remains
  bit-consistent (`1.7e-13` over pilus types 1-5 in the regression check) and
  `frozen` pili still start at exact zero-moment equilibrium.
- **`pili_building.m` — attachment-site maps are proper rotations.**
  `pili_location_single` now applies `p_body = R * p_template + SPoint` with
  `R = [c1, c2, c1 x c2]` (`c1` = outward site direction = image of the
  template tangent `+x`, `c2` = image of the template bending direction
  `+y`) and asserts `det R = +1`.  The legacy per-site coordinate swaps had
  `det = -1` at sites **1, 2, 5, 7, 10, 11, 14** and mirrored every
  out-of-plane feature: a right-handed conical-spiral template (type 5)
  became left-handed at those seven sites, and the sense of the
  `deflection_x` twist was inverted, so the same `Pili_Matrix` column
  produced opposite chirality depending only on the attachment site.  The
  new rotations agree with the legacy maps on the template bending plane
  `z = 0`, so all planar templates (types 1-4 with `deflection_x = 0`) are
  placed **identically** to before (regression difference 0.0 at all 14
  sites); only the out-of-plane component changes sign at the seven
  affected sites.
- **`pili_type_building.m` / `custom_pili_building.m` — right-handed `RX`.**
  The deflection matrix about the template tangent (`gama` / `deflection_x`)
  used a left-handed sign convention (documented as "project-internal").  It
  is now the standard right-handed rotation, consistent with `RZ`,
  `pili_directors` and `quat_to_rotmat`; `RM = RX * RZ` is a rotation with
  `det = +1`.  Verified: a parabola bending toward `+y` with
  `deflection_x = +pi/2` now bends toward `+z`.
- **`verify_migration.m`** — expected file lists extended with
  `background_flow_velocity.m` (added in 8.2.0) and
  `compute_pili_curvature.m`.
- **Documentation** — `docs/theory.md` gains a "Discrete curvature stencil"
  subsection and a paragraph on the attachment-site rotations; version
  strings bumped to 8.3.0 in `project_version.m`, READMEs, user guides and
  the theory document.

### Behaviour changes / migration notes

- Runs with **type-5 (conical spiral) pili** at sites 1, 2, 5, 7, 10, 11 or
  14, or with **non-zero `deflection_x`** at any site, will see the mirror
  image of the pre-8.3.0 initial geometry.  The pre-8.3.0 geometry was the
  unintended reflection; if a specific handedness was tuned against it,
  flip the sign of spiral parameter `c` (type 5) or of `deflection_x`.
- All planar configurations (types 1-4, `deflection_x = 0`) and all
  `frozen` / `zero` / `constant` / `helix` intrinsic-curvature settings are
  unaffected in their initial state.  Elastic moments away from the
  reference configuration change at `O(ds)` because of the improved stencil;
  trajectories at coarse `seg_per_um` will differ slightly from 8.2.0 and
  should be regarded as more accurate.

## [8.2.0] — 2026-09-04

Physics-correctness release addressing the "secondary" findings of the full
mathematical/physical code review (motor torque-speed law, flow-field
background flow, shape-parameter validation, unit documentation).  The
hydrodynamic kernels, elastic-rod discretisation and Brownian stepping are
untouched in this release; the major review findings (free-torsion root
moment double counting, body-only Brownian kick kink, half-segment curvature
stencil offset, attachment-site handedness) are tracked separately.

### Changed

- **`solve_motor_self_consistent.m` / `ode_coupled_rhs.m` — signed, asymmetric
  linear torque-speed law in the relative rotor speed.**  The former law
  `tau_eff = tau_motor*(1 - |omega|/omega_max)` with `omega = dpsi_1` was
  symmetric in the rotation sense (a rotor forced backwards saw the same
  torque drop as one spinning forwards) and measured the speed of segment 1
  in the lab basis with the `cos(theta_1)*dphi_1` part of the axial rate
  omitted.  The motor now follows
  `tau_eff = tau_motor*(1 - omega_rel/omega_max)` with
  `omega_rel = sign(tau_motor) * [(dpsi_1 + cos(theta_1) dphi_1) - Omega_body . d3_1]`,
  i.e. the full axial angular velocity of the rotor segment relative to the
  cell body, positive in the driving direction (Berg & Turner 1993; Berg
  2003).  Forced reverse rotation raises the torque above stall and speeds
  beyond `max_omega_Hz` make the motor a brake.  Because the law is linear,
  the reduced motor problem is solved exactly as a small linear system for
  any number of motors (the scalar piecewise root / damped fixed point are
  removed).  `ode_coupled_rhs` builds one observation vector per motor
  (`motor_obs`) instead of indexing `dpsi_1` columns.
- **Motor torque sign convention defined (`ode_coupled_rhs.m` Step 9).**
  `R_global` stores the negative of the applied generalised force, so the
  motor torque now enters its `psi_1` row as `-tau`: a positive
  `torque_pN_nm` drives the rotor in the right-handed sense about the
  outward root tangent `d3` (counter-clockwise viewed from the filament tip
  toward the cell).  Before 8.2.0 the sense was the opposite and undocumented;
  configurations that relied on the old rotation direction must flip the
  sign of `torque_pN_nm`.  The `|omega|` law made this sign irrelevant to the
  speed limiter, which is why it went unnoticed.
- **`compute_plot_flowfield.m` — background flow honoured.**  Force and
  torque densities are now solved from the velocity of the structure
  *relative to the ambient fluid* (`v - u_bg` for the translational block and
  the segment twist rates), matching `ode_coupled_rhs`; the ambient velocity
  is added back to every exterior grid point so the plotted field is the
  total fluid velocity.  Previously both steps ignored `background_flow`,
  yielding wrong force densities and a disturbance-only picture whenever a
  `uniform` / `shear_z` flow was active.
- **`coupled_background_flow.m`** delegates the translational profile to
  the new shared evaluator (behaviour unchanged).

### Added

- **`src/physics/background_flow_velocity.m`** — single evaluator of the
  translational background-flow profile (and the fluid angular velocity of
  the shear layer) at arbitrary points, used by both the ODE RHS and the
  flow-field post-processor so the simulated and visualised ambient flow
  cannot diverge.
- **`validate_simulation_params.m` — `shape_params` checks (local
  `validate_pili_shape_params`).**  Parametric pili of type 2-5 now fail
  validation with an explicit message when `Pili_Matrix` rows 5-7 are
  missing (NaN from `parse_pili_config`) or out of range: parabolic needs
  `K > 0` and `K*length_um <= 4` (fitted range of the arc-length
  polynomial); circular needs a subtended angle in `(0, 2*pi]`; sinusoidal
  needs finite `[a, b]`; spiral needs finite `[a, b, c]` and not `a = b = 0`.
  Previously such inputs surfaced only as an opaque "dimension mismatch"
  concatenation error inside `pili_type_building.m`.

### Documentation

- **`USER_GUIDE.md` / `USER_GUIDE.zh-CN.md` — units of pilus curvature and
  moment inputs.**  New subsection stating that pilus elasticity is
  nondimensionalised by the pilus length: `intrinsic_curvature.kappa`,
  `kappa_d2`, `kappa_d3` are in units of `1/L_pili` (helix conversion
  formulas given), while `active_moment.amplitude_Nm`, `torque_pN_nm` and
  `max_omega_Hz` are converted to `EI/L_pili` and `1/T_0` units.  §7 documents
  the new torque-speed law and the torque sign convention.
- `README.md` / `README.zh-CN.md` version notes; `project_version` -> 8.2.0.

## [8.1.0] — 2026-07-25

Minor feature release for batch post-processing kinematics: 3-D
translational MSD fitting with a \(0\sim T/2\) window, plus XZ-plane
rotational diffusion from MSAD.

### Added

- **`post_processing/postprocess_cases.m` — rotational diffusion (`D_rot`)**:
  mean-square angular displacement (MSAD) of the body long-axis angle in
  the XZ plane. Prefer \(\phi=\mathrm{atan2}(e_z,e_x)\) from `results.q`
  in `simulation_results.mat`; fall back to `tilt_rad` from
  `timeseries.csv`. Fit \(\mathrm{MSAD}(\tau)=2D_\mathrm{rot}\tau\)
  (no intercept) on \(0<\tau\le T/2\); export `D_rot_rad2_s` and
  `MSAD_R2` to Excel and the batch summary figure.

### Changed

- **3-D translational MSD fit window**: TAMSD still uses full observation
  lags \(\tau=dt,\ldots,\sim T\), but the ballistic–diffusive model
  \(\mathrm{MSD}=6D\tau+V^2\tau^2\) is fitted only on \(0<\tau\le T/2\)
  (replacing the former \(\sim T/4\) fit cutoff). Reported `D_trans` /
  `V_msd` remain non-negative-clamped.
- Batch summary figure layout expanded to \(2\times3\) panels to include
  \(D_\mathrm{rot}\).

### Documentation

- Update `README.md`, `README.zh-CN.md`, `USER_GUIDE.md`,
  `USER_GUIDE.zh-CN.md`, `requirements.txt`, both theory documents, the
  project website, presentation, runtime banners, and `project_version.m`
  to version 8.1.0 dated 2026-07-25.

### Compatibility

- No change to the hydrodynamic / elastic / Brownian simulation core.
- No new MATLAB, toolbox, Python, Conda, or pip dependency.
- Post-process Excel schema gains two columns (`D_rot_rad2_s`, `MSAD_R2`);
  consumers that hard-code column indices should be updated.

---

## [8.0.1] — 2026-07-21

Patch release hardening runtime numerical diagnostics after production
Brownian batch failures on dense mobility operators (`N_body ~ 300`).

### Fixed

- Dense weighted-mobility dissipation checks no longer call
  `eigs(...,'smallestreal')`. Dense operators use full `eig`; iterative
  `eigs` is retained only for sparse matrices with automatic fallback to
  dense `eig` on non-convergence (instead of aborting).
- Isolated `rcond < 1e-14` no longer hard-aborts. It warns as
  `NearSingularMatrix` and is recorded in `results.diagnostics`.
- Residual WARN-band checks no longer escalate into a singularity abort.
  `SingularMatrix` requires the joint gate:
  backward residual `> 1e-6` **and** `rcond < 1e-14`.
- Initial / post-basis health checks probe a gravity-or-unit RHS through
  `check_linear_solve`, so singularity is gated jointly with solve quality.
- Diagnostic snapshot field renamed from `rcond_error_threshold` to
  `rcond_near_singular_threshold`; log banners label it as near-singularity /
  joint-abort context rather than a standalone error threshold.
- `verify_migration.m` no longer expects a `tests/` tree (offline MATLAB
  `functiontests` suite is not shipped in this release layout).

### Documentation

- Update `README.md`, `README.zh-CN.md`, `USER_GUIDE.md`,
  `USER_GUIDE.zh-CN.md`, `requirements.txt`, both theory documents, the
  project website, presentation, runtime banners, and `project_version.m`
  to version 8.0.1 dated 2026-07-21.

### Compatibility

- **Runtime contract:** `numerical_diagnostics:SingularMatrix` now means the
  joint singularity gate failed (large residual **and** near-singular
  `rcond`). Isolated near-singularity warns via
  `numerical_diagnostics:NearSingularMatrix`.
- No new MATLAB, toolbox, Python, Conda, or pip dependency.
- Physical hydrodynamic / elastic / Brownian algorithms are unchanged from
  8.0.0 aside from the diagnostic gate behaviour above.

---

## [8.0.0] — 2026-07-13

Major mathematical-physics correctness and runtime-contract release following
a first-principles review of the coupled hydrodynamic, elastic, Brownian,
geometric, motor, timeout, and numerical-diagnostics algorithms.

### Corrected

- Enforce line-measure-weighted reciprocity for translation/twist coupling
  and unequal-length twist segments; use the same operator in simulation and
  flow-field reconstruction.
- Remove the non-functional per-pilus stiffness interfaces.  `Pili_Matrix`
  is exactly eight rows and all pili use the global material modulus and
  radius.
- Preserve full material-director objectivity under arbitrary initial body
  rotations, and derive custom-pilus reference geometry, directors, and
  intrinsic curvature from one equal-chord canonical centerline.
- Apply Brownian SPD repair in dimensionless coordinates, include the final
  remainder step, and isolate Brownian random draws from global/cache RNG
  state.  The stochastic model remains intentionally body-only.
- Reject non-orthogonal `shear_z` flow/shear directions as incompatible with
  incompressibility.
- Treat any integration-segment timeout as a normal partial termination:
  truncate preallocated arrays, save MAT/CSV/config/log data, return
  `termination_reason='timeout'`, and exclude the trajectory from batch
  steady-state statistics.
- Add finite-value, linear backward-residual, weighted reciprocity,
  dissipation, conditioning, and multi-motor convergence diagnostics.

### Validation

- Add MATLAB `functiontests` for hydrodynamic invariants, canonical geometry,
  initial-frame objectivity, uniform-stiffness contracts, Brownian remainder
  stepping, timeout persistence, strict shear validation, and numerical
  guards.
- Pass 38 MATLAB tests, five end-to-end families (sphere/standard,
  motor-driven, custom-geometry, seeded Brownian, and HPC), JSON
  configuration validation, Windows/HPC deterministic-entry agreement,
  migration verification, and whitespace checks.

### Documentation

- Update `README.md`, `README.zh-CN.md`, `USER_GUIDE.md`,
  `USER_GUIDE.zh-CN.md`, `requirements.txt`, both theory documents, the
  project website, presentation, runtime banners, and `project_version.m` to
  version 8.0.0 dated 2026-07-13.
- Document the body-only Brownian model boundary, weighted hydrodynamic
  operators, equal-chord canonical geometry, strict incompressible shear
  contract, timeout return lifecycle, and numerical diagnostic fields.

### Compatibility

- **Breaking:** `Pili_Matrix` must have exactly eight rows. The former ninth
  row and custom/named `bending_stiffness` overrides are rejected; all pili
  derive stiffness from global `material.E` and `r_pili`.
- **Breaking:** `shear_z.direction` and `shear_z.shear_axis` must be finite,
  non-zero, and orthogonal within `1e-12`; invalid configurations now fail
  pre-flight validation.
- **Runtime contract:** a segment timeout no longer aborts without a usable
  return. It saves the achieved trajectory and returns
  `results.completed=false` with `termination_reason='timeout'`; first-party
  batch and post-processing code excludes it from steady-state statistics.
- No new MATLAB, toolbox, Python, Conda, or pip dependency.

---

## [7.10.1] — 2026-07-12

Documentation-accuracy and API-contract consistency patch following a complete review of all first-party MATLAB scripts, functions, entry points, and frozen legacy snapshots. This release does not change the hydrodynamic, elastic, Brownian, motor, or geometry numerical algorithms.

### Fixed

- **`src/utils/odetpbar.m` return contract:** initialise `status = 0` at function entry so the documented success value is returned for initialization, normal-update, and completion callbacks.
- **`src/utils/verify_migration.m` inventory:** recognise `run_simulation_HPC_Linux.m` and `src/utils/project_version.m`; document the repository-root working-directory requirement and destructive `close all force; clear; clc` side effects.
- **Stable capsule-mesh error identifiers:** use `capsule_sphere_cap_points:InvalidSPN` and `capsule_sphere_mesh_cache:UnsupportedMeshOutput` with separate messages.

### Documentation

- Corrected and expanded MATLAB help headers across 49 first-party files: function signatures, state and matrix dimensions, SI/dimensionless units, lab/body/basis frames, root boundary conditions, Brownian-only symmetrisation and diagnostics, motor solve ordering, configuration round-trip limitations, file outputs, randomness, side effects, warnings, errors, cache paths, and legacy/live-API boundaries.
- Corrected current-demo metadata, JSON fallback wording, combined-cap full-sphere `SPN` semantics, sphere-cache paths, custom-geometry coordinate contracts, and Linux HPC `.avi`/heartbeat documentation.
- Updated `README.md`, `README.zh-CN.md`, `USER_GUIDE.md`, `USER_GUIDE.zh-CN.md`, `requirements.txt`, theory documents, website, presentation, runtime banners, and `project_version.m` to version 7.10.1 dated 2026-07-12.

### Validation

- Parsed every first-party MATLAB file and verified non-empty MATLAB help text.
- Passed configuration loading/validation for all three shipped JSON examples, JSON save/load round trip, quaternion/director/coordinate identities, body geometry, coupled-system assembly, finite RHS, Brownian SPD covariance, background flow, stable error identifiers, and PNG/PDF export.
- Passed five end-to-end smoke cases: standard active-pilus, free-torsion motor, custom geometry, seeded Brownian body-only, and `run_simulation_HPC_Linux`.
- `verify_migration` passed all eight directory checks; `git diff --check` reported no whitespace errors.

### Compatibility

- **No physical or numerical simulation algorithm changes.** Existing numerical results are unaffected except that previously failing utility/error paths now return their documented status or stable identifiers.
- No new MATLAB, toolbox, Python, Conda, or pip dependency.

---

## [7.10.0] — 2026-07-03

New HPC/Linux cluster execution capability, plus batch-driver bug fixes surfaced during the corresponding review. A dedicated review of `batch_repeat_brownian.m` for deployment on a Linux HPC partition (single-node-exclusive, distributed/network storage, MATLAB R2023a) found several Windows-only API calls that either hard-error or silently mis-size the worker pool on Linux, an `ensure_parpool` `JobStorageLocation` default that is a documented source of pool-startup races on clusters with network home directories, and no visibility into a still-running repeat's progress. This release adds a maintained fork — `run_simulation_HPC_Linux.m` / `batch_sweep/batch_repeat_brownian_HPC_Linux.m` — that resolves all of these for Linux while leaving the original Windows-oriented scripts completely untouched, and separately fixes two reproducibility/robustness bugs found in `batch_repeat_brownian.m` and `batch_run_from_json.m` during the same review that affect both platforms.

### Added

- **New files: `run_simulation_HPC_Linux.m`, `batch_sweep/batch_repeat_brownian_HPC_Linux.m`** — a maintained fork of `run_simulation.m` / `batch_repeat_brownian.m`, hardened for headless, single-node-exclusive Linux HPC cluster execution. Not auto-synced with the originals; the originals are unchanged and remain the Windows/interactive-oriented reference implementation.
- **Platform-aware `VideoWriter` profile (`run_simulation_HPC_Linux.m`)**: Windows/Mac keep `'MPEG-4'` (`.mp4`); Linux uses `'Motion JPEG AVI'` (`.avi`), MATLAB's built-in, dependency-free profile. `'MPEG-4'` is documented by MathWorks as Windows/Mac-only and throws `"The specified profile is not valid."` on any Linux system regardless of MATLAB release or display availability — with `generate_animation = true` this previously meant every repeat failed at the video-encoding step on Linux, after the physics integration had already completed. The output-directory "has prior data" file-extension check also recognises `.avi`.
- **Linux-native resource probing (`decide_pool_size` in `batch_repeat_brownian_HPC_Linux.m`)**: physical-core count now falls back to `lscpu -p=Core,Socket` (counts unique Core/Socket pairs, excludes hyper-thread siblings) when the Windows-only `wmic`/CIM/`NUMBER_OF_PROCESSORS` probes are unavailable; available memory now falls back to parsing `/proc/meminfo`'s `MemAvailable` field when the Windows-only `memory()` function is unavailable. Previously, on any Linux system, both probes always fell through to their weakest fallback: `memory()` throws immediately, defaulting the pool-sizing calculation to "assume 4 GB available" regardless of actual system RAM (silently capping the parallel pool at 1 worker on a 512 GB node), and the physical-core fallback incorrectly divided `feature('numCores')` by 2 even though that function already returns the physical (not logical) core count.
- **`cpu_resource_frac` / `mem_resource_frac` (replaces the single `max_resource_frac`)**: CPU and memory headroom are now capped independently (defaults `0.90` / `0.65`), since per-run memory is harder to predict than CPU usage and running out of RAM is the failure mode that can hang the whole batch (see OOM item below), not just slow it down.
- **Heavy-output sampling (`heavy_output_sample_idx`)**: `generate_animation` / `compute_flowfield` are forced off for every repeat except a small sample (first, middle, and last by default), regardless of what the inline config or JSON specifies. These features hold significant extra per-worker memory (frame buffers, video encoding, flow-field grids) and were implicated in an OOM-killer incident observed on a real HPC run (worker process count in `top` dropping abruptly mid-batch, followed by an indefinitely hung `parfor` — `parfor` has no built-in timeout or dead-worker recovery). Sampling reduces the likelihood of hitting this; it does not eliminate it, since `parfor` itself was kept (a `parfeval` + `fetchNext` timeout-based task queue would close this gap fully but was deliberately deferred).
- **Node-local `JobStorageLocation` (`ensure_parpool` in `batch_repeat_brownian_HPC_Linux.m`)**: on Linux, the local parallel pool's internal coordination/lock files are relocated from MATLAB's default location (which resolves under the user's home directory) to `tempdir()` — honouring `$TMPDIR` when the job scheduler sets it to node-local scratch — in a job-unique subfolder named from `SLURM_JOB_ID` / `PBS_JOBID` (falling back to a random UUID). The MATLAB default is a documented cause of slow pool startup and, on clusters where the home directory is distributed/network storage, race conditions or pool-start failures when concurrent jobs share it (MathWorks Answers #437544; multiple HPC-center deployment guides). Falls straight back to the original `'Processes'`/`'local'` string-profile calls if the relocation attempt itself fails, so this can never leave pool creation worse off than before. Windows behaviour is completely unchanged.
- **Live sub-case progress heartbeat**: `run_simulation_HPC_Linux.m` accepts optional `cfg.progress_queue` / `cfg.progress_tag` and, when set by the batch driver, sends a throttled (≥ 15 s apart per repeat) progress message through the same `parallel.pool.DataQueue` used for completion messages, reusing the `progress_pct` / `elapsed_ode` already computed for the existing `[ODE] Progress/ETA` log line (no extra computation cost). This is the only way to see a still-running repeat's progress while `cfg.silent = true` suppresses all of its normal console echo. `dispatch_progress_msg` (renamed from `on_run_complete`, which it now calls internally) routes heartbeat vs. completion messages; heartbeats are console-only and are never written to the persistent summary file.
- **`batch_run_from_json.m` — FAILED-run `Outdir`/`Log hint` display**: the serial catch block, the parallel `on_progress_impl` real-time log, and the final per-run summary table now print the recovered output directory and, when present, a `log_hint` pointing at that run's partial `simulation.log`, for any failed run where `capture_partial_outdir` (see below) could recover a path.

### Fixed

- **`capture_partial_outdir` in `batch_repeat_brownian.m` (and the equivalent addition to `batch_run_from_json.m`) — FAILED-run `outdir`/`log_hint` were always empty**: MATLAB does not propagate a function's output arguments to the caller when that function throws an uncaught error, even if the value (here, `run_simulation`'s `output_dir`) was already assigned internally before the error — verified empirically. Since every failed repeat's `output_dir` return value was therefore always `''`, the `outdir`/`log_hint` fields documented in v7.9.5 never actually populated for failures. Fixed by reconstructing the expected path from the repeat's deterministic `case_name` and confirming it exists on disk (`isfolder`) before reporting it, instead of relying on the unusable return value.
- **`batch_repeat_brownian_HPC_Linux.m` — invalid `char(datetimeObj, 'Format', fmt)` call**: `char()` on a `datetime` object only accepts a positional format string, not a `'Format'` name/value pair (that pair is only valid on the `datetime()` constructor itself). The previous call made MATLAB attempt to parse the literal text `'Format'` as the date format pattern, failing on the unsupported token `'F'`. Reproduced on both MATLAB R2023a (the target cluster) and R2026a — confirmed as a latent, version-independent bug rather than a platform/version compatibility issue — and fixed by switching to the positional-format signature. `finished_human`'s existing `char(datetime('now', 'Format', fmt))` call was already correct (the name/value pair there targets the `datetime()` constructor, not `char()`) and is unchanged.

### Changed

- **Outer `try`/`catch` around the whole `batch_repeat_brownian.m` / `batch_repeat_brownian_HPC_Linux.m` execution section**: a fatal, uncaught escape from the `parfor`/serial dispatch (e.g. a worker process crash not caught by the existing per-iteration `try`/`catch`) now falls through to the final report instead of aborting the script with no summary at all. Verified empirically that MATLAB's `parfor` commits sliced-output writes incrementally as each iteration returns rather than deferring all of them to loop completion, so iterations that finished before the abort are not lost. New `batch_aborted` / `batch_abort_identifier` / `batch_abort_message` footer fields record whether this path was taken.
- **`mem_per_worker_gb` documentation corrected (`batch_repeat_brownian.m`)**: the comment previously described a stale default (~1.5 GB) that no longer matched the actual in-effect value (2.5 GB).

### Compatibility

- **No change to any physical/numerical simulation output.** Every change in this release is confined to batch orchestration, platform-specific resource probing, parallel-pool management, and console/log display; the underlying physics, discretisation, and numerical results are bit-for-bit unaffected on both the original scripts and the new HPC_Linux fork.
- **`run_simulation.m` and `batch_repeat_brownian.m` (Windows-oriented originals) are unmodified except for the two reproducibility bug fixes above**, which apply identically regardless of platform. `run_simulation_HPC_Linux.m` / `batch_repeat_brownian_HPC_Linux.m` are new, independent files; using them is opt-in and has no effect on any existing workflow that continues to call the original scripts.
- **`batch_repeat_brownian_HPC_Linux.m` is not a drop-in replacement**: `generate_animation`/`compute_flowfield` default to **off** for all but the sampled repeats (first/middle/last), a deliberate memory-safety default for cluster batches that differs from the original script's "every repeat" behaviour. Set `heavy_output_sample_idx = 1:n_repeats` to restore full output for every repeat once memory headroom is confirmed safe for your configuration.
- **`mem_per_worker_gb`'s default (2.5 GB) is not empirically calibrated** for any specific configuration; the file header now says so explicitly and recommends measuring real peak RSS for your exact case before trusting the auto-sized pool on a real cluster allocation.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`run_simulation_HPC_Linux.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.10.0**, date **2026-07-03**.

---

## [7.9.6] — 2026-07-03

Console-output and error-resilience correctness release. A focused review of `run_simulation.m`'s logging/display pipeline uncovered a severe polarity bug that silenced ~100 `logf` console lines on every normal (non-batch) interactive run since v7.9.3, plus a `trisurf` axes-targeting crash that made `geometry_mode='custom'` runs fail deterministically at the summary-plot step. This release fixes both, adds a run-wide `try/catch` safety net so `simulation.log` always records the failure reason instead of truncating silently, and closes several smaller gaps (dual-channel warning omissions, `cfg.silent` coverage gaps, log encoding/timestamp consistency, Brownian log-flood throttling).

### Fixed

- **`run_simulation.m` — `dual_log` console-echo polarity inverted (critical)**: the `logf` closure forwarded `silent` un-negated to `dual_log`'s trailing `echo_console` flag. Since v7.9.3 this meant `cfg.silent=false` (the default, all interactive/non-batch runs) echoed **nothing** to the console through `logf` (~100 call sites: every `[INIT]`/`[PARAMS]`/`[ODE]`/`[POST]`/`[SAVE]`/`[PLOT]`/`[ANIM]`/`[FLOW]` line), while `cfg.silent=true` (batch/parallel workers) echoed everything — the exact opposite of the documented contract. `simulation.log` itself was never affected (the file-write branch is unconditional). Fixed by forwarding `~silent`. Verified with a live console-output smoke test.
- **`run_simulation.m` — custom-geometry `trisurf` axes-targeting crash**: both custom-body `trisurf(ax, F, x, y, z, ...)` calls (static summary-plot snapshot loop and per-frame animation loop) crashed with `MATLAB:catenate:dimensionMismatch` on any `geometry_mode='custom'` run, because MATLAB's `trisurf(tri, varargin)` binds its first positional argument unconditionally to the face matrix and never checks it for an axes handle (unlike `surf`/`plot3`, whose `axescheck`-based dispatch does support this). This made the `[PLOT] Generating summary figure` step — which cannot be disabled — fail deterministically for every custom-geometry run. Fixed by targeting the axes via the `'Parent'` name-value pair instead. Root-caused with an isolated MATLAB `trisurf` reproduction and verified end-to-end with a custom sphere body + pilus, including `generate_animation = true`.

### Added

- **`run_simulation.m` — run-wide `try/catch` with `[FATAL]` diagnostics**: the full run body (from logger setup to the closing `fclose`) is now wrapped in `try/catch`. Any unhandled error is recorded in `simulation.log` as a boxed `[FATAL]` block (identifier, message, full `getReport` stack trace) before being `rethrow`n, so a crash no longer leaves `simulation.log` truncated with no indication of the failure. The `catch` block also makes a best-effort attempt to close any `fig_preview` / `fig_summary` / `fig_anim` / `VideoWriter` handles left open by the failure.
- **`run_simulation.m` — `cleanup_log_fid` `onCleanup` guard**: guarantees `simulation.log`'s file handle is released on every exit path, including Ctrl-C interruption (which `try/catch` cannot intercept). Idempotent alongside the explicit `fclose` calls on the normal and `catch` paths.
- **`run_simulation.m` — `silent` + `preview_initial` deadlock guard**: `cfg.silent = true` with `preview_initial = true` and `auto_confirm = false` now raises `run_simulation:SilentPreviewDeadlock` immediately instead of blocking forever on an unanswerable `questdlg` on a headless/parallel worker. `fig_preview` visibility also now follows `silent` as a second line of defence.
- **`src/utils/project_version.m`**: new single-source-of-truth function for the runtime version string, referenced by `run_simulation.m`'s log banner.

### Changed

- **`src/visualization/compute_plot_flowfield.m` — new optional `silent` argument**: the chunked XZ/YZ flow-field evaluation progress bars (`render_progress_bar`) are now suppressed under `cfg.silent`, matching the ODE and animation progress bars. Previously these were the only progress bars not gated by `silent`, so batch/parallel workers with `compute_flowfield = true` still produced interleaved console output.
- **Dual-channel (`logf_h` + `warning`) consistency for 9 previously `warning`-only diagnostics**: `init_coupled_params.m` (`ConflictingOrientation`, `CustomAlpha0Assumption`, `DuplicateAttachSite`), `custom_body_building.m` (`DuplicatePoints`, `VolumeLengthFallback`, `TriangulationFailed`, `NonClosedSurface`, `UnexpectedFaceCount` — the latter two required adding a `logf_h` parameter to the local `check_closed_surface` helper), and `custom_pili_building.m` (`DeflectionAxisSwap`, `DuplicateCenterlinePoints` — required adding a `logf_h` parameter to the local `remove_consecutive_duplicates` helper). Previously these warnings were silenced with zero trace under `cfg.silent = true` (via `warning('off','all')`) since they never wrote to `logf_h`; they now always reach `simulation.log` regardless of `silent`.
- **`run_simulation.m` — `simulation.log` / `brownian_rng_seed.txt` opened with explicit `'n','UTF-8'` encoding**: matches the encoding already used by `save_sim_config.m` and the batch summary logs; removes a latent locale-dependent mojibake risk on non-English-locale Windows.
- **`run_simulation.m` — log timestamps now include the date** (`uuuu-MM-dd HH:mm:ss`, previously `HH:mm:ss`): avoids ambiguous/backwards-looking timestamps for multi-hour or overnight runs.
- **`run_simulation.m` — Brownian segment-log throttling**: `[ODE] Segment N: integrating...` / `[ODE] Segment N completed...` announcements are now only logged for the first 3 segments and every 50th segment thereafter when `enable_brown = true` (unthrottled for non-Brownian runs, which have few segments). A typical long Brownian run splits into thousands of short segments (one per `dt_B` node); without throttling this flooded `simulation.log` with tens of thousands of near-duplicate lines. Per-segment timing statistics (`seg_walltimes`, `seg_dtstar`) are unaffected — only the announcement text is throttled.
- **`src/config/validate_simulation_params.m` — aggregated validation errors now logged before throwing**: `logf_h('[VALIDATE] ERROR: ...')` is called immediately before the `error(...)` that reports aggregated parameter-validation failures, so the failure reason reaches `simulation.log` even for direct callers that do not wrap the call in `try/catch`.

### Compatibility

- **No change to any physical/numerical simulation output.** Every fix in this release is confined to console/log display, error handling, and a graphics-API workaround for custom-geometry rendering; the underlying physics, discretisation, and numerical results are bit-for-bit unaffected. `simulation.log` content for successful runs is unchanged except for the added date field in timestamps and (for Brownian runs) fewer `[ODE] Segment` announcement lines.
- `compute_plot_flowfield.m`'s new `silent` parameter is optional (defaults to `false`) and backward-compatible with any existing direct callers.
- Custom-geometry runs (`geometry_mode = 'custom'`) that previously crashed at `[PLOT] Generating summary figure` will now complete normally; this is a bug fix, not a behavioural change for previously-working configurations.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.9.6**, date **2026-07-03**.

---

## [7.9.5] — 2026-07-03

Silent-mode completeness for custom geometry and batch-summary reproducibility release. v7.9.4 closed the four main console-output channels under `cfg.silent` but custom-geometry builders still used bare `fprintf`, and the summary figure remained `Visible='on'`. Separately, a focused review of `batch_sweep/batch_repeat_brownian.m` identified gaps in summary metadata (no config snapshot, no pending/failed roll-up, inconsistent parallel empty-result messages, and lost partial output paths on mid-run crashes). **No simulation outputs change**; only console verbosity under `silent` and batch audit-trail completeness improve.

### Added

- **`src/geometry/custom_body_building.m`, `custom_pili_building.m` — optional `logf_h` argument**: when provided (via `init_coupled_params.m`), `[Custom Body]` / `[Custom Pili]` diagnostics route through the same `logf` closure as the rest of `run_simulation`, so they respect `cfg.silent` (file-only under silent; console echo when not silent).
- **`batch_sweep/batch_repeat_brownian.m` — reproducibility header fields**: `batch_id`, human-readable `started`, `config_source_type` / `config_json` or `config_source_note`, per-batch `batch_config_snapshot_<batch_id>.json` (replaces a single overwrite-prone filename), `project_root`, `script_path`, `brown_rng_policy`, `parallel_requested`, `summary_record_order`, and `seeds_precomputed_unused` when parallel seeds were pre-computed but serial execution ran instead.
- **`batch_sweep/batch_repeat_brownian.m` — footer roll-up fields**: `finished`, `failed_count`, sorted `failed_indices`, per-failure `failed_runs` block (`k`, `seed`, `id`, `msg`, `outdir`), `pending_count`, and `pending_indices` for runs still marked `PENDING` after an abnormal exit.
- **`batch_sweep/batch_repeat_brownian.m` — per-run FAILED audit lines**: on failure, `outdir` and `log_hint` (path to `simulation.log` when present) are written even when the run did not complete cleanly; OK lines use `path_for_log` for stable relative paths.

### Changed

- **`run_simulation.m` — summary figure hidden under `silent`**: `fig_summary` is created with `Visible='off'` when `cfg.silent=true`, matching the animation figure behaviour and preventing spurious figure windows during parallel batch workers.
- **`run_simulation.m` — header comment documents `cfg.silent`**: clarifies that it is a runtime-only batch flag (not stored in JSON), lists the four suppressed console channels, and notes that `simulation.log` is always written and errors are never suppressed.
- **`batch_sweep/batch_repeat_brownian.m` — serial fallback still writes resource probe**: when parallel is requested but falls back to serial, the summary header still records `[POOL] Resource probe` fields, `brown_rng_policy`, `parallel_requested`, and `seeds_precomputed_unused` so the log reflects what actually ran.
- **`batch_sweep/batch_repeat_brownian.m` — parallel empty-results message unified**: the post-`parfor` guard now uses the same `run_errmsg` text as the serial branch instead of a divergent placeholder string.
- **`batch_sweep/batch_repeat_brownian.m` — partial output path capture on error**: the `catch` block calls `capture_partial_outdir` so a FAILED record can still cite the run folder when `run_simulation` created one before aborting.

### Fixed

- **Custom-geometry console leak under `silent` (v7.9.4 compatibility gap)**: parametric and custom batch runs are now equally silent on workers; the v7.9.4 note that custom geometry still `fprintf`s to the console is resolved.

### Compatibility

- **Bitwise-identical to v7.9.4 for all simulation outputs.** `cfg.silent` defaults to `false`; JSON configs are unchanged. Summary-file per-run record format is backward-compatible (new header/footer keys are additive). The snapshot filename now includes `<batch_id>` — update any external scripts that assumed a fixed `batch_config_snapshot.json` name.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.9.5**, date **2026-07-03**.

---

## [7.9.4] — 2026-07-02

Layered-silence completeness release. v7.9.3 introduced `cfg.silent` but the silence gate covered only `dual_log` (~100 `logf` call sites), leaving five direct `fprintf` calls in `run_simulation.m` (notably `Log saved:` on every run and the `[OUTPUT]` rename/clear messages) plus **all** `src/` sub-function `warning` calls (~24 distinct IDs across 16 files) and two in-place `\b`-backspace progress bars (`odetpbar`, `render_progress_bar`) un-gated. In parallel mode these still produced thousands of interleaved console lines and backspace-character garbage. This release closes all four console-output channels via a layered design that keeps `src/` sub-functions untouched and restores caller `warning` settings via `onCleanup`. **No simulation outputs change**; the only behavioural difference is console verbosity in `silent` mode.

### Fixed

- **`run_simulation.m` — direct `fprintf` calls now gated by `silent`**: the six `fprintf` calls bypassing `dual_log` (four `[OUTPUT]` directory-conflict messages in the auto-confirm/rename/overwrite/abort branches, the preview-abort message, and the final `Log saved:` line) are each wrapped in `if ~silent`. The `[OUTPUT]` and preview messages fire before `log_fid` is open (or after it is closed), so they cannot be rerouted through `logf`; `if ~silent` is the correct gate. `Log saved:` is the most impactful fix — it was emitted on every run and dominated parallel stdout when N=10 workers.
- **`run_simulation.m` — `odetpbar` progress bar no longer mounted under `silent`**: `progressFunc` is set to `''` when `silent=true`, so `odeset('OutputFcn', '')` disables the ODE progress bar entirely. This eliminates the `\b`-backspace in-place overwrite that produced unreadable interleaved garbage across N workers.
- **`run_simulation.m` — `render_progress_bar` no longer called under `silent`**: the animation rendering progress bar is wrapped in `if ~silent`. Same rationale as `odetpbar`.
- **`run_simulation.m` — sub-function `warning` calls suppressed via `warning('off','all')` + `onCleanup`**: when `silent=true`, all `warning` display is toggled off at function entry and restored via an `onCleanup` object that fires on any exit path (normal return, error, or Ctrl-C). Enumerating the ~24 warning IDs across 16 `src/` files would be fragile and require touching every sub-function; `warning('off','all')` is the MATLAB-recommended temporary-suppress idiom and the `onCleanup` guard guarantees the caller's settings are never lost. Errors are NOT suppressed (only `warning`-level messages).

### Compatibility

- **Bitwise-identical to v7.9.3 for all simulation outputs.** `cfg.silent` defaults to `false`, so all non-batch (interactive) callers see identical console + warning behaviour as before. The `onCleanup` restoration ensures no caller-side `warning` setting is permanently altered even if `run_simulation` errors out under `silent=true`.
- Custom-geometry batch runs (`geometry_mode='custom'`) still emit `custom_body_building` / `custom_pili_building` `fprintf` calls to the console under `silent=true` (these sub-functions are not yet `silent`-aware); parametric batch runs are fully silent. A future release may thread `silent` into the custom-geometry builders if needed.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.9.4**, date **2026-07-02**.

---

## [7.9.3] — 2026-07-02

Batch-logging robustness and parallel stdout-hygiene release. Two issues identified during a focused review of `batch_sweep/batch_repeat_brownian.m` are addressed: (1) the per-run summary records were written in a single batch at the end of the loop, so a mid-batch crash lost **all** per-run status/timing/error information (only the header survived on disk); (2) in parallel mode each worker echoed its full `run_simulation` log (~100–250 lines per run, dominated by `[ODE] Segment N completed` lines) to the shared stdout, producing ~15000 interleaved and unreadable lines for a 100-run × 10-worker batch and burying the useful `[k/100] OK` progress lines. **No simulation outputs change**; this release only strengthens batch logging and console readability.

### Added

- **`run_simulation.m` — `cfg.silent` field (default `false`)**: when `true`, the `dual_log` helper suppresses the console `fprintf` echo while still writing every line to the run's `simulation.log` file. The `logf` closure now forwards `silent` as a trailing flag to `dual_log`, which strips it (detected via `islogical`) and gates the console write. Backward-compatible: absent or `false` preserves the legacy dual-channel (console + file) behaviour; `logf_file` (file-only variant) is unchanged.
- **`batch_sweep/batch_repeat_brownian.m` — incremental per-run summary writes**: each run's record (`NNN seed=... status=OK/FAILED wall=...` plus, on failure, `Identifier` / `Message` / `Location`, or `out:` on success) is now appended to `batch_repeat_summary_*.txt` immediately after the run completes and the file is `fclose`-d (forced flush). Serial branch: `append_run_summary` is called at the end of each loop iteration. Parallel branch: the `afterEach` DataQueue callback (`on_run_complete`) calls `append_run_summary` on the main thread (no concurrent-write race). The previous end-of-batch aggregate write loop is removed; the footer (`elapsed_total_s`, `success_count`, parallel speedup) is now written in append mode. **Result**: a mid-batch crash now leaves the summary file with the header plus all completed-run records on disk; only the footer is lost.
- **`batch_sweep/batch_repeat_brownian.m` — `cfg_k.silent = true` in the parfor branch**: workers no longer echo their per-run logs to the shared stdout. Each run's full detail remains in `output/<case>/<NNN>/simulation.log` for post-hoc inspection. Serial branch is unchanged (verbose) so interactive single-run debugging is unaffected.
- **`batch_sweep/batch_repeat_brownian.m` — `append_run_summary` and `on_run_complete` helper functions**: the former opens the summary in append mode, writes one record, and `fclose`s (no-op when `sum_path` is empty, `warning` on open failure); the latter is the DataQueue callback wrapping `print_progress` + `append_run_summary` and defensively extracts `seed` / `outdir` / `stack_name` / `stack_line` via `isfield`.

### Changed

- **`batch_sweep/batch_repeat_brownian.m` — summary file open-failure now warns**: previously `fid_sum == -1` silently set `sum_path = ''`; a `warning('batch_repeat_brownian:SummaryOpenFailed', ...)` is now emitted so the user knows per-run records are not being persisted.
- **`batch_sweep/batch_repeat_brownian.m` — removed duplicate `[POOL] Resource probe` console block**: `decide_pool_size` previously printed an 8-line `[POOL] Resource probe:` block (physical cores, memory, caps, binding, BLAS threads) that duplicated information already shown in the main batch header. The block is removed; `info.core_probe` (the probe-source label, e.g. `wmic` / `CIM`) is still recorded in the summary file's `core_probe:` field, so no diagnostic metadata is lost. The `blas_threads` parameter of `decide_pool_size` is retained for API stability (annotated `%#ok<NASGU>`; not currently read inside the function).

### Compatibility

- **Bitwise-identical to v7.9.2 for all simulation outputs.** The `cfg.silent` field is optional with default `false`; `dual_log` defaults to `echo_console = true` when the trailing flag is absent (preserves any non-`logf` callers). The summary-file format (header, per-run records, footer) is byte-for-byte identical to v7.9.2 when the batch completes normally; only the write timing changes (incremental vs end-of-batch). Parallel-mode per-run records are written in completion order (the leading `NNN` index preserves `k` for sorting), while serial-mode records remain in `k` order.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.9.3**, date **2026-07-02**.

---

## [7.9.2] — 2026-06-30

Numerical-robustness and documentation-accuracy release. A front-to-back critical re-review of the underlying mathematical physics surfaced four actionable items: (1) a defensive quaternion normalisation guard against `ode15s` internal trial-step drift, (2) a corrected Euler-angle convention label (the directors follow Walker (2020)'s `R = Rz(phi)·Ry(-theta)·Rz(psi)`, NOT textbook ZXZ as previously documented), (3) a documented declaration of the deterministic/stochastic operator-splitting approximation in the Brownian path, with a quantitative `rel_asym` diagnostic, and (4) removal of a redundant symmetrisation in the Brownian mobility dimensionalisation. **No simulation outputs change** for configurations with a unit quaternion and reciprocal hydrodynamic discretisation; the quaternion guard only activates on numerical drift, and the symmetrisation removal is algebraically exact (Hadamard division by a symmetric outer product preserves symmetry).

### Fixed

- **`src/geometry/reconstruct_geometry.m` — defensive unit-quaternion normalisation**: `ode15s` internal trial steps may produce states with `|q| != 1`; the rotation-matrix formula in `quat_to_rotmat` is only orthonormal for `|q| = 1`, so a non-unit quaternion would silently deform the body (non-rigid surface points, non-conserved volume). The body quaternion is now normalised immediately after `unpack_state`, with a zero-norm fallback to identity orientation and a diagnostic warning. Only the local copy is touched; the solver's state vector and error control are not disturbed.
- **`src/math/quat_to_rotmat.m` — defensive normalisation at the leaf**: a second normalisation guard at the function entry protects any other callers (test scripts, post-processing, visualisation helpers) against `|q| != 1`. Cost: one 4-vector norm + division per call — negligible.

### Changed

- **`src/physics/compute_effective_body_brownian.m` — remove redundant symmetrisation of `M_eff_dim`**: the dimensional mobility `M_eff_dim = M_eff_dimless ./ (S * S')` is already symmetric because `M_eff_dimless` was symmetrised earlier (Step 3) and Hadamard division by the symmetric outer product `S * S'` preserves symmetry. The explicit `0.5*(M + M')` line is removed and replaced with a derivation comment. Algebraically exact; no numerical change.
- **`src/physics/compute_effective_body_brownian.m` — quantitative asymmetry diagnostic**: in addition to the existing `rel_asym > 1e-3` warning, a new `rel_asym > 1e-2` log line reports the current `eps_body` and `eps_pili` values and points the user to `theory.md §10 "Known approximation"` and §5.1 (RSM leak tolerance). This makes the deterministic/stochastic operator-splitting bias visible and actionable without changing the numerical path.
- **Euler-angle convention label corrected (comment-only)**: the original `pili_euler_angles.m` docstring labelled the rotation sequence "ZXZ intrinsic rotation (`R_z(phi) → R_x(theta) → R_z(psi)`)". Re-analysis against the actual `pili_directors.m` formulas shows the convention is Walker (2020)'s `R = Rz(phi)·Ry(-theta)·Rz(psi)` (the transpose of the standard ZYZ active rotation), NOT textbook ZXZ. The `d3` director `d3 = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))` is the spherical-coordinate tangent; the clamped-BC Jacobian `J_euler_basis` in `build_coupled_lin_sys.m` is derived consistently from this convention. Affected files: `pili_euler_angles.m`, `pili_directors.m`, `build_coupled_lin_sys.m` (new derivation comment on `J_euler_basis`), `init_coupled_state.m` (ZYZ → Walker wording), `docs/index.html`, `CHANGELOG.md`. Comment-only correction; bitwise identical outputs.

### Documentation

- **`docs/theory.md`, `docs/theory.zh-CN.md` — new "Known approximation" subsection in §10 (Brownian Dynamics)**: declares that the deterministic solve uses the raw (structurally non-symmetric when `eps_body != eps_pili`) hydrodynamic mobility `M_weighted`, while the Brownian covariance path uses the symmetrised effective body mobility `M_eff = 0.5*(M + M')` to satisfy the FDT requirement of a symmetric diffusion tensor. Documents the `rel_asym` diagnostic, its typical magnitude (`~1e-5..1e-3` for mixed-epsilon configurations, `<1e-10` for pure RSM), the Schur-complement attenuation, and user guidance for `rel_asym > 1e-2`.
- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.9.2**, date **2026-06-30**.

### Compatibility

- **NOT bitwise-identical to v7.9.1** only in pathological cases where `ode15s` produced a non-unit quaternion that was previously passed through to `quat_to_rotmat` unnormalised. For all well-behaved runs the quaternion guard is a no-op (the quaternion is already unit to machine precision after the RHS projection in `ode_coupled_rhs.m`), and the symmetrisation removal is algebraically exact, so outputs are identical to v7.9.1.
- The Euler-angle convention correction and Brownian documentation are comment/doc-only; no numerical change.

---

## [7.9.1] — 2026-06-22

Post-review bug fixes — addresses two issues found in a front-to-back data-flow review: a non-physical negative diffusion coefficient in the batch post-processor, and missing error stack traces in three batch sweep scripts. **No simulation outputs change**; this release only fixes post-processing and batch error reporting.

### Fixed

- **`post_processing/postprocess_cases.m` — clamp `D_trans` to non-negative (Bug 2)**: the MSD least-squares fit `MSD(τ) = 6Dτ + V²τ²` could produce a negative `coeff(1)` (and hence a negative `D_trans_um2_s`) when the true motion is nearly ballistic (`D → 0`) or when finite-sampling noise dominates. `D_trans_um2_s` is now clamped with `max(0, coeff(1)) / 6`, mirroring the existing `max(0, coeff(2))` guard on `V_msd_um_s`. This eliminates non-physical negative diffusion coefficients in the Excel summary report.
- **`batch_sweep/batch_sweep_alpha.m`, `batch_sweep/batch_sweep_pili_length.m`, `batch_sweep/batch_repeat_brownian.m` — capture error stack traces (Bug 6)**: the three batch scripts previously stored only `ME.message` on failure, losing the `ME.identifier` and `ME.stack(1)` location. All catch blocks now capture and write the error identifier, message, and top-of-stack location (`file (line N)`) to both the console and the summary log file, matching the existing style of `batch_sweep/batch_run_from_json.m`. In `batch_repeat_brownian.m`, the parallel branch additionally forwards `identifier` / `stack_name` / `stack_line` through the `DataQueue` so the `print_progress` callback can display them on the console in real time; the "empty results" failure path now sets a `"N/A (empty results)"` identifier placeholder instead of leaving a blank `Identifier :` line.

### Compatibility

- **Bitwise-identical to v7.9.0** for all simulation outputs. The MSD clamp only affects the post-processing Excel report (and only when the unclamped fit was already negative); the batch error-reporting changes only affect failure-path logging.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.9.2**, date **2026-06-30**.

---

## [7.9.0] — 2026-06-17

Parameter-validation deep refactor — consolidates all pre-flight parameter checks into `validate_simulation_params.m`, fixes a custom-mode `alpha_0` regression caught by static review, and removes dead code in `init_coupled_params.m`. **No numerical change to simulation outputs**; this release strengthens input validation and error reporting only.

### Changed

- **`run_simulation.m` — pass the full `cfg` struct to the validator**: the manual partial-`val_cfg` construction is replaced by `validate_simulation_params(cfg, logf)`. All field accesses inside the validator are guarded by `isfield()`, so adding a new user-configurable parameter no longer requires updating the call site (A2).
- **`src/config/validate_simulation_params.m` — full-cfg acceptance (A2)**: every field access is now guarded by `isfield()` / type checks, so the validator accepts the complete `cfg` struct from `run_simulation` / `load_sim_config` / batch sweeps without requiring a pre-curated subset. Missing optional fields are silently skipped.

### Added

- **`src/config/validate_simulation_params.m` — `pili_ext` count consistency (P1)**: in both parametric and custom geometry modes, an explicit error is raised when `numel(cfg.pili_ext)` does not match the number of defined pili (`size(Pili_Matrix, 2)` or `numel(custom_geometry.pili)`), preventing the silent fallback in `init_coupled_params` from masking user config errors.
- **`src/config/validate_simulation_params.m` — `delta` range check (P2)**: `delta <= 0` now errors and `delta >= pi/4` warns, since extreme values break the Euler-angle singularity event in `run_simulation`.
- **`src/config/validate_simulation_params.m` — `N_body` upper-bound warning (P4)**: `N_body > 2000` warns about the `O(N^2)` distance matrix in `init_coupled_params` with an estimated memory footprint.
- **`src/config/validate_simulation_params.m` — background-flow orthogonality warning (P5)**: for `shear_z` flow, `direction` nearly parallel to `shear_axis` (`|cos(theta)| > 0.99`) warns that vorticity is ~0 and pili experience no shear-driven rotation despite a non-zero translational flow.
- **`src/config/validate_simulation_params.m` — `alpha_0` range check (P6/A1)**: the `[-pi/2, pi/2]` range check for `alpha_0` was moved here from `init_coupled_params.m` so all parameter errors are reported in one pass. Enforced in parametric mode and in custom mode with empty `initial_orientation`; skipped in custom mode with non-empty `initial_orientation` (where `alpha_0` is ignored by `init_coupled_params`).
- **`src/config/validate_simulation_params.m` — unified string/struct acceptance (P7)**: `intrinsic_curvature`, `active_moment`, and `motor_torque` now consistently accept bare strings (e.g. `'off'`, `'frozen'`) by coercing them to `struct('type', char(...))`, matching the `parse_pili_config` / `init_coupled_params` convention. Unknown string types now produce a type error instead of being silently dropped.
- **`src/config/init_coupled_params.m` — `motor_torque` string acceptance (P8)**: both parametric and custom branches now coerce string `motor_torque` inputs to `struct('type', char(...))`, mirroring `intrinsic_curvature` and `active_moment`.

### Fixed

- **`src/config/validate_simulation_params.m` — `alpha_0` range check lost in custom mode (Bugbot)**: the P6/A1 refactor initially guarded the `alpha_0` check with `strcmp(geometry_mode, 'parametric')`, which let out-of-range `alpha_0` pass silently in custom mode when `initial_orientation` is empty (the fallback path where `init_coupled_params` still consumes `alpha_0` to build `q_0_init`). The guard now also fires in custom mode with empty `initial_orientation`.
- **`src/config/validate_simulation_params.m` — `Temper` type robustness**: the Stokes-Einstein diagnostic now checks `isnumeric(cfg.Temper)` before `isfinite`, guarding against non-numeric `Temper` values from misconfigured JSON.
- **`src/config/validate_simulation_params.m` — `AbsTol`/`RelTol` warning wording (P9)**: the warnings now explicitly state the recommended range `1e-6 ~ 1e-3`, so the default `1e-4` is not mistaken for "out of range".
- **`src/config/validate_simulation_params.m` — error aggregation performance (P10)**: the aggregated error message is now built with a single loop + `strjoin` instead of `cellfun` + `strjoin`, reducing intermediate cell allocation when the error list is long (e.g. per-pilus errors across many pili).
- **`src/config/init_coupled_params.m` — dead-code removal (L314)**: removed an `if/else` whose two branches both assigned `Pili_Matrix = zeros(8, 0)`; the custom-mode pilus count is recomputed later from `numel(custom_pili_cfg)`, so the branch was a no-op. Eliminates a MATLAB lint ERROR.

### Compatibility

- **Bitwise-identical to v7.8.0** for all simulation outputs. This release only tightens pre-flight validation and fixes a latent custom-mode `alpha_0` acceptance bug; no ODE rhs, linear-system, or post-processing code changed.
- **Stricter input rejection**: configurations that previously passed validation with mismatched `pili_ext` counts (P1), out-of-range `delta` (P2), or out-of-range `alpha_0` in custom mode (Bugbot fix) will now error. These were always misconfigurations; they now fail fast with a clear message instead of silently producing unintended geometry.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.9.0**, date **2026-06-17**.

---

## [7.8.0] — 2026-06-11

Reference-curvature and motor-torque nondimensionalisation correctness fixes — **changes numerical outputs for curved/helical `frozen`-curvature pili and for `free_torsion` constant-motor configurations**.

### Fixed

- **`src/config/init_coupled_params.m` — intrinsic curvature `kb0` discretisation consistency**: the reference curvature `kb0` is now built with the same node-director operator (`compute_node_directors`, SVD-projected segment averaging) and the same finite-difference stencil used by the run-time curvature in `src/physics/compute_pili_elastic_moments.m`, replacing the former "duplicate-end" nodal-director construction (`d1NC0 = [d1_0; d1_0(end,:)]`). A `frozen`-curvature pilus in its reference configuration is now an exact zero-moment equilibrium (`kappa_live - kb0 == 0` to machine precision), removing the spurious initial elastic restoring moment that previously acted on curved or helical pili at `t = 0`. Applied to both the parametric and custom geometry branches.
- **`src/config/init_coupled_params.m` — constant motor-torque nondimensionalisation**: the motor torque enters the segment-1 `psi_1` twist moment-balance row, whose moment scale is `EI / L_pili` under the per-pilus normalised arc length `s' in [0,1]` (`ds = 1/N_k`). The nondimensionalisation now divides by `EI / L_pili` (with `L_pili = ell * L_0`) instead of `EI / L_0`, matching the active-moment convention. This restores the intended applied motor-torque magnitude, which was previously off by the factor `ell = L_pili / L_0`. The stall-speed nondimensionalisation (`max_omega_dimless`) was already correct and is unchanged. Applied to both the parametric and custom geometry branches.

### Compatibility

- **Not bitwise-identical to v7.7.0** for: (a) pili with curved or helical `frozen` intrinsic curvature — straight pili are unaffected, since their reference curvature is zero under both schemes; (b) `free_torsion` pili driven by a `constant` motor torque, where the effective torque changes by the factor `ell`. All other configurations (pili-free, straight-pili, `zero`/`constant`/`helix` intrinsic-curvature overrides, and active-moment-only driving) are numerically unchanged.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.8.0**, date **2026-06-11**.

---

## [7.7.0] — 2026-06-10

Strict Brownian Ito drift and physics-correctness cleanup — **changes numerical outputs for Brownian, pili-bearing torque-balance, and time-dependent active-moment configurations**.

### Changed

- **`run_simulation.m` — strict body-only Ito Brownian step**: Brownian kicks now use the explicit drift term `kBT*div(M_eff)` computed by centered finite differences in the six body coordinates, plus correlated noise from the current SPD-projected diffusion tensor.
- **`src/physics/compute_effective_body_brownian.m` — SPD covariance projection**: diffusion tensors are projected through symmetric eigenvalue clamping before Cholesky sampling, replacing the previous unit-diagonal regularisation fallback.
- **Brownian configuration cleanup**: removed the former user-facing midpoint drift option from the MATLAB config path, JSON load/save path, validators, example JSON files, batch scripts, website, and user documentation. Brownian dynamics now expose one stochastic interpretation: the strict Ito path.

### Fixed

- **`src/core/build_coupled_lin_sys.m` — Kirchhoff torque-density scaling**: force-density columns keep the `ell_k^3` scaling while torque-density columns now use `ell_k^2`, matching the line-torque dimensional scaling.
- **`src/visualization/compute_plot_flowfield.m` and `run_simulation.m` — final flow-field phase**: final flow-field force extraction now calls `ode_coupled_rhs` at `T_achieved` instead of hard-coding `t=0`, so time-dependent active moments keep the correct final phase.
- **Brownian validation**: Brownian-enabled runs now require a finite positive temperature before diffusion and drift computation.

### Compatibility

- **Not bitwise-identical to v7.6.0** for Brownian-enabled cases, torque-active pili cases, or cases with time-dependent active moments shown in final flow-field panels.
- JSON files containing the removed midpoint Brownian field are still loadable because unknown JSON fields are ignored, but newly saved configs and examples no longer write that field.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.7.0**, date **2026-06-10**.

---

## [7.6.0] — 2026-06-10

Non-local pili twist-coupling physics upgrade — **changes numerical outputs for all pili-bearing configurations**.

### Changed

- **`src/core/build_global_resistance.m` — replaced local RFT axial twist resistance with non-local `M_twist` mobility**: pili twist now uses a segment-integral axial rotation mobility assembled from Cortez (2018) regularized-segment recurrences, with the scalar axial torque density embedded through `D_twist`.
- **`src/core/build_coupled_lin_sys.m` and `src/core/ode_coupled_rhs.m` — Schur-complement twist solve**: the coupled force/twist system now solves through `M_weighted`, `M_twist`, and the reciprocal twist cross block instead of the former local `A_twist_total` path.
- **Physical rotlet normalization**: the stored Cortez (2018) segment rotlet is used with a `0.5` factor in the cross-coupling block, because Cortez's rotlet kernel is the full curl of the regularized Stokeslet and is twice the standard physical rotlet.

### Fixed

- **`src/core/build_coupled_lin_sys.m` — line-integrated net pili torque now enters the system torque balance (rows 4-6)**: the global angular-momentum balance uses `sum(r x f) + sum(L_s m_s)`, not bare `sum(m)`, because the rotlet strength is a torque density along each segment. With this fix a force/torque-free swimmer remains torque-free, and a motor-driven `free_torsion` pilus produces the physically correct body counter-rotation about the spin axis.
- **`src/visualization/compute_plot_flowfield.m` — flow-field reconstruction follows the same non-local twist solve and physical `0.5*M_rotlet` scaling** used by the dynamics.

### Compatibility

- **NOT bitwise-identical to v7.4.6** for any configuration containing pili. The non-local twist solve, physical standard-rotlet scaling, and line-integrated `sum(L_s m_s)` torque term can change the solved body angular velocity and full trajectory. The effect is small for weak twist activity but can be O(1) for fast motor-driven or strongly spinning pili. **Pili-free configurations are unaffected.**

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**, **`docs/presentation.html`**: version **7.6.0**, date **2026-06-10**.

---

## [7.4.6] — 2026-06-03

Code-quality and documentation-consistency pass — **no numerical change** to `run_simulation` outputs.

### Changed

- **`src/math/pili_spatial_coords.m`**: replaced the per-segment cumulative `for` loop with a vectorised `cumsum` (mathematically identical node positions; output shape `[N_k+1 x 3]` unchanged).
- **`src/config/validate_simulation_params.m`**: corrected the rotational relaxation time used in the `dt_B` pre-flight heuristic from `1/(6·D_R)` to `1/(2·D_R)` (standard first-order orientational correlation time). This is a warning threshold only and does not affect time integration.

### Fixed

- **`src/geometry/pili_building.m`**: removed the dead `compute_curvature_binormal` call whose `kb`/`kb0` output was never read downstream (curvature is recomputed by `init_coupled_params` as `params.kb0{k}`, `[N+1 x 3]`). The struct fields are now `zeros(3,0)` placeholders, matching the established `custom_pili_building.m` pattern.
- **`run_simulation.m`**: replaced the invalid `%#ok<NFFT>` lint directive on the `logf` closure with a plain comment; added `%#ok<AGROW>` to the intentional in-loop array growth in the initial-pose preview, snapshot, and animation bounding-box passes.

### Documentation

- **`src/core/build_coupled_lin_sys.m`**, **`CHANGELOG.md`**: removed stale references to non-existent "regression tests"; condition-number health is now described as monitored via runtime diagnostics (`compute_effective_body_brownian.m`) and the `validate_simulation_params` pre-flight checks.
- **`main_coupled_simulation.m`**: corrected the "README defaults" wording to point at the actual source of the `T_obs` / `dt` / `rho_bac` fallback values (`load_sim_config.m`).
- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**: version **7.4.6**, date **2026-06-03**.

### Compatibility

- **Bitwise identical numerical outputs vs. v7.4.5** for all configs. Every change is a comment, a lint directive, an equivalent vectorisation, dead-code removal, or a warning-only heuristic constant — no solver, geometry, or physics path is altered.

---

## [7.4.5] — 2026-05-25

Parallel JSON batch runner and batch-log polish — **no numerical change** to single-run `run_simulation` outputs.

### Added

- **`batch_sweep/batch_run_from_json.m` — optional `parfor` parallel batch over a folder of JSON configs**:
  - User switches: `use_parallel` (default `true`), `mem_per_worker_gb` (default `3.0`), `max_resource_frac` (default `0.75`), `worker_blas_threads` (default `1`).
  - **Pool sizing:** same probe policy as `batch_repeat_brownian.m` — `n_pool = min(floor(n_phys × max_frac), floor(available_GB × max_frac / mem_per_worker_gb), n_run)`.
  - **Parallel Computing Toolbox** required for the parallel branch; missing license or `parpool` failure **falls back to serial**.
  - **Pre-scan deduplication:** non-empty duplicate `case_name` values (alphabetical first file wins) are marked `SKIPPED` and never queued; empty `case_name` entries are not deduplicated.
  - **`parallel.pool.DataQueue`** streams completion lines to the console and **`batch_summary_<timestamp>.log`** in the JSON folder (`fflush` after each event).
  - Worker init: `DefaultFigureVisible` off + `maxNumCompThreads` cap via `parfevalOnAll` with **`wait`**; pool lifetime comment documents `delete(gcp('nocreate'))` to restore defaults.
  - Final summary reports **Total JSON files**, **Queued to run**, **Succeeded** / **Failed** / **Skipped**, optional parallel speedup estimate, and per-file `V_sink` table.
  - Brownian RNG, animation, and flow-field flags follow each JSON (only `auto_confirm` and `preview_initial` are forced headless).

### Changed

- **`batch_sweep/batch_run_from_json.m` — logging consistency (P0+P2 polish):** parallel progress uses the shared `blog` closure; execution-mode header lines are mirrored to console and batch log; unified run header `Run queue kk/nn  [JSON ii/tt]: filename`; itemized batch-complete counts.

### Documentation

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**: version **7.4.5**, date **2026-05-25**; JSON parallel batch runner and PCT scope notes.

### Compatibility

- **Bitwise identical numerical outputs vs. v7.4.4** for identical single-run configs. Only the batch driver `batch_run_from_json.m` gains scheduling, deduplication, and logging code; `use_parallel = false` preserves serial JSON batch behaviour for the queued files.

---

## [7.4.4] — 2026-05-22

Configuration validation hardening — **no numerical change** to simulation outputs for configs that previously passed validation.

### Fixed

- **`validate_simulation_params.m` — parametric/custom no-mixing (reverse direction)**: when `geometry_mode = 'parametric'`, a non-empty `custom_geometry` block is now rejected with a clear error, symmetric to the existing custom-mode rejection of parametric fields. Aligns runtime behaviour with `docs/theory.md` / `docs/theory.zh-CN.md` §1.
- **`validate_simulation_params.m` — `helix` intrinsic curvature fields**: `intrinsic_curvature.type = 'helix'` now requires finite scalar `.kappa_d2` and `.kappa_d3` (including rejection of bare string `'helix'`), mirroring the existing `'constant'` / `.kappa` check and preventing obscure init-time field errors.
- **`run_simulation.m`**: `val_cfg` now forwards a non-empty `custom_geometry` in all geometry modes so the reverse no-mixing rule is enforced on the main entry path.

### Changed

- **`README.md`**, **`README.zh-CN.md`**, **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**, **`requirements.txt`**, **`run_simulation.m`**, **`main_coupled_simulation.m`**, **`docs/theory.md`**, **`docs/theory.zh-CN.md`**, **`docs/index.html`**: version **7.4.4**, date **2026-05-22**; validation hardening and geometry-mixing / `helix` field notes in user guide.

### Compatibility

- **Bitwise identical numerical outputs vs. v7.4.3** for identical inputs that previously passed validation. Only invalid or mixed configurations that previously slipped through (or failed late in init) are affected.

---

## [7.4.3] — 2026-05-22

Documentation and path alignment — **no numerical change** to simulation outputs vs. v7.4.2.

### Fixed

- **`docs/theory.md`**, **`docs/theory.zh-CN.md`**: corrected `src/` links for `reconstruct_geometry.m`, `compute_pili_elastic_moments.m`, and `compute_plot_flowfield.m`; added sinking-force, root-gap `dis`, custom-mode `L_0`, lab-frame `\Omega`, and `\varepsilon_k` notes; English §1 now states bidirectional custom/parametric mixing rejection.
- **`README.md`**, **`README.zh-CN.md`**: JSON examples now use `configs/` paths; corrected `configs/` git-tracking wording; documented `batch_sweep/` (drivers) vs `batch sweep/` (sweep output); Chinese core-files table uses `src/` paths.
- **`USER_GUIDE.md`**, **`USER_GUIDE.zh-CN.md`**: JSON schema by geometry mode, configuration enumerations, corrected paths and batch folder convention (Scheme A).
- **`docs/index.html`**: v7.4.2 release card content aligned with CHANGELOG; batch driver paths use `batch_sweep/`; JSON examples use `configs/`.
- **`requirements.txt`**: batch driver path and release version.
- **`main_coupled_simulation.m`**: JSON example comment uses `configs/config_custom_example.json`.
- **`build_coupled_lin_sys.m`**: removed stale reference to non-existent `tests/test_motor_reaction_torque.m` (comment only).

### Compatibility

- **Bitwise identical numerical outputs vs. v7.4.2** for identical inputs that previously passed validation. No executable simulation algorithm was modified.

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
- **`pili_euler_angles.m`**: added Euler angle convention specification. NOTE: the original entry labelled this "ZXZ intrinsic rotation (R_z(phi) → R_x(theta) → R_z(psi))"; re-analysis against the actual `pili_directors.m` formulas shows the convention is Walker (2020)'s R = Rz(phi)·Ry(−theta)·Rz(psi) (the transpose of standard ZYZ active rotation), NOT textbook ZXZ. The d3 director d3 = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)) is the spherical-coordinate tangent; the clamped-BC Jacobian `J_euler_basis` in `build_coupled_lin_sys.m` is derived consistently from this convention. Comment-only correction; bitwise identical outputs.

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
- **`validate_simulation_params.m`**: Added a top-level `geometry_mode` switch. For `'custom'` mode, calls `validate_custom_geometry.m` and includes explicit checks to error if any parametric-specific fields (`Pili_Matrix`, `minor_axis`, `major_axis`, `N_body`) are set. Common validations apply regardless of geometry mode. Brownian validation now correctly handles custom pili. The Stokes-Einstein `dt_B` warning applies only in parametric mode where `major_axis` is defined.
- **`run_simulation.m`**: Unpacks `geometry_mode` and `custom_geometry_cfg` from `cfg`. Parameter validation passes these new fields. Logging conditionally prints geometry details. Visualization routines (initial pose preview, summary plots, and animation) use `trisurf()` with the automatically generated mesh for custom bodies, replacing `surf()` for parametric capsules. Axis limit calculations adapt for both geometry types.
- **`main_coupled_simulation.m`**: Added comments explaining the `geometry_mode` options and a commented-out inline configuration example for custom mode. Inline parameter structure refined with an `if/else` block for parametric vs. custom-specific configuration, with common parameters defined outside the block. Explicit initialization of `cfg.alpha_0`, `cfg.Pili_Matrix`, and `cfg.pili_ext` in the custom branch ensures these top-level fields are always defined.

### Changed

- **"No mixing" enforcement**: When `geometry_mode = 'custom'`, any parametric-specific configuration fields (`Pili_Matrix`, `minor_axis`, `major_axis`, `N_body`) are explicitly rejected by `validate_simulation_params.m` with clear error messages.
- **`custom_pili_building.m`**: Removed redundant `kb`/`kb0` (curvature binormal) computation since `init_coupled_params` recomputes all curvature information from scratch. These fields are now initialized to `zeros(3,0)` as placeholders.
- **Custom body visualization**: The automatic surface triangulation (via `alphaShape` + `convhulln` fallback) generates `body_faces`/`body_verts` stored in `body_meta`, enabling `trisurf` rendering in `run_simulation.m` without requiring a separate mesh file.

### Fixed

- **Inline custom mode**: Missing `cfg.alpha_0`, `cfg.Pili_Matrix`, and `cfg.pili_ext` fields in the custom branch of `main_coupled_simulation.m` caused runtime "Reference to non-existent field" errors. All three fields are now explicitly initialized.
- **Custom body volume `NaN` propagation**: If `custom_body_building.m` failed to triangulate, `V_body = NaN` would silently propagate to physical parameters. `init_coupled_params.m` now errors explicitly on invalid volume.
- **Brownian warning in custom mode**: The custom-mode Brownian validation path now also considers `cfg.custom_geometry.pili`.
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

- **`Brownian_Implementation_Plan.md` (tracked companion)**: New repository document with **§1.5 — Scope & Simplifying Assumptions** (modelled body 6 DOF + pili-induced drag through the coupled mobility; deliberate omission of pili internal thermal modes; extension references) and **§2.4.1 — Stochastic Interpretation & Missing Drift**. The file was previously listed in `.gitignore` and therefore missing from a clean clone; it is now versioned so cross-references in `validate_simulation_params.m` resolve.
- **`validate_simulation_params.m`**: Two new Brownian-related warnings. One covered configuration-dependent Brownian drift in pili-bearing cases, and one warns when `dt_B` exceeds 10% of a Stokes-Einstein back-of-envelope estimate of `min(tau_trans, tau_rot)` (using the major semi-axis as a conservative radius). Both warnings are inert when `enable_brown = false`.
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
  - `'CheckCondition' = false` disables MATLAB's per-solve automatic `rcond` estimate (`O(N²)` saved per invocation); overall condition health is validated via runtime conditioning diagnostics.
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
- **Brownian drift option**: Added an early optional midpoint approximation for configuration-dependent diffusion drift. This option has since been removed in favor of the strict explicit Ito drift path.

### Changed

- **`ode_coupled_rhs.m`**: Refactored to call `build_coupled_lin_sys` for Steps 1-7, retaining only elastic moment computation and state-rate packing.
- **`run_simulation.m`**: Brownian kick block replaced: uses per-step dynamically computed `L_eff` (6×6 Cholesky) from full-system effective friction and generates 6D correlated noise (translation-rotation coupling included).
- **`init_coupled_params.m`**: Removed one-time body-only `compute_body_brownian_params` call at initialisation; Brownian diffusion is now computed dynamically per step.
- **`load_sim_config.m`**: Added Brownian JSON parsing for the full-system effective friction path.
- **`config_example.json`**: Added Brownian configuration fields to the `brownian` block.
- **`validate_simulation_params.m`**: Added validation for Brownian configuration parameters.
- **`Brownian_Implementation_Plan.md`**: Added upgrade notice documenting the transition from body-only to full-system effective friction.
- **`README.md`**: Updated section 16, file listings, architecture diagram, and feature descriptions.

### Removed

- **`compute_body_brownian_params.m`**: Deleted. Old body-only RSM friction matrix implementation fully superseded by `compute_effective_body_brownian.m`.
- **`test_effective_brownian.m`**: Deleted. Verification tests passed (7/7); script no longer needed.

### Fixed

- **`main_coupled_simulation.m`**: Added missing Brownian configuration wiring to the inline configuration block.
- **`configs/OnePili-P1-L3.json`**: Added missing Brownian fields to the `brownian` block.

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
