# Theory and Numerical Model

English | [中文](theory.zh-CN.md)  
[Back to README](../README.md)

> **Author:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **Updated:** 2026-05-25
> **Version:** 7.4.5

This document contains the model details that were formerly embedded in the long README. The repository entry page is now intentionally short; practical run instructions are in [`../USER_GUIDE.md`](../USER_GUIDE.md), and release history is in [`../CHANGELOG.md`](../CHANGELOG.md).

## 1. Physical Problem

The simulator models a bacterial body with flexible pili or flagella settling in a Stokes-flow environment:

- **Body:** rigid capsule, sphere-like body, or imported custom 3D point cloud.
- **Pili / flagella:** one or more inextensible Kirchhoff elastic rods.
- **Fluid:** zero-Reynolds-number Stokes flow.
- **Coupling:** hydrodynamic interactions between all body and pili points are retained through a global mobility solve.

The default parametric body is a capsule with two hemispherical caps and a cylindrical midsection. In custom geometry mode, both the body and pili are imported from user-provided point clouds; mixing custom and parametric geometry in either direction is rejected by `validate_simulation_params.m` when the mixed fields reach the validator (inline or programmatic `cfg`, or JSON loaded into `cfg` and passed through `run_simulation.m`). Note: in JSON parametric mode, `load_sim_config.m` does **not** read a top-level `custom_geometry` block into `cfg`, so such stray keys are silently ignored at load time rather than validated—remove them from the JSON file or switch to `geometry_mode: "custom"`.

## 2. Main Methods

| Component | Method | Main source |
| --- | --- | --- |
| Body hydrodynamics | Regularized Stokeslet Method (RSM) on body surface points | Cortez et al. (2005) |
| Pili hydrodynamics | Regularized Stokeslet Segments (RSS) with exact segment integration | Cortez (2018) |
| Pili twist-flow coupling | Regularized rotlet segments | Cortez (2018) |
| Pili elasticity | Kirchhoff elastic rod with Euler-angle directors | Walker (2020) |
| Body orientation | Quaternion kinematics | Project implementation |
| Time integration | MATLAB `ode15s` with event-driven segmenting | Project implementation |

The key numerical advantage is that RSS integrates the regularized Stokeslet kernel along each straight pilus segment exactly. This separates segment length from regularization radius, allowing roughly `10` segments per micrometer for typical pili while keeping hydrodynamic interactions resolved.

## 3. Nondimensionalization

The body semi-major axis is used as the reference length:

| Quantity | Symbol | Definition |
| --- | --- | --- |
| Length | $L_0$ | body `major_axis` (parametric); `L_0_override` or body-derived suggestion (custom) |
| Bending stiffness | $EI$ | $E\pi r_{\mathrm{pili}}^4/4$ |
| Time | $T_0$ | $8\pi\mu L_0^4/EI$ |
| Velocity | $U_0$ | $L_0/T_0$ |
| Force | $F_0$ | $EI/L_0^2$ |

Important dimensionless parameters are:

- Elastic-gravity number: $E_g = F_{\mathrm{sink}}L_0^2/EI$.
- Pili slenderness: $\varepsilon = r_{\mathrm{pili}}/L_0$.
- Poisson ratio: $\sigma$, supplied by the user.

The sinking force (dimensional) is

```math
F_{\mathrm{sink}} = (\rho_{\mathrm{bac}} - \rho_{\mathrm{fluid}})\, V_{\mathrm{body}}\, g,
```

where $V_{\mathrm{body}}$ is the capsule volume in parametric mode or the alpha-shape volume in custom mode. Gravity enters the ODE right-hand side as $\mathbf{R}_{1:3} = -E_g\,\hat{\mathbf{g}}$ with $\hat{\mathbf{g}} = [0,0,-1]^\mathsf{T}$.

The root–body gap (dimensional) is

```math
\mathrm{dis} = 1.1\,(\varepsilon_{\mathrm{body}} + \varepsilon_{\mathrm{pili}})\,L_0,
```

passed to `pili_building` / `custom_pili_building` to keep regularized kernels separated at the pilus root.

Additional model inputs use the same reference scales:

- Background flow velocity: $\mathbf{u}_{\mathrm{bg}}^\ast = \mathbf{u}_{\mathrm{bg}}/U_0$.
- Integrated active moments: $\mathbf{M}^\ast = \mathbf{M}_{\mathrm{phys}}L_{\mathrm{pili}}/EI$.
- Intrinsic curvature: stored on the same nondimensional arclength scale as the Euler-angle curvature fields.

## 4. State Vector and Geometry

The full state contains body position, body orientation, and each pilus rod state:

```math
\mathbf{Z} =
\left[
\mathbf{X}_{\mathrm{cm}},
\mathbf{q},
\boldsymbol{\theta}_1,
\boldsymbol{\phi}_1,
\boldsymbol{\psi}_1,
\ldots,
\boldsymbol{\theta}_M,
\boldsymbol{\phi}_M,
\boldsymbol{\psi}_M
\right].
```

The body quaternion $\mathbf{q}$ represents the body-to-lab rotation and is advanced with:

```math
\frac{d\mathbf{q}}{dt}
= \frac{1}{2}
\begin{bmatrix}
0 \\
\boldsymbol{\Omega}
\end{bmatrix}
\otimes \mathbf{q}.
```

Here $\boldsymbol{\Omega}$ is the **body angular velocity in the lab frame**, taken from `z_dot(4:6)` of the coupled solve (`quat_kinematics.m`).

For each pilus, Euler angles define the material directors `(d1, d2, d3)`, where `d3` is the segment tangent. [`reconstruct_geometry.m`](../src/geometry/reconstruct_geometry.m) maps `Z` to lab-frame body points, pilus nodes, directors, and auxiliary geometry used by the hydrodynamic solve.

## 5. Hydrodynamics

### 5.1 Body RSM

Body surface forces are represented by the Cortez regularized Stokeslet kernel:

```math
\mathbf{S}^{\varepsilon}(\mathbf{r})
=
\frac{
\left(\lVert\mathbf{r}\rVert^2 + 2\varepsilon^2\right)\mathbf{I}
+ \mathbf{r}\mathbf{r}^{\mathsf{T}}
}{
\left(\lVert\mathbf{r}\rVert^2 + \varepsilon^2\right)^{3/2}
}.
```

The body regularization length defaults to:

```math
\varepsilon_{\mathrm{body}} = h_{\mathrm{mean}}/2.6667.
```

where `h_mean` is the mean nearest-neighbor spacing of the body surface points in nondimensional units. This default was calibrated by the `N_body x eps_body` Stokes-sphere sweep and can be overridden by `body.eps_body` in JSON or `cfg.eps_body` in MATLAB.

### 5.2 Pili RSS

For a straight segment from $\mathbf{y}_0$ to $\mathbf{y}_1$, RSS evaluates the velocity contribution of a linearly varying force density:

```math
\mathbf{u}(\mathbf{x})
=
\mathbf{M}_1\mathbf{f}_0
+ \mathbf{M}_2\mathbf{f}_1.
```

The matrices `M1` and `M2` include exact segment integration weights. No additional trapezoidal length factor is applied to the RSS segment columns. The required `T_{n,q}` integrals are evaluated by the Cortez (2018) recurrence.

For pili, the regularization length is physical:

```math
\varepsilon_{\mathrm{pili}} = r_{\mathrm{pili}}/L_0.
```

Unlike body RSM regularization, this value represents the finite radius of the pilus and is not a numerical parameter to send toward zero.

### 5.3 Rotlet Segments

Twist-induced flow is included through regularized rotlet segments:

```math
\mathbf{u}(\mathbf{x})
=
-L
\left(
c_a[\mathbf{x}_0]_{\times}
+ c_b[\mathbf{v}]_{\times}
\right)
\boldsymbol{\tau}.
```

where `tau` is the segment torque density. For typical bacterial pili, this coupling is very small because it scales with high powers of slenderness, but it keeps the formulation theoretically complete.

## 6. Coupled Linear System

At each ODE right-hand-side call, the simulator assembles:

```math
\mathbf{L}\dot{\mathbf{z}} = \mathbf{R},
```

with:

```math
\mathbf{L} = -\mathbf{B}\mathbf{A}\mathbf{Q},
```

```math
\mathbf{A}\mathbf{Q}
=
\begin{bmatrix}
\mathbf{M}_{\mathrm{weighted}}^{-1}\mathbf{Q}_{\mathrm{trans,corr}} \\
\widetilde{\mathbf{A}}_{\mathrm{twist}}\mathbf{Q}_{\mathrm{twist}}
\end{bmatrix},
```

```math
\mathbf{Q}_{\mathrm{trans,corr}}
=
\mathbf{Q}_{\mathrm{trans}}
- \mathbf{M}_{\mathrm{rotlet}}
\widetilde{\mathbf{A}}_{\mathrm{twist}}
\mathbf{Q}_{\mathrm{twist}}.
```

The main blocks are:

- `Q_global`: maps generalized rates to body and pilus node velocities.
- `M_weighted`: full hydrodynamic mobility matrix, including body RSM and pili RSS blocks.
- `A_twist`: torsional resistance operator.
- `B_global`: maps distributed forces and torques into global force, moment, and pilus moment-balance equations.
- `R`: gravitational force, elastic moments, active moments, background-flow corrections, and motor torque terms.

Solves use MATLAB backslash / LU factorization rather than explicit matrix inversion.

## 7. Elastic Rod Moments

Pili are treated as Kirchhoff rods. The Darboux vector is:

```math
\mathbf{U}
=
\kappa_1\mathbf{d}_1
+ \kappa_2\mathbf{d}_2
+ \kappa_3\mathbf{d}_3,
```

with director evolution:

```math
\mathbf{d}_{\alpha}' = \mathbf{U}\times\mathbf{d}_{\alpha}.
```

The curvature components used in [`compute_pili_elastic_moments.m`](../src/physics/compute_pili_elastic_moments.m) are:

```math
\kappa_1
=
\mathbf{d}_3\cdot\mathbf{d}_2'
- \kappa_{1,0},
```

```math
\kappa_2
=
\mathbf{d}_1\cdot\mathbf{d}_3'
- \kappa_{2,0},
```

```math
\kappa_3
=
\frac{
\mathbf{d}_2\cdot\mathbf{d}_1'
- \kappa_{3,0}
}{1+\sigma}.
```

Intrinsic curvature can be configured per pilus as `frozen`, `zero`, `constant`, or `helix`. Active moments can be disabled, prescribed as a traveling wave, or supplied by custom configuration.

## 8. Root Coupling and Motor Torque

Each pilus root supports two boundary modes:

- `clamped`: root bending and twist angles are locked to the body frame.
- `free_torsion`: root bending angles remain locked, but the root twist angle `psi_1` evolves dynamically.

In `free_torsion`, the root `d3` moment row may receive elastic torsion, active `d3` moment, and constant motor torque. As of v7.3.0, no explicit motor reaction term is added to `R_global(4:6)`. Those rows are the angular-momentum balance for the entire coupled body-pilus hydrodynamic system, so a motor torque exchanged internally between body and pilus cancels in the system sum. The body angular response is still present, but it emerges implicitly from the coupled Stokes solve and the root kinematic constraints.

## 9. Background Flow

Background flow is evaluated as a prescribed steady velocity field and incorporated through the hydrodynamic correction:

```math
\mathbf{R}
\leftarrow
\mathbf{R}
- \mathbf{B}\mathbf{A}^{-1}\mathbf{u}_{\mathrm{bg}}.
```

Supported flow types are:

- `off`: no background flow.
- `uniform`: constant velocity with user-defined direction.
- `shear_z`: finite-thickness clipped shear layer with user-defined flow direction, shear axis, and thickness.

## 10. Brownian Dynamics

When Brownian motion is enabled, the deterministic elastohydrodynamic solve is operator-split from stochastic body motion:

1. Integrate the deterministic coupled system with `ode15s`.
2. At each Brownian time node, compute the full-system effective body mobility.
3. Apply a correlated 6D translational-rotational Brownian increment.

The effective mobility is computed by probing the coupled linear system with six unit generalized body loads:

```math
\mathbf{M}_{\mathrm{eff}}^{\ast}(:,i)
=
\left(\mathbf{L}^{-1}(-\mathbf{e}_i)\right)_{1:6},
\qquad i=1,\ldots,6.
```

After dimensional conversion:

```math
\mathbf{D}_{\mathrm{eff}}
=
k_B T\,\mathbf{M}_{\mathrm{eff,dim}},
```

```math
\mathbf{L}_{\mathrm{eff}}
\mathbf{L}_{\mathrm{eff}}^{\mathsf{T}}
=
2\mathbf{D}_{\mathrm{eff}},
```

```math
\delta\mathbf{U}
=
\mathbf{L}_{\mathrm{eff}}
\boldsymbol{\eta}
\sqrt{\Delta t_B},
\qquad
\boldsymbol{\eta}\sim\mathcal{N}(\mathbf{0},\mathbf{I}_6).
```

The optional `fixman_midpoint` mode recomputes the mobility at a stochastic midpoint to approximate the configuration-dependent Ito drift. This doubles the mobility-evaluation cost per Brownian step.

## 11. Flow-Field Visualization

[`compute_plot_flowfield.m`](../src/visualization/compute_plot_flowfield.m) reconstructs the end-state geometry, evaluates the flow on XZ and YZ planes around the body center, masks body interiors, and exports a 2 x 2 panel figure:

- 3D XZ plane view.
- 3D YZ plane view.
- 2D XZ velocity magnitude and normalized direction arrows.
- 2D YZ velocity magnitude and normalized direction arrows.

Custom bodies use stored triangulation and alpha-shape based interior masking where available.

## 12. Walker Scaling Consistency

When the coupled problem is reduced to a single free elastic filament without the body, the nondimensional structure is designed to remain consistent with the Walker-style Kirchhoff-rod formulation:

- The spatial $Q$ matrix has the same role: it maps body-frame or rod-frame generalized rates into translational and angular velocities along the filament.
- The force-balance rows use trapezoidal quadrature weights, while the pilus moment-balance rows use the per-rod length scaling.
- Elastic and active moments are stored in the per-pilus arclength coordinate $s_k\in[0,1]$.

For torsional resistance, the dimensional segment relation is proportional to the local segment length:

```math
\tau_{\mathrm{twist},j}
\sim
c_t\,\Delta s_k\,\dot{\psi}_j,
\qquad
c_t \sim 4\pi\mu r_{\mathrm{pili}}^2.
```

Using the global project scales gives the nondimensional torsional coefficient:

```math
\frac{c_t\,\Delta s_k/T_0}{EI/L_0}
\sim
\frac{\varepsilon^2}{2}\,\Delta \hat{s}_k.
```

The implementation stores each pilus on its own normalized arclength, so the conversion between the global body scale and the per-pilus Walker scale introduces the pilus length ratio $\ell_k = L_{\mathrm{pili},k}/L_0$. Per-pilus slenderness is $\varepsilon_k = \varepsilon/\ell_k = r_{\mathrm{pili}}/L_{\mathrm{pili},k}$ (`build_global_resistance.m`).

```math
\widetilde{\mathbf{A}}_{\mathrm{twist},k}
=
\frac{\ell_k\varepsilon_k^2}{2}
\widehat{\mathbf{A}}_k
=
\frac{\varepsilon^2}{2\ell_k}
\widehat{\mathbf{A}}_k.
```

The two right-hand forms are algebraically equivalent because $\ell_k\varepsilon_k^2 = \varepsilon^2/\ell_k$.

The complementary $\ell_k^3$ factor used in each pilus moment-balance row converts the per-rod moment equations back to the global scale, keeping the assembled system

```math
\mathbf{L}\dot{\mathbf{z}} = \mathbf{R}
```

dimensionally homogeneous.

## 13. Validation Scope

The codebase has been checked against the following model-level invariants:

- Quaternion convention and body orientation kinematics.
- RSM kernel scaling and body-volume consistency.
- RSS recurrence terms and segment matrix symmetry checks.
- Rotlet anti-symmetry and recurrence behavior.
- Walker-style Euler-angle `Q` matrix and moment-balance scaling.
- Clamped and free-torsion root boundary constraints.
- Multi-pilus state layout, geometry reconstruction, and matrix dimensions.
- Full-system Brownian effective mobility and Cholesky positive-definiteness checks.
- Output units for MAT, CSV, figure, and log files.

Detailed bug-fix history belongs in [`../CHANGELOG.md`](../CHANGELOG.md), not in this theory document.

## 14. References

- Cortez, R. (2005). The method of regularized Stokeslets.
- Cortez, R. (2018). Regularized Stokeslet segments.
- Walker, B. J. (2020). Elastic filament / bacterial appendage formulation used as the rod-model reference.

## 15. License

MIT License. See [`../LICENSE`](../LICENSE).
