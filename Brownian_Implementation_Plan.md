# Brownian dynamics — implementation notes

**Author:** Lucien (lucien-6@qq.com)  
**Created:** 2026-04-21

This note complements `README.md` §16 and the Brownian-related entries in `CHANGELOG.md`. It records modelling scope, stochastic interpretation, and where warnings are emitted in code.

---

## 1.5 Scope & simplifying assumptions

**What is modelled**

- **Translational and rotational Brownian motion of the rigid bacterial body** (6 degrees of freedom) using an effective mobility tensor obtained from the **full coupled body–pili linear system** (probing / Schur-complement viewpoint in `compute_effective_body_brownian.m`).
- **Hydrodynamic drag on pili** enters only through their effect on the coupled Stokes system that defines the effective body mobility. There is **no** explicit Langevin noise on internal pili degrees of freedom in the current implementation.

**What is not modelled**

- **Thermal fluctuations of pili shape** (bending/twist internal modes) are omitted. The pili are still deterministic elastic filaments driven by elasticity, clamped boundary conditions, and optional active moments; only the **body** receives the stochastic kicks.
- Analyses that require **fluctuating filament conformations** (e.g. tension–extension statistics along the pilus) are outside the present scope.

**Extensions (literature pointers)**

- Filament-level Brownian elastohydrodynamics and related coupling ideas appear, for example, in Westwood & Keaveny (2022) and Maxian & Donev (2023). The present codebase does not implement those extensions.

---

## 2.4.1 Stochastic interpretation & missing drift

**Overdamped Langevin form**

In configuration space, the **exact Itô** overdamped Langevin equation for the body with a configuration-dependent mobility tensor \( \mathbf{M}_{\mathrm{eff}}(\mathbf{q}) \) includes a drift correction proportional to the divergence of the mobility:

\[
\mathrm{d}\mathbf{q} = \mathbf{M}_{\mathrm{eff}}\,\mathbf{f}\,\mathrm{d}t + k_B T\,\nabla_{\mathbf{q}}\cdot \mathbf{M}_{\mathrm{eff}}\,\mathrm{d}t + \sqrt{2 k_B T}\,\mathbf{L}_{\mathrm{eff}}\,\mathrm{d}\mathbf{W},
\]

where \( \mathbf{M}_{\mathrm{eff}} = \mathbf{L}_{\mathrm{eff}}\mathbf{L}_{\mathrm{eff}}^{\mathsf{T}} \) (Cholesky factor used in the integrator).

**`fixman_midpoint = false`**

The simple **Euler–Maruyama** step implemented in this mode uses the Cholesky factor at the beginning of the step but **drops** the \( k_B T\,\nabla\cdot \mathbf{M}_{\mathrm{eff}} \) drift term. When pili strongly modulate \( \mathbf{M}_{\mathrm{eff}} \), that omission can bias **long-time** diffusion statistics (effective diffusion coefficients, MSD fits) even if short-time trajectories look reasonable.

**`fixman_midpoint = true`**

The **Fixman midpoint** correction evaluates mobility-related quantities at an intermediate point in the step so that the discrete scheme recovers the correct **thermal drift** in the overdamped limit (at the cost of an extra linear-system / mobility evaluation per Brownian sub-step).

**Practical decision matrix (rule of thumb)**

| Situation | Recommendation |
|----------|------------------|
| Brownian off | `fixman_midpoint` irrelevant. |
| Brownian on, **no** pili (or very stiff / weakly coupled pili) | `fixman_midpoint = false` often sufficient for exploratory runs. |
| Brownian on, **multiple** pili or large-amplitude pili motion; long \(T_{\mathrm{obs}}\); MSD / diffusion metrics matter | Prefer `fixman_midpoint = true`. |

See `validate_simulation_params.m` for the concrete warning when `enable_brown` is true, pili are present, and `fixman_midpoint` is false.

---

## 7 Related source files (validator & batch tooling)

**`validate_simulation_params.m`**

- Emits a **warning** when Brownian is enabled, `Pili_Matrix` is non-empty, and `fixman_midpoint` is false (missing Itô drift in simple EM; see §2.4.1).
- Emits a **warning** when `dt_B` is large compared to a Stokes–Einstein spherical estimate of \( \min(\tau_{\mathrm{trans}}, \tau_{\mathrm{rot}}) \) (accuracy of Brownian statistics).

**`batch sweep/batch_repeat_brownian.m`**

- Repeats the same Brownian-on configuration with controlled RNG policy (`brown_rng_seed_mode`: `timestamp` vs `fixed`).

---

## Upgrade note (v4.0.0)

Brownian dynamics was upgraded from a **body-only** effective friction approximation to the **full-system** effective mobility obtained from the coupled body–pili linear system (`build_coupled_lin_sys.m` + `compute_effective_body_brownian.m`). Legacy `compute_body_brownian_params.m` was removed; see `CHANGELOG.md` [4.0.0].
