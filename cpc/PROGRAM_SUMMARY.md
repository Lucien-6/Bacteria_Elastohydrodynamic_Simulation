# CPC Program Summary — BacEHD

Ready-to-paste *Computer Physics Communications* (CPC) Program Summary metadata for **BacEHD** (`v10.1.0`).

---

- **Program title:** BacEHD
- **CPC Library link to program files:** *(to be assigned by Elsevier CPC upon acceptance)*
- **Developer's repository link:** *(insert GitHub / Zenodo DOI link before submission)*
- **Code Ocean capsule:** *(optional)*
- **Licensing provisions:** MIT License (see `LICENSE`; bundled third-party S2 Sampling Toolbox by Anton Semechko is also MIT licensed — see `MATLAB Add-On/Suite of functions to perform uniform sampling of a sphere/LICENSE.md`).
- **Programming language:** MATLAB (R2022a or later; developed and tested on R2026a on Windows 11 and R2023a on Linux HPC). No Python or external compiled MEX dependencies are required to run simulations.
- **Supplementary material:** Reviewer verification suite (`cpc/run_cpc_tests.m`), external hydrodynamic/elastic benchmark drivers (`validation/`), performance scaling driver (`validation/cpc_performance_benchmark.m`), example JSON configurations (`configs/`), and *Chlamydomonas* 3-D helical swimmer module (`chlamydomonas/`).
- **Journal reference of previous version:** *(Not applicable — first release in CPC)*
- **Does the new version supersede the previous version?:** *(Not applicable)*
- **Nature of problem:**
  Simulating the low-Reynolds-number three-dimensional elastohydrodynamics and Brownian motion of a rigid microbial cell body coupled to one or more flexible, inextensible slender appendages (pili or flagella) in an unbounded Newtonian Stokes fluid. Applications span gravitational sedimentation of multi-pilus bacteria, steady/shear background flows, motor-driven bacterial propulsion (*Escherichia coli*), and active-moment-driven eukaryotic swimming (*Chlamydomonas reinhardtii*), covering four orders of magnitude in bending stiffness (\(EI \sim 2\times 10^{-26}\) to \(8.4\times 10^{-22}\,\mathrm{N\cdot m^2}\)).
- **Solution method:**
  The cell body surface is discretized by regularized Stokeslets (RSM, Cortez et al. 2005) while each appendage is modeled as an inextensible Kirchhoff rod coupled to the fluid by exact regularized Stokeslet and rotlet segment integrals (RSS, Cortez 2018) with a calibrated hydrodynamic radius factor (\(\beta = \varepsilon_{\mathrm{pili}}/r_{\mathrm{pili}} = \sqrt{2}\) matching no-slip cylinder axial twist drag) and weighted-reciprocal non-local twist mobility. Second-order node-collocated SVD directors evaluate Kirchhoff bending and twist moments; Euler-angle coordinate singularities are avoided via adaptive proper-rotation basis changes. Root boundary conditions support clamped attachments and free-torsion rotary motors with an exact in-call self-consistent signed linear torque–speed solve. Thermal fluctuations are realized either as a body-only 6-DOF constant generalized thermal force per Brownian interval (producing implicit Itô drift without root kinks) or jointly with the slow elastohydrodynamic pilus bending modes using constrained-probe centerline force densities and balancing axial torque densities (`brown_pili_tau`) with held-load drift compensation. Time integration uses MATLAB's stiff BDF/NDF solver `ode15s`, accelerated by rigid-body body-frame Block Schur complement precomputation (`build_body_precomp`, `build_pili_mobility_blocks`), a 3-colour frozen-hydrodynamics analytical/sparse Jacobian (`ode_coupled_jacobian`), and a lagged Jacobian cache (`LaggedJacobian`) across Brownian restarts.
- **Additional comments including restrictions and unusual features:**
  - *Restrictions:* Assumes an unbounded quiescent or prescribed linear/uniform Newtonian Stokes fluid (no wall boundaries, no fluid inertia, and no steric contact/excluded-volume repulsion between filaments and cell body; clearance is monitored a posteriori in swimming analyses). The coupled pilus thermal fluctuation scheme (`brownian.pili_thermal`) is formulated and validated for passive, naturally straight, clamped filaments satisfying adiabatic scale separation (\(\Delta t_B \le \tau_{\mathrm{fast}}/4\)).
  - *External routines/libraries:* Capsule surface mesh generation uses the bundled MIT-licensed S2 Sampling Toolbox by Anton Semechko (`ParticleSampleSphere`, `SpiralSampleSphere`, `IcosahedronMesh`, `QuadCubeMesh`, `SubdivideSphericalMesh`). Optional multi-worker batch sweeps (`batch_sweep/`) use MATLAB Parallel Computing Toolbox when available and fall back to serial execution otherwise. Flow-field plotting uses pure base-MATLAB percentile interpolation (no Statistics Toolbox dependency).
  - *Unusual features:* Built-in runtime numerical diagnostics (scale-independent linear-solve backward residuals, condition-number monitoring, weighted reciprocity and dissipation checks), self-contained JSON configuration snapshots, platform-aware MP4/AVI 3-D animation export, and 6-panel 3-D/2-D streamline flow-field reconstruction from saved states.
- **Running time:**
  Measured on MATLAB R2026a (Windows 11, single-threaded BLAS unless noted):
  - Full 18-check reviewer verification suite (`cpc/run_cpc_tests('SkipExternal', true)`): **~22 s** (the 12 core physics checks take **~3.2 s**).
  - External hydrodynamic/elastic validation drivers (`validation/run_cpc_external_checks`): **~1–2 min** (helix propulsive matrix + cantilever bending relaxation).
  - Operator & end-to-end performance benchmark (`validation/cpc_performance_benchmark`): **~18 s**. Single RHS evaluation takes **0.09 ms** (bare capsule, \(N_{\mathrm{body}}=306\), 99× faster than dense LU), **5.1 ms** (1 pilus, 30 segments, 3.8×), and **30.7 ms** (4 pili, 120 segments, 1.5×); single Jacobian evaluation takes **18.6 ms** (1 pilus, 75× faster than full FD) to **67.6 ms** (4 pili, \(N_{\mathrm{body}}=600\), 517× faster). Short 0.5 s deterministic settling (1 pilus, \(N_{\mathrm{body}}=306\)) completes in **0.19 s** (9.9× speedup vs legacy dense+FD); 0.2 s Brownian trajectory completes in **1.66 s** (5.3× speedup).

---

Author: Lucien (lucien-6@qq.com)
Created: 2026-11-24
