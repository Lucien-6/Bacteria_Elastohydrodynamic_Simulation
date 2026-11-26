# CPC Program Summary — BacEHD

Ready-to-paste *Computer Physics Communications* (CPC) Program Summary metadata for **BacEHD** (`v10.2.2`).

---

- **Program title:** BacEHD
- **CPC Library link to program files:** *(to be assigned by Elsevier CPC upon acceptance)*
- **Developer's repository link:** https://github.com/Lucien-6/Bacteria_Elastohydrodynamic_Simulation
- **Code Ocean capsule:** *(optional)*
- **Licensing provisions:** MIT License (see `LICENSE`; bundled third-party S2 Sampling Toolbox by Anton Semechko is also MIT licensed — see `MATLAB Add-On/Suite of functions to perform uniform sampling of a sphere/LICENSE.md`).
- **Programming language:** MATLAB (R2022a or later; developed and tested on R2026a on Windows 11 and R2023a on Linux HPC). No Python or external compiled MEX dependencies are required to run simulations.
- **Supplementary material:** Reviewer verification suite (`cpc/run_cpc_tests.m`), user manual (`USER_GUIDE.md`), external hydrodynamic/elastic benchmark drivers (`validation/`), performance and parallel/GPU scaling drivers (`validation/cpc_performance_benchmark.m`, `validation/cpc_parallel_and_flowfield_benchmark.m`), example JSON configurations (`configs/`), and *Chlamydomonas* 3-D helical swimmer module (`chlamydomonas/`).
- **Journal reference of previous version:** *(Not applicable — first release in CPC)*
- **Does the new version supersede the previous version?:** *(Not applicable)*
- **Nature of problem:**
  Simulating the low-Reynolds-number three-dimensional elastohydrodynamics and Brownian motion of a rigid microbial cell body coupled to one or more flexible, inextensible slender appendages (pili or flagella) in an unbounded Newtonian Stokes fluid. Applications span gravitational sedimentation of multi-pilus bacteria, steady/shear background flows, motor-driven bacterial propulsion (*Escherichia coli*), and active-moment-driven eukaryotic swimming (*Chlamydomonas reinhardtii*), covering four orders of magnitude in bending stiffness (\(EI \sim 2\times 10^{-26}\) to \(8.4\times 10^{-22}\,\mathrm{N\cdot m^2}\)).
- **Solution method:**
  The cell body surface is discretized by regularized Stokeslets (RSM, Cortez et al. 2005) while each appendage is modeled as an inextensible Kirchhoff rod coupled to the fluid by exact regularized Stokeslet and rotlet segment integrals (RSS, Cortez 2018) with a calibrated hydrodynamic radius factor (\(\beta = \varepsilon_{\mathrm{pili}}/r_{\mathrm{pili}} = \sqrt{2}\) matching no-slip cylinder axial twist drag) and weighted-reciprocal non-local twist mobility. Second-order node-collocated SVD directors evaluate Kirchhoff bending and twist moments; Euler-angle coordinate singularities are avoided via adaptive proper-rotation basis changes. Root boundary conditions support clamped attachments and free-torsion rotary motors with an exact in-call self-consistent signed linear torque–speed solve. Thermal fluctuations are realized either as a body-only 6-DOF constant generalized thermal force per Brownian interval (producing implicit Itô drift without root kinks) or jointly with the slow elastohydrodynamic pilus bending modes using constrained-probe centerline force densities and balancing axial torque densities (`brown_pili_tau`) with held-load drift compensation. Time integration uses MATLAB's stiff BDF/NDF solver `ode15s`, accelerated by rigid-body body-frame Block Schur complement precomputation (`build_body_precomp`, `build_pili_mobility_blocks`), a 3-colour frozen-hydrodynamics analytical/sparse Jacobian (`ode_coupled_jacobian`), and a lagged Jacobian cache (`LaggedJacobian`) across Brownian restarts.
- **Additional comments including restrictions and unusual features:**
  - *Restrictions:* Assumes an unbounded quiescent or prescribed linear/uniform Newtonian Stokes fluid (no wall boundaries, no fluid inertia, and no steric contact/excluded-volume repulsion between filaments and cell body; clearance is monitored a posteriori in swimming analyses). The coupled pilus thermalisation (`brownian.pili_thermal`) is formulated for passive, naturally straight, clamped filaments with \(\tau_{\mathrm{det}} \lesssim \Delta t_B \lesssim \tau_{\mathrm{slow}}/4\).
  - *External routines/libraries:* Capsule surface mesh generation uses the bundled MIT-licensed S2 Sampling Toolbox by Anton Semechko (`ParticleSampleSphere`, `SpiralSampleSphere`, `IcosahedronMesh`, `QuadCubeMesh`, `SubdivideSphericalMesh`). Optional multi-worker batch sweeps (`batch_sweep/`) use MATLAB Parallel Computing Toolbox when available and fall back to serial execution otherwise. Flow-field plotting uses pure base-MATLAB percentile interpolation (no Statistics Toolbox dependency).
  - *Unusual features:* Built-in runtime numerical diagnostics (scale-independent linear-solve backward residuals, condition-number monitoring, weighted reciprocity and dissipation checks), self-contained JSON configuration snapshots, platform-aware MP4/AVI 3-D animation export, and 6-panel 3-D/2-D streamline flow-field reconstruction from saved states.
- **Running time:**
  Measured on MATLAB R2026a Update 5 (Windows 11, 24-core workstation, NVIDIA GeForce RTX 5080 for optional GPU flow-field post-processing):
  - Full 18-check reviewer verification suite (`cpc/run_cpc_tests('SkipExternal', true)`): **22.5 s** (the 12 core physics checks take ~3.2 s).
  - External hydrodynamic/elastic validation drivers (`validation/run_cpc_external_checks`): **~1–2 min** (helix propulsive matrix + cantilever bending relaxation).
  - Operator & end-to-end performance benchmark (`validation/cpc_performance_benchmark`): **17.2 s**. Single RHS evaluation takes **0.16 ms** (bare capsule, \(N_{\mathrm{body}}=306\), **48.8×** faster than dense LU), **5.13 ms** (1 pilus, 30 segments, **4.13×**), **29.1 ms** (4 pili, 120 segments, **1.46×**), and **49.6 ms** (4 pili, \(N_{\mathrm{body}}=600\), **1.81×**); single Jacobian evaluation takes **18.5 ms** (1 pilus, **72.5×** faster than full FD) to **59.3 ms** (4 pili, \(N_{\mathrm{body}}=600\), **538.6×** faster). Short 0.5 s deterministic settling (1 pilus, \(N_{\mathrm{body}}=306\)) completes in **0.19 s** (**9.60×** speedup vs legacy dense+FD); 0.2 s Brownian trajectory completes in **1.67 s** (**5.13×** speedup).
  - Multi-worker Brownian ensemble (`validation/cpc_parallel_and_flowfield_benchmark`, 32 independent 0.25 s runs, 1 BLAS thread/worker): **77.7 s** (1 worker) → **13.7 s** (16 workers, **5.65×** speedup). 12-frame 3-plane *E. coli* flow-field reconstruction (150,750 exterior points): **6.62 s** (serial CPU) → **1.52 s** (RTX 5080 FP64 GPU + 8-thread pool, **4.36×** speedup; FP64 relative error \(4.3\times 10^{-14}\)).

---

Author: Lucien (lucien-6@qq.com)
Created: 2026-11-24
