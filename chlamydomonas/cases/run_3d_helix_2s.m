%RUN_3D_HELIX_2S  Realistic 3-D helical Chlamydomonas swimming (2.0 s)
%
%  Description:
%    Open this script in MATLAB and press Run (or call run('chlamydomonas/cases/run_3d_helix_2s.m')).
%    Simulates 2.0 s (~115.5 beat cycles, ~2.2 full left-handed superhelical
%    revolutions) of biflagellate Chlamydomonas reinhardtii driven by the
%    calibrated 3-D active moment:
%      * Measured-gait dual-k 1w+2w in-plane drive (cell 1p00 A, 57.75 Hz,
%        chlamydomonas/data/fit_1p00A.mat)
%      * C2-chiral beat-plane tilt: tilt_deg = [+17, +17] deg
%      * Proximal 1w+2w out-of-plane bending: r_perp = 0.25, r_perp2 = 0.05,
%        s_w = 0.20, delta = delta2 = +pi/2, perp_sign = [+1, -1]
%      * Physiological 5% cis/trans amplitude asymmetry: drive_amp = [1, 0.95]
%      * High-accuracy ODE tolerances: RelTol = AbsTol = 1e-4, dt = 1 ms
%      * Live progress written to output/Chlamy-3D-Helix-2s/progress.txt
%      * 3-D motion MP4 animation enabled; single-segment wall limit = 120 h
%
%    Expected steady-state kinematics (from 3-beat tol=1e-4 benchmark):
%      U_fwd ~ 77.6 um/s | spin ~ -1.09 Hz | left-handed helix R ~ 1.5 um,
%      pitch ~ -70 um, period ~ 0.92 s | clearance f_body_min ~ 1.124 (> 1).
%
%  Output directory:
%    output/Chlamy-3D-Helix-2s/
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

this_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(this_dir));
s = chlamy_run_preset('HELIX3D_2S');
