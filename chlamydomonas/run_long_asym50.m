%RUN_LONG_ASYM50  Long closed-loop case: asymmetric biflagellate, 50 beats (open in MATLAB and press Run)
%
%  Description:
%    Symmetry-broken cell 1p00 A at 57.8 Hz for 50 beat periods (~0.865 s),
%    with the 3-D motion animation: both beat planes tilted by 17 deg about
%    their root tangents (same sense) and the drive amplitude of flagellum 2
%    scaled by 0.85. Results go to output/Chlamy-SB-ASYM50/ .
%    The script locates the project root from its own path, so it can be run
%    from any current folder. Needs output/chlamy_inverse/fit_1p00A.mat.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
chlamy_run_symbreak('ASYM50');
