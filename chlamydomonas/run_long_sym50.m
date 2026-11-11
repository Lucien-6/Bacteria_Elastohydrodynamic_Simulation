%RUN_LONG_SYM50  Long closed-loop case: symmetric biflagellate, 50 beats (open in MATLAB and press Run)
%
%  Description:
%    Symmetric reference (no beat-plane tilt, identical drive of both flagella)
%    of cell 1p00 A at 57.8 Hz for 50 beat periods (~0.865 s), with the 3-D
%    motion animation. Results go to output/Chlamy-SB-SYM50/ .
%    The script locates the project root from its own path, so it can be run
%    from any current folder. Needs output/chlamy_inverse/fit_1p00A.mat.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
chlamy_run_symbreak('SYM50');
