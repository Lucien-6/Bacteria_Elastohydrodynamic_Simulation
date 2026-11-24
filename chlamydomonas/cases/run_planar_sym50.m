%RUN_PLANAR_SYM50  Coplanar symmetric Chlamydomonas reference (50 beats)
%
%  Description:
%    Coplanar symmetric breaststroke reference (no beat-plane tilt, no
%    out-of-plane drive, identical mirror-symmetric flagella) at 57.75 Hz
%    for 50 beat periods (~0.866 s) with 3-D motion animation.
%    Results go to output/Chlamy-SB-SYM50/.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-11
%  Updated: 2026-11-24  refactored into chlamydomonas/cases/ via chlamy_run_preset

this_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(this_dir));
s = chlamy_run_preset('PLANAR_SYM50');
