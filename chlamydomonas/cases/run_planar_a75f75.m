%RUN_PLANAR_A75F75  Coplanar Chlamydomonas with flagellum-2 amp & freq x0.75 (1 s)
%
%  Description:
%    Coplanar control case (no beat-plane tilt, no out-of-plane drive):
%    flagellum 1 (+Y) keeps the measured gait; flagellum 2 (-Y) uses 0.75x
%    oscillating amplitude and 0.75x beat frequency for 1.0 s (dt = 1 ms).
%    Results go to output/Chlamy-SB-A75F75/.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-14
%  Updated: 2026-11-24  refactored into chlamydomonas/cases/ via chlamy_run_preset

this_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(this_dir));
s = chlamy_run_preset('PLANAR_A75F75');
