function s = chlamy_run_symbreak(tag, n_threads)
%CHLAMY_RUN_SYMBREAK  Closed-loop symmetry-breaking cases of the biflagellate Chlamydomonas model
%
%  Description:
%    Runs one closed-loop case driven by the analytic 2-harmonic fit of the
%    measured gait of cell 1p00 A (57.8 Hz), with a prescribed asymmetry,
%    and evaluates swimming speed, spin about the long axis and the screw
%    (helical) path of the one-beat window.
%      A  beat-plane tilt: both flagella rotated by gamma about their own
%         root tangents, same sense (chiral, produces spin)
%      B  cis/trans in-plane asymmetry: drive amplitude / phase of flagellum 2
%
%  Syntax:
%    s = chlamy_run_symbreak(tag)
%
%  Inputs:
%    tag - Case tag: 'T3','T6','T10' (tilt 3/6/10 deg), 'A85' (flagellum-2
%          amplitude x0.85), 'P30' (flagellum-2 phase +30 deg),
%          'T6A85' (tilt 6 deg and amplitude 0.85), 'T17' (tilt 17 deg, the
%          measured orbital tilt of Wan, PRL 126, 088003), 'T17A85' (tilt
%          17 deg and amplitude 0.85), 'T17A85L' (same as T17A85, 20 beats
%          and the 3-D animation), 'SYM50' (symmetric reference, 50 beats,
%          animation), 'ASYM50' (tilt 17 deg and amplitude 0.85, 50 beats,
%          animation)
%
%    n_threads - (Optional) computational threads for this MATLAB session,
%          limits oversubscription when several cases run in parallel
%
%  Output:
%    s - Struct from chlamy_swim_analysis (saved with the case outputs)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-10

if nargin >= 2 && ~isempty(n_threads), maxNumCompThreads(n_threads); end
root = fileparts(fileparts(mfilename('fullpath')));
cd(root);
addpath(genpath('src')); addpath('chlamydomonas'); warning('off', 'all');
S = load('output/chlamy_inverse/fit_1p00A.mat');

ov = struct('case_name', ['Chlamy-SB-' tag], 'T_obs', 0.053, 'moment_table', S.fit);
switch tag
    case 'T3',    ov.tilt_deg = [3 3];
    case 'T6',    ov.tilt_deg = [6 6];
    case 'T10',   ov.tilt_deg = [10 10];
    case 'A85',   ov.drive_amp = [1 0.85];
    case 'P30',   ov.drive_phase = [0 pi/6];
    case 'T17',   ov.tilt_deg = [17 17];
    case 'T17A85', ov.tilt_deg = [17 17]; ov.drive_amp = [1 0.85];
    case 'T17A85L'   % 20 beats, with the 3-D motion animation
        ov.tilt_deg = [17 17]; ov.drive_amp = [1 0.85];
        ov.T_obs = 20 / S.fit.f_Hz; ov.animation = true;
    case 'SYM50'     % symmetric reference, 50 beats, with the 3-D motion animation
        ov.T_obs = 50 / S.fit.f_Hz; ov.animation = true;
    case 'ASYM50'    % tilt 17 deg + flagellum-2 amplitude x0.85, 50 beats, with animation
        ov.tilt_deg = [17 17]; ov.drive_amp = [1 0.85];
        ov.T_obs = 50 / S.fit.f_Hz; ov.animation = true;
    case 'T6A85', ov.tilt_deg = [6 6]; ov.drive_amp = [1 0.85];
    otherwise,    error('Unknown tag %s', tag);
end

[res, ~, out_dir] = chlamy_run_case('configs/config_chlamy_biflagellate_gz.json', ov);
s = chlamy_swim_analysis(res, struct('f_Hz', S.fit.f_Hz, 'n_skip', 2));
save(fullfile(out_dir, 'swim_analysis.mat'), 's');
end
