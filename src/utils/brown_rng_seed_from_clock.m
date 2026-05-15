function seed = brown_rng_seed_from_clock()
%BROWN_RNG_SEED_FROM_CLOCK  Default Brownian RNG seed from local wall-clock time
%
%  Description:
%    Returns a nonnegative integer suitable for rng(..., 'twister') from
%    the current local time formatted as MMddHHmmss (month, day, hour,
%    minute, second), e.g. 0410210635 for April 10 at 21:06:35.
%
%  Syntax:
%    seed = brown_rng_seed_from_clock()
%
%  Inputs:
%    (none)
%
%  Outputs:
%    seed - Nonnegative integer derived from wall-clock time (MMddHHmmss)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-10

seed = str2double(char(datetime('now', 'Format', 'MMddHHmmss')));

end
