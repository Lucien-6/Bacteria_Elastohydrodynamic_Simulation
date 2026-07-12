function Vtool = sphere_particle_sample_cached(SPN)
%SPHERE_PARTICLE_SAMPLE_CACHED  Cached ParticleSampleSphere for SPN in (100,1000)
%
%  Syntax:
%    Vtool = sphere_particle_sample_cached(SPN)
%
%  Input:
%    SPN - Integer number of sphere surface points (expected: 100 < SPN < 1000)
%
%  Output:
%    Vtool - [SPN x 3] Unit-sphere points in toolbox z-axis frame
%            (before tool_to_body_frame rotation)
%
%  Returns the N-by-3 Vtool matrix (toolbox z-axis frame, before
%  tool_to_body_frame) either from disk cache or by running the
%  Riesz s-energy optimisation and saving the result for future use.
%
%  Cache location : <project_root>/src/geometry/cache/sphere_particle/
%  Cache filename : spherical_points_NXXXX.mat
%
%  Each cache file stores:
%    Vtool  — N-by-3 double, unit-sphere points in toolbox frame
%    meta   — struct with fields used to validate cache consistency:
%               meta.N          point count
%               meta.Etol       energy tolerance (toolbox default)
%               meta.Dtol_deg   displacement tolerance in degrees (toolbox default)
%               meta.Nitr       max iterations (toolbox default)
%               meta.Vo_method  initial-position strategy identifier
%               meta.format     file format version (bump on any constant change)
%               meta.created    creation timestamp (informational)
%
%  If the stored meta does not match the constants defined below the file is
%  treated as a miss: the optimisation runs again and the file is overwritten.
%
%  On any I/O failure (unreadable file, full disk, read-only path …) the
%  function falls back to direct computation and the simulation continues
%  uninterrupted.
%
%  Called by capsule_sphere_cap_points when 100 < SPN < 1000.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Date:    2026-04-18
%  Updated: 2026-07-12

% --- Cache key constants — bump FORMAT_VERSION on any change ---------------
CACHE_ETOL      = 1e-5;          % ParticleSampleSphere toolbox default
CACHE_DTOL_DEG  = 1e-4;          % ParticleSampleSphere toolbox default (degrees)
CACHE_NITR      = 1e4;           % ParticleSampleSphere toolbox default
CACHE_VO_METHOD = 'SpiralSampleSphere';  % deterministic initial-position strategy
FORMAT_VERSION  = 1;
% ---------------------------------------------------------------------------

cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache', 'sphere_particle');
fpath     = fullfile(cache_dir, sprintf('spherical_points_N%04d.mat', SPN));

% --- Attempt cache hit ---
if isfile(fpath)
    Vtool = try_load(fpath, SPN, ...
        CACHE_ETOL, CACHE_DTOL_DEG, CACHE_NITR, CACHE_VO_METHOD, FORMAT_VERSION);
    if ~isempty(Vtool)
        return
    end
end

% --- Cache miss or invalid: run optimisation ---
Vtool = run_particle(SPN);

% --- Attempt cache save (fail silently on any error) ---
try_save(fpath, cache_dir, SPN, Vtool, ...
    CACHE_ETOL, CACHE_DTOL_DEG, CACHE_NITR, CACHE_VO_METHOD, FORMAT_VERSION);

end


% ==========================================================================
function Vtool = try_load(fpath, SPN, etol, dtol_deg, nitr, vo_method, fmt_ver)
%TRY_LOAD  Load and validate a cache file; return [] on any failure.
%  Deletes the file on any failure so a stale or corrupt entry does not
%  persist across runs if the subsequent try_save also fails.

Vtool = [];

try
    s = load(fpath, 'Vtool', 'meta');
catch
    warning('sphere_particle_sample_cached:loadFail', ...
        'Cache file unreadable, will recompute: %s', fpath);
    try; delete(fpath); catch; end
    return
end

% Validate structure and meta in one pass; delete and bail on any mismatch.
valid = false;
if isfield(s, 'Vtool') && isfield(s, 'meta')
    V = s.Vtool;
    m = s.meta;
    required_fields = {'N', 'Etol', 'Dtol_deg', 'Nitr', 'Vo_method', 'format'};
    valid = isnumeric(V) && ismatrix(V) && ...
            size(V, 1) == SPN  && size(V, 2) == 3 && ...
            all(abs(vecnorm(V, 2, 2) - 1) <= 1e-9) && ...
            isstruct(m) && all(isfield(m, required_fields)) && ...
            m.N == SPN && m.Etol == etol && m.Dtol_deg == dtol_deg && ...
            m.Nitr == nitr && strcmp(m.Vo_method, vo_method) && ...
            m.format == fmt_ver;
end

if ~valid
    try; delete(fpath); catch; end
    return
end

Vtool = V;

end


% ==========================================================================
function Vtool = run_particle(SPN)
%RUN_PARTICLE  Generate sphere points via ParticleSampleSphere (no cache).

Vo    = SpiralSampleSphere(SPN, false);
Vtool = ParticleSampleSphere('Vo', Vo, 'qdlg', false, 'upd', false);

end


% ==========================================================================
function try_save(fpath, cache_dir, SPN, Vtool, etol, dtol_deg, nitr, vo_method, fmt_ver)
%TRY_SAVE  Write Vtool and meta to cache file; fail silently on any error.

tmp_path = '';
try
    if ~isfolder(cache_dir)
        mkdir(cache_dir);
    end

    meta = struct( ...
        'N',          SPN, ...
        'Etol',       etol, ...
        'Dtol_deg',   dtol_deg, ...
        'Nitr',       nitr, ...
        'Vo_method',  vo_method, ...
        'format',     fmt_ver, ...
        'created',    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')) ...
    );

    tmp_path = sprintf('%s.tmp_%d', fpath, round(rand() * 1e9));
    save(tmp_path, 'Vtool', 'meta', '-v7');
    movefile(tmp_path, fpath);

catch err
    warning('sphere_particle_sample_cached:saveFail', ...
        'Cache save failed (%s), continuing without cache.', err.message);
    if ~isempty(tmp_path) && isfile(tmp_path)
        try; delete(tmp_path); catch; end
    end
end

end
