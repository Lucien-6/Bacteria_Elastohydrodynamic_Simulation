function V = capsule_sphere_cap_points(SPN)
%CAPSULE_SPHERE_CAP_POINTS  Unit-sphere cap samples for capsule body (S2 toolbox)
%
%  SPN must already equal the intended cap count (call snap_sphere_cap_count
%  in capsule_body_building when adjusting allocation to Ico/Quad native N).
%
%  Decision order:
%    1) If SPN exactly matches Ico or Quad native vertex count for some k in
%       0..4 (counts from capsule_sphere_mesh_cache), use all vertices of
%       that mesh (no subsampling). If Ico and Quad share the same N, Ico
%       wins, then lower k.
%    2) Else if SPN >= 1000, SpiralSampleSphere(SPN).
%    3) Else if SPN >= 14, ParticleSampleSphere (see below).
%    4) Else SpiralSampleSphere(SPN).
%
%  --- ParticleSampleSphere (name-value overrides only where needed) ---
%  Initial positions: 'Vo' = SpiralSampleSphere(SPN,false) (deterministic, no RNG).
%  'Etol', 'Dtol', 'Nitr' are not passed — toolbox defaults apply
%  (Etol=1e-5, Dtol=1e-4 deg, Nitr=1e4 per ParticleSampleSphere/VerifyInputArgs).
%    'qdlg'   false — no confirmation dialog for large N
%    'upd'    false — suppress iteration print to command window
%  'asym' is omitted (false). Riesz s-energy uses s=1 per toolbox default.
%
%  --- Particle result cache ---
%  When 100 < SPN < 1000 the result is read from / written to disk via
%  sphere_particle_sample_cached (cache/sphere_particle/spherical_points_NXXXX.mat).
%  For 14 <= SPN <= 100 ParticleSampleSphere is called directly (no cache).
%
%  Toolbox z-axis pole -> body x-axis pole (capsule).
%
%  Requires: Anton Semechko S2 Sampling Toolbox on the MATLAB path.
%
%  Author: Lucien (lucien-6@qq.com)
%  Date:   2026-04-17
%  Updated: 2026-04-18 — ParticleSampleSphere: toolbox Etol/Dtol/Nitr defaults;
%                        deterministic Vo via SpiralSampleSphere(SPN,false);
%                        disk cache for 100 < SPN < 1000 via sphere_particle_sample_cached

if SPN < 1
    error('capsule_sphere_cap_points: SPN must be at least 1.');
end

cache = capsule_sphere_mesh_cache();

candidates = [];
for k = 0:4
    Ni = cache.Nico(k + 1);
    if SPN == Ni
        candidates(end + 1, :) = [0, k]; %#ok<AGROW>
    end
end
for k = 0:4
    Nq = cache.Nquad(k + 1);
    if SPN == Nq
        candidates(end + 1, :) = [1, k]; %#ok<AGROW>
    end
end

if ~isempty(candidates)
    [~, ord] = sortrows(candidates, [1, 2]);
    sel = candidates(ord(1), :);
    k = sel(2);
    if sel(1) == 0
        Vtool = cache.ico{k + 1};
    else
        Vtool = cache.quad{k + 1};
    end
    V = tool_to_body_frame(Vtool);
    return
end

if SPN >= 1000
    Vtool = SpiralSampleSphere(SPN, false);
elseif SPN >= 14
    if SPN > 100 && SPN < 1000
        Vtool = sphere_particle_sample_cached(SPN);
    else
        Vo    = SpiralSampleSphere(SPN, false);
        Vtool = ParticleSampleSphere('Vo', Vo, 'qdlg', false, 'upd', false);
    end
else
    Vtool = SpiralSampleSphere(SPN, false);
end

V = tool_to_body_frame(Vtool);

end

function V = tool_to_body_frame(Vtool)
    V = [Vtool(:, 3), Vtool(:, 2), Vtool(:, 1)];
    V = V ./ vecnorm(V, 2, 2);
end
