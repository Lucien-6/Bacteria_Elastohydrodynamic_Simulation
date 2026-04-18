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
%    3) Else if SPN >= 14, ParticleSampleSphere with the specification below.
%    4) Else SpiralSampleSphere(SPN).
%
%  --- ParticleSampleSphere specification (name-value overrides) ---
%  These values replace toolbox defaults for silent batch use. Change only
%  deliberately; document revisions for reproducibility.
%    'Etol'   Absolute energy convergence tolerance ............ 1e-4
%    'Dtol'   Max particle displacement tolerance .............. 1e-3 deg
%    'Nitr'   Maximum iterations ................................... 3500
%    'qdlg'   false — no confirmation dialog for large N
%    'upd'    false — suppress iteration print to command window
%  Riesz s-energy minimisation parameters follow ParticleSampleSphere.m (s=1).
%
%  Toolbox z-axis pole -> body x-axis pole (capsule).
%
%  Requires: Anton Semechko S2 Sampling Toolbox on the MATLAB path.
%
%  Author: Lucien (lucien-6@qq.com)
%  Date:   2026-04-17

% --- ParticleSampleSphere public tuning (keep in sync with header block) ---
SPHERE_PARTICLE_ETOL = 1e-4;
SPHERE_PARTICLE_DTOL_DEG = 1e-3;
SPHERE_PARTICLE_NITR = 3500;

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
    particle_opts = {'qdlg', false, 'upd', false, ...
        'Etol', SPHERE_PARTICLE_ETOL, 'Dtol', SPHERE_PARTICLE_DTOL_DEG, ...
        'Nitr', SPHERE_PARTICLE_NITR};
    Vtool = ParticleSampleSphere('N', SPN, particle_opts{:});
else
    Vtool = SpiralSampleSphere(SPN, false);
end

V = tool_to_body_frame(Vtool);

end

function V = tool_to_body_frame(Vtool)
    V = [Vtool(:, 3), Vtool(:, 2), Vtool(:, 1)];
    V = V ./ vecnorm(V, 2, 2);
end
