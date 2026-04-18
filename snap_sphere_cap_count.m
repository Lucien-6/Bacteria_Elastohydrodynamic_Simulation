function SPN = snap_sphere_cap_count(SPN0)
%SNAP_SPHERE_CAP_COUNT  Snap cap count to Ico/Quad native N when close enough
%
%  If |SPN0 - N| <= 5 for some Ico or Quad mesh at subdivision level k in
%  0..4, set SPN to that N (smallest |SPN0-N|; ties: Ico over Quad, then
%  lower k). Native vertex counts come from capsule_sphere_mesh_cache
%  (actual SubdivideSphericalMesh sizes for the installed toolbox).
%
%  Author: Lucien (lucien-6@qq.com)
%  Date:   2026-04-17

tol = 5;
cache = capsule_sphere_mesh_cache();
candidates = [];

for k = 0:4
    Ni = cache.Nico(k + 1);
    if abs(SPN0 - Ni) <= tol
        candidates(end + 1, :) = [abs(SPN0 - Ni), 0, k, Ni]; %#ok<AGROW>
    end
end
for k = 0:4
    Nq = cache.Nquad(k + 1);
    if abs(SPN0 - Nq) <= tol
        candidates(end + 1, :) = [abs(SPN0 - Nq), 1, k, Nq]; %#ok<AGROW>
    end
end

if isempty(candidates)
    SPN = SPN0;
else
    [~, ord] = sortrows(candidates, [1, 2, 3]);
    SPN = candidates(ord(1), 4);
end

end
