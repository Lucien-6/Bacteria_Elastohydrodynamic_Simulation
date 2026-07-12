function SPN = snap_sphere_cap_count(SPN0)
%SNAP_SPHERE_CAP_COUNT  Snap combined-cap full-sphere count to native mesh N
%
%  Syntax:
%    SPN = snap_sphere_cap_count(SPN0)
%
%  Input:
%    SPN0 - Requested full-sphere sample count for the two combined capsule
%           caps (any positive integer)
%
%  Output:
%    SPN  - Snapped combined-cap full-sphere count: if |SPN0 - N| <= 5 for
%           some Ico/Quad mesh at subdivision level k in 0..4, returns that N
%           (smallest |SPN0-N|; ties: Ico over Quad, then lower k); otherwise
%           returns SPN0 unchanged.
%
%  Native vertex counts come from capsule_sphere_mesh_cache
%  (actual SubdivideSphericalMesh sizes for the installed toolbox).
%
%  Author: Lucien (lucien-6@qq.com)
%  Date:   2026-04-17
%  Updated: 2026-07-12

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
