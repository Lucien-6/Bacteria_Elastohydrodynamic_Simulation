function cache = capsule_sphere_mesh_cache()
%CAPSULE_SPHERE_MESH_CACHE  Cached Ico/Quad subdivided meshes from S2 toolbox
%
%  Description:
%    Pre-computes and caches IcosahedronMesh and QuadCubeMesh subdivisions
%    (levels 0-4) from the S2 toolbox using SubdivideSphericalMesh.  The
%    results are stored in persistent variables so that subsequent calls
%    return instantly without re-running the toolbox.  Vertex counts
%    Nico(k+1), Nquad(k+1) for k = 0..4 are taken from the actual mesh
%    sizes after subdivision so they always match the installed toolbox
%    (no closed-form formula is used).
%
%  Syntax:
%    cache = capsule_sphere_mesh_cache()
%
%  Output:
%    cache - Struct with fields:
%      .Nico   - [1 x 5] Vertex counts for Ico subdivisions k=0..4
%      .Nquad  - [1 x 5] Vertex counts for Quad subdivisions k=0..4
%      .ico    - {1 x 5} Cell array of [N x 3] Ico vertex matrices
%      .quad   - {1 x 5} Cell array of [N x 3] Quad vertex matrices
%
%  Errors:
%    Throws capsule_sphere_mesh_cache:UnsupportedMeshOutput if the installed
%    S2 toolbox returns neither a triangulation nor a supported mesh struct.
%
%  See also: capsule_sphere_cap_points, snap_sphere_cap_count,
%            sphere_particle_sample_cached
%
%  Author: Lucien (lucien-6@qq.com)
%  Created: 2026-04-17
%  Updated: 2026-07-12

persistent icoPool quadPool Nico Nquad ready
if isempty(ready)
    icoPool = cell(1, 5);
    quadPool = cell(1, 5);
    Nico = zeros(1, 5);
    Nquad = zeros(1, 5);
    for k = 0:4
        TR = SubdivideSphericalMesh(IcosahedronMesh, k);
        P = mesh_vertices(TR);
        icoPool{k + 1} = P;
        Nico(k + 1) = size(P, 1);
        fv = SubdivideSphericalMesh(QuadCubeMesh, k);
        Q = mesh_vertices(fv);
        quadPool{k + 1} = Q;
        Nquad(k + 1) = size(Q, 1);
    end
    ready = true;
end

cache = struct();
cache.Nico = Nico;
cache.Nquad = Nquad;
cache.ico = icoPool;
cache.quad = quadPool;

end

function P = mesh_vertices(TR)
    if isa(TR, 'triangulation')
        P = TR.Points;
    elseif isstruct(TR) && isfield(TR, 'vertices')
        P = TR.vertices;
    elseif isstruct(TR) && isfield(TR, 'Points')
        P = TR.Points;
    else
        error('capsule_sphere_mesh_cache:UnsupportedMeshOutput', ...
            'Unsupported spherical mesh output type: %s.', class(TR));
    end
end
