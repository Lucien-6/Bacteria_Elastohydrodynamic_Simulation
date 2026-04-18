function cache = capsule_sphere_mesh_cache()
%CAPSULE_SPHERE_MESH_CACHE  Cached Ico/Quad subdivided meshes from S2 toolbox
%
%  Vertex counts Nico(k+1), Nquad(k+1) for k = 0..4 are taken from
%  size(mesh points) after SubdivideSphericalMesh(...) so they match the
%  installed toolbox (no closed-form formula).
%
%  Author: Lucien (lucien-6@qq.com)
%  Date:   2026-04-17

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
        error('capsule_sphere_mesh_cache:unsupportedMeshOutput %s', class(TR));
    end
end
