function u_bg = coupled_background_flow(X_all, geom, params)
%COUPLED_BACKGROUND_FLOW  Evaluate steady background flow at all hydro points
%
%  Description:
%    Computes the dimensionless background flow velocity at every
%    hydrodynamic evaluation point (body surface + pili nodes).
%    Returns a vector compatible with the DOF layout used by
%    build_coupled_lin_sys (translation block followed by twist block).
%
%    Supported flow profiles:
%      'off'      - Zero flow everywhere (default)
%      'uniform'  - Spatially uniform velocity
%      'shear_z'  - Piecewise-linear shear along a specified axis
%
%  Syntax:
%    u_bg = coupled_background_flow(X_all, geom, params)
%
%  Inputs:
%    X_all  - [N_total x 3] Lab-frame positions of all points
%    geom   - Geometry struct from reconstruct_geometry
%    params - Parameter structure (from init_coupled_params)
%
%  Outputs:
%    u_bg - [(3*N_total + sum(N_segs)) x 1] Dimensionless background
%           velocity.  Layout:
%             u_bg(1 : 3*N_total)  - translational velocity (xyz per point)
%             u_bg(3*N_total+1 : end) - twist rates (dot(d3, Omega) per segment)
%
%  References:
%    [1] Walker et al. (2020), Physical Review E.
%    [2] 3D-Flexible-Filaments/background_flow.m (single-filament version)
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-04-03

N_total = size(X_all, 1);
N_segs  = params.N_segments;
M       = params.M_pili;
n_twist = sum(N_segs);
n_total_dof = 3 * N_total + n_twist;

bg = params.bg_flow;

if strcmp(bg.type, 'off')
    u_bg = zeros(n_total_dof, 1);
    return
end

U_max  = bg.U_max_dimless;
d_flow = bg.direction(:)';

switch bg.type
    case 'uniform'
        u_vec = U_max * d_flow;
        u_trans = repmat(u_vec, N_total, 1);
        omega_lab = [0, 0, 0];

    case 'shear_z'
        s_axis = bg.shear_axis(:)';
        H      = bg.shear_height_dimless;

        h_proj = X_all * s_axis';
        profile = U_max * max(0, min(h_proj / H, 1));

        u_trans = profile .* d_flow;

        omega_lab_dir = cross(s_axis, d_flow);
        omega_lab_mag = U_max / H;
        omega_lab = omega_lab_mag * omega_lab_dir;

    otherwise
        error('coupled_background_flow:UnknownType', ...
            'Unknown background_flow type: ''%s''.', bg.type);
end

u_bg = zeros(n_total_dof, 1);

for i = 1:N_total
    u_bg(3*(i-1)+1 : 3*i) = u_trans(i, :)';
end

twist_offset = 3 * N_total;
seg_count    = 0;
pt_offset    = params.N_body;

for k = 1:M
    N_k   = N_segs(k);
    d3_k  = geom.d3{k};

    if strcmp(bg.type, 'shear_z')
        s_axis = bg.shear_axis(:)';
        X_k    = geom.X_pili{k};
        Z_seg  = 0.5 * (X_k(1:N_k,:) + X_k(2:N_k+1,:)) * s_axis';
        H      = bg.shear_height_dimless;
        inside = (Z_seg > 0) & (Z_seg < H);

        for j = 1:N_k
            if inside(j)
                Omega_half = 0.5 * omega_lab;
            else
                Omega_half = [0, 0, 0];
            end
            u_bg(twist_offset + seg_count + j) = dot(d3_k(j,:), Omega_half);
        end
    else
        for j = 1:N_k
            u_bg(twist_offset + seg_count + j) = dot(d3_k(j,:), 0.5 * omega_lab);
        end
    end

    seg_count  = seg_count + N_k;
    pt_offset  = pt_offset + N_k + 1;
end

end
