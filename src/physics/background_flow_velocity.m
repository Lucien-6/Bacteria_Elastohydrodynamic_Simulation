function [u_trans, omega_lab] = background_flow_velocity(X, params)
%BACKGROUND_FLOW_VELOCITY  Translational background flow at arbitrary points
%
%  Description:
%    Evaluates the dimensionless steady background flow velocity at an
%    arbitrary set of lab-frame points.  This is the single source of truth
%    for the translational flow profile; coupled_background_flow (RHS
%    assembly) and compute_plot_flowfield (post-processing) both call it so
%    the simulated and the visualised ambient flow can never diverge.
%
%    Supported flow profiles (params.bg_flow.type):
%      'off'      - Zero flow everywhere
%      'uniform'  - Spatially uniform velocity  U_max * direction
%      'shear_z'  - Piecewise-linear Couette layer
%                   u(z) = U_max * max(0, min(z/H, 1)) * direction,
%                   z = X . shear_axis  (see coupled_background_flow for the
%                   vorticity conventions of the plateau regions)
%
%  Syntax:
%    u_trans = background_flow_velocity(X, params)
%    [u_trans, omega_lab] = background_flow_velocity(X, params)
%
%  Inputs:
%    X      - [N x 3] Lab-frame positions (dimensionless, units of L_0)
%    params - Parameter structure (from init_coupled_params); only
%             params.bg_flow is used
%
%  Outputs:
%    u_trans   - [N x 3] Dimensionless background velocity at each point
%    omega_lab - [1 x 3] Lab-frame fluid angular velocity (half vorticity)
%                of the linear region; zero for 'off' and 'uniform'.  For
%                'shear_z' it applies only inside 0 < z < H; callers must
%                apply the layer mask themselves.
%
%  See also: coupled_background_flow, compute_plot_flowfield
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-09-04

N = size(X, 1);
bg = params.bg_flow;

if ~isfield(bg, 'type') || strcmp(bg.type, 'off')
    u_trans   = zeros(N, 3);
    omega_lab = [0, 0, 0];
    return
end

U_max  = bg.U_max_dimless;
d_flow = bg.direction(:)';

switch bg.type
    case 'uniform'
        u_trans   = repmat(U_max * d_flow, N, 1);
        omega_lab = [0, 0, 0];

    case 'shear_z'
        s_axis = bg.shear_axis(:)';
        H      = bg.shear_height_dimless;

        h_proj  = X * s_axis';
        profile = U_max * max(0, min(h_proj / H, 1));
        u_trans = profile .* d_flow;

        % Fluid angular velocity = 0.5 * vorticity = 0.5 * (U_max/H) (s x d)
        omega_lab = 0.5 * (U_max / H) * cross(s_axis, d_flow);

    otherwise
        error('background_flow_velocity:UnknownType', ...
            'Unknown background_flow type: ''%s''.', bg.type);
end

end
