function validate_simulation_params(cfg)
%VALIDATE_SIMULATION_PARAMS  Pre-flight parameter sanity checks
%
%  Description:
%    Performs comprehensive validation of user-configurable simulation
%    parameters before the (potentially long) initialisation and ODE
%    integration.  Raises errors for critical violations and warnings
%    for suboptimal values.  All errors are aggregated and reported
%    together so the user can fix them in one pass.
%
%  Syntax:
%    validate_simulation_params(cfg)
%
%  Input:
%    cfg - struct with fields matching user-configurable parameters:
%      T_obs, dt, minor_axis, major_axis, N_body, seg_per_um,
%      Pili_Matrix, E, r_pili, mu, rho_bac, rho_fluid,
%      AbsTol, RelTol, tlim, enable_brown, dt_B,
%      animation_fps, sigma_p, fixman_midpoint
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17
%  Updated: 2026-03-17 — v3.2: Fixed duplicate site detection indexing bug
%  Updated: 2026-03-30 — allow major_axis == minor_axis (spherical body)

errs  = {};
warns = {};

%% ---- Critical errors ----

% Simulation time
if cfg.T_obs <= 0
    errs{end+1} = 'T_obs must be > 0.';
end
if cfg.dt <= 0
    errs{end+1} = 'dt (output time step) must be > 0.';
elseif cfg.dt > cfg.T_obs
    errs{end+1} = sprintf('dt = %g > T_obs = %g: output step exceeds total time.', ...
        cfg.dt, cfg.T_obs);
end

% Body geometry
if cfg.minor_axis <= 0
    errs{end+1} = 'minor_axis must be > 0.';
end
if cfg.major_axis <= 0
    errs{end+1} = 'major_axis must be > 0.';
elseif cfg.major_axis < cfg.minor_axis
    errs{end+1} = sprintf('major_axis (%.2e) must be >= minor_axis (%.2e).', ...
        cfg.major_axis, cfg.minor_axis);
end
if cfg.N_body < 10
    errs{end+1} = sprintf('N_body = %d is too small for capsule geometry (minimum 10).', ...
        cfg.N_body);
end

% Material
if cfg.E <= 0
    errs{end+1} = 'Young''s modulus E must be > 0.';
end
if cfg.r_pili <= 0
    errs{end+1} = 'Pili radius r_pili must be > 0.';
end
if cfg.sigma_p < 0 || cfg.sigma_p >= 0.5
    errs{end+1} = sprintf('Poisson ratio sigma_p = %g outside valid range [0, 0.5).', ...
        cfg.sigma_p);
end

% Fluid
if cfg.mu <= 0
    errs{end+1} = 'Dynamic viscosity mu must be > 0.';
end

% Numerics
if cfg.tlim <= 0
    errs{end+1} = 'ODE timeout tlim must be > 0.';
end
if cfg.animation_fps <= 0
    errs{end+1} = 'animation_fps must be > 0.';
end

% Pili matrix (only if non-empty)
PM = cfg.Pili_Matrix;
if ~isempty(PM)
    if size(PM, 1) ~= 8
        errs{end+1} = sprintf('Pili_Matrix must have exactly 8 rows, got %d.', size(PM, 1));
    else
        for kk = 1:size(PM, 2)
            if PM(2, kk) <= 0
                errs{end+1} = sprintf('Pili %d: length (row 2) must be > 0, got %g.', ...
                    kk, PM(2, kk));
            end
            if PM(8, kk) < 1 || PM(8, kk) > 14 || floor(PM(8, kk)) ~= PM(8, kk)
                errs{end+1} = sprintf('Pili %d: attachment site (row 8) must be integer 1-14, got %g.', ...
                    kk, PM(8, kk));
            end
            if PM(1, kk) < 1 || PM(1, kk) > 5 || floor(PM(1, kk)) ~= PM(1, kk)
                errs{end+1} = sprintf('Pili %d: type (row 1) must be integer 1-5, got %g.', ...
                    kk, PM(1, kk));
            end
        end
    end
end

% Brownian
if cfg.enable_brown
    if cfg.dt_B <= 0
        errs{end+1} = 'Brownian time step dt_B must be > 0 when enable_brown = true.';
    end
end

% Fixman midpoint
if isfield(cfg, 'fixman_midpoint') && cfg.fixman_midpoint && ~cfg.enable_brown
    warns{end+1} = 'fixman_midpoint = true has no effect when enable_brown = false.';
end

%% ---- Warnings ----

if cfg.N_body >= 10 && cfg.N_body < 100
    warns{end+1} = sprintf('N_body = %d < 100: hydrodynamic accuracy may be poor.', cfg.N_body);
end

if cfg.seg_per_um < 5
    warns{end+1} = sprintf('seg_per_um = %d < 5: RSS accuracy may be severely degraded.', ...
        cfg.seg_per_um);
end

if cfg.rho_bac <= cfg.rho_fluid
    warns{end+1} = sprintf('rho_bac (%.1f) <= rho_fluid (%.3f): bacterium will not sink.', ...
        cfg.rho_bac, cfg.rho_fluid);
end

if cfg.AbsTol <= 0 || cfg.AbsTol > 1
    warns{end+1} = sprintf('AbsTol = %g outside typical range (0, 1].', cfg.AbsTol);
end
if cfg.RelTol <= 0 || cfg.RelTol > 1
    warns{end+1} = sprintf('RelTol = %g outside typical range (0, 1].', cfg.RelTol);
end

if cfg.enable_brown && cfg.dt_B > 0 && cfg.dt_B >= cfg.dt
    warns{end+1} = sprintf('dt_B = %g >= dt = %g: Brownian step >= output step.', ...
        cfg.dt_B, cfg.dt);
end

% Pili duplicate attachment sites
if ~isempty(PM) && size(PM, 1) >= 8 && size(PM, 2) > 1
    sites = PM(8, :);
    [u_sites, ~, ic] = unique(sites);
    counts = accumarray(ic, 1);
    if any(counts > 1)
        dup_sites = u_sites(counts > 1);
        warns{end+1} = sprintf('Duplicate pili attachment sites detected: %s.', ...
            mat2str(dup_sites));
    end
end

%% ---- Report ----

for i = 1:numel(warns)
    warning('validate_params:SuboptimalValue', '[VALIDATE] WARNING: %s', warns{i});
end

if ~isempty(errs)
    msg = strjoin(cellfun(@(s) ['  - ' s], errs, 'UniformOutput', false), newline);
    error('validate_params:InvalidParams', ...
        'Parameter validation failed with %d error(s):\n%s', numel(errs), msg);
end

end
