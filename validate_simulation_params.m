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
%      AbsTol, RelTol, tlim, enable_brown, dt_B, Temper,
%      animation_fps, sigma_p, fixman_midpoint
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17

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

%  Fixman recommended when Brownian is on AND pili are present.
%  The simple Euler-Maruyama scheme drops the kB*T * div(M_eff) Ito drift,
%  which becomes non-negligible when pili conformation strongly modulates
%  the effective body mobility.  See Brownian_Implementation_Plan.md
%  section 2.4.1 for a quantitative discussion and a decision matrix.
if cfg.enable_brown && isfield(cfg, 'fixman_midpoint') && ~cfg.fixman_midpoint ...
        && ~isempty(cfg.Pili_Matrix)
    warns{end+1} = ['Brownian enabled with pili present but fixman_midpoint = false: ', ...
        'the kB*T * div(M_eff) Ito drift term is omitted (simple Euler-Maruyama). ', ...
        'For long-time diffusion statistics or large-amplitude pili motion, ', ...
        'consider enabling fixman_midpoint (see Brownian_Implementation_Plan.md section 2.4.1).'];
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

%  dt_B vs. Stokes-Einstein relaxation times (back-of-envelope spherical
%  estimate; a factor-of-2 accuracy is sufficient for a warning).  We use
%  the major semi-axis as a conservative radius so the warning triggers
%  only when dt_B is clearly too coarse for accurate Brownian statistics.
if cfg.enable_brown && cfg.dt_B > 0 && isfield(cfg, 'Temper') ...
        && isfinite(cfg.Temper) && cfg.Temper > 0 && cfg.mu > 0 ...
        && cfg.major_axis > 0 && cfg.minor_axis > 0
    kB    = 1.380649e-23;
    a_eff = max(cfg.major_axis, cfg.minor_axis);
    D_T   = kB * cfg.Temper / (6 * pi * cfg.mu * a_eff);
    D_R   = kB * cfg.Temper / (8 * pi * cfg.mu * a_eff^3);
    tau_trans = a_eff^2 / D_T;
    tau_rot   = 1 / (6 * D_R);
    tau_min   = min(tau_trans, tau_rot);
    if isfinite(tau_min) && tau_min > 0 && cfg.dt_B > 0.1 * tau_min
        warns{end+1} = sprintf(['dt_B = %.3g s exceeds 10%% of estimated ', ...
            'min(tau_trans, tau_rot) ~ %.3g s (Stokes-Einstein with a = %.2g m). ', ...
            'Consider dt_B <= %.3g s for accurate Brownian statistics.'], ...
            cfg.dt_B, tau_min, a_eff, 0.1 * tau_min);
    end
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

% Background flow
if isfield(cfg, 'background_flow') && isstruct(cfg.background_flow)
    bg = cfg.background_flow;
    valid_bg = {'off', 'uniform', 'shear_z'};
    if ~any(strcmp(bg.type, valid_bg))
        errs{end+1} = sprintf('background_flow.type must be one of {off, uniform, shear_z}, got ''%s''.', bg.type);
    end
    if ~strcmp(bg.type, 'off')
        if ~isfield(bg, 'U_max') || bg.U_max <= 0
            errs{end+1} = 'background_flow.U_max must be > 0 when flow is enabled.';
        end
        if isfield(bg, 'direction') && numel(bg.direction) ~= 3
            errs{end+1} = 'background_flow.direction must be a 3-element vector.';
        end
        if isfield(bg, 'shear_axis') && numel(bg.shear_axis) ~= 3
            errs{end+1} = 'background_flow.shear_axis must be a 3-element vector.';
        end
        if strcmp(bg.type, 'shear_z')
            if ~isfield(bg, 'shear_height_m') || bg.shear_height_m <= 0
                errs{end+1} = 'background_flow.shear_height_m must be > 0 for shear_z flow.';
            end
        end
    end
    if ~strcmp(bg.type, 'off') && cfg.enable_brown
        warns{end+1} = ['Background flow + Brownian: flow is deterministic; ' ...
            'diffusion tensor computed without flow contribution.'];
    end
end

% Per-pili: intrinsic curvature and active moment
M_pili = size(PM, 2);
if isfield(cfg, 'pili_ext') && M_pili > 0
    valid_ic = {'frozen', 'zero', 'constant', 'helix'};
    valid_am = {'off', 'traveling_wave', 'custom'};
    for kk = 1:min(numel(cfg.pili_ext), M_pili)
        pe = cfg.pili_ext(kk);
        % Intrinsic curvature
        if isfield(pe, 'intrinsic_curvature') && ~isempty(pe.intrinsic_curvature)
            ic_val = pe.intrinsic_curvature;
            if ischar(ic_val) || isstring(ic_val)
                ic_type = char(ic_val);
            elseif isstruct(ic_val) && isfield(ic_val, 'type')
                ic_type = ic_val.type;
            else
                ic_type = '';
            end
            if ~any(strcmp(ic_type, valid_ic))
                errs{end+1} = sprintf('Pili %d: intrinsic_curvature type must be one of {frozen, zero, constant, helix}, got ''%s''.', kk, ic_type);
            end
            if strcmp(ic_type, 'constant')
                if isstruct(ic_val) && isfield(ic_val, 'kappa')
                    if numel(ic_val.kappa) ~= 3
                        errs{end+1} = sprintf('Pili %d: intrinsic_curvature.kappa must be 3-element, got %d.', kk, numel(ic_val.kappa));
                    end
                else
                    errs{end+1} = sprintf('Pili %d: intrinsic_curvature ''constant'' requires .kappa field.', kk);
                end
            end
        end
        % Active moment
        if isfield(pe, 'active_moment') && isstruct(pe.active_moment)
            am_val = pe.active_moment;
            am_type = '';
            if isfield(am_val, 'type'), am_type = am_val.type; end
            if ~any(strcmp(am_type, valid_am))
                errs{end+1} = sprintf('Pili %d: active_moment type must be one of {off, traveling_wave, custom}, got ''%s''.', kk, am_type);
            end
            if strcmp(am_type, 'traveling_wave')
                if ~isfield(am_val, 'amplitude_Nm') || am_val.amplitude_Nm < 0
                    errs{end+1} = sprintf('Pili %d: active_moment.amplitude_Nm must be >= 0.', kk);
                end
                if ~isfield(am_val, 'omega_Hz') || am_val.omega_Hz <= 0
                    errs{end+1} = sprintf('Pili %d: active_moment.omega_Hz must be > 0.', kk);
                end
                if isfield(am_val, 'wavenumber') && am_val.wavenumber < 0
                    warns{end+1} = sprintf(['Pili %d: active_moment.wavenumber = %.4g is negative. ' ...
                        'Wave propagates in -s direction (from tip to root).'], kk, am_val.wavenumber);
                end
                if isfield(am_val, 'd3_amp') && am_val.d3_amp ~= 0
                    warns{end+1} = sprintf(['Pili %d: active_moment.d3_amp = %.4g (non-zero). ' ...
                        'This drives torsion rather than bending, which is unusual for standard pili beating.'], ...
                        kk, am_val.d3_amp);
                end
                if isfield(am_val, 'phase_rad')
                    if ~isnumeric(am_val.phase_rad) || ~isscalar(am_val.phase_rad) || ~isfinite(am_val.phase_rad)
                        errs{end+1} = sprintf('Pili %d: active_moment.phase_rad must be a finite scalar.', kk);
                    end
                end
            end
        end
    end
end

if M_pili == 0 && isfield(cfg, 'pili_ext') && ~isempty(cfg.pili_ext)
    warns{end+1} = 'pili_ext is non-empty but no pili are defined (Pili_Matrix is empty).';
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
