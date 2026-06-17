function validate_simulation_params(cfg, logf_h)
%VALIDATE_SIMULATION_PARAMS  Pre-flight parameter sanity checks
%
%  Description:
%    Performs comprehensive validation of user-configurable simulation
%    parameters before the (potentially long) initialisation and ODE
%    integration.  Raises errors for critical violations and warnings
%    for suboptimal values.  All errors are aggregated and reported
%    together so the user can fix them in one pass.
%
%    All field accesses are guarded by isfield(), so the function accepts
%    the full cfg struct from run_simulation / load_sim_config / batch
%    sweeps without requiring a pre-curated subset.  Missing optional
%    fields are silently skipped.
%
%  Syntax:
%    validate_simulation_params(cfg)
%    validate_simulation_params(cfg, logf_h)
%
%  Input:
%    cfg    - struct with fields matching user-configurable parameters:
%      T_obs, dt, minor_axis, major_axis, N_body, seg_per_um,
%      Pili_Matrix, pili_ext, E, r_pili, mu, rho_bac, rho_fluid,
%      AbsTol, RelTol, tlim, enable_brown, dt_B, Temper,
%      animation_fps, sigma_p, delta, alpha_0,
%      geometry_mode, custom_geometry, initial_orientation,
%      background_flow
%    logf_h - (Optional) Logger function handle for dual-channel output.
%             If omitted or empty, falls back to warning-only.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-03-17
%  Updated: 2026-06-17 — v7.9.0 parameter-validation deep refactor:
%                         A1/A2 — alpha_0/delta checks moved here from
%                         init_coupled_params; full-cfg acceptance (isfield-guarded);
%                         P1 — pili_ext count consistency; P2 — delta range;
%                         P4 — N_body upper bound; P5 — background-flow orthogonality;
%                         P7 — unified string/struct acceptance for
%                         intrinsic_curvature/active_moment/motor_torque;
%                         P9 — AbsTol/RelTol wording; Bugbot-fixed alpha_0
%                         custom-mode regression (empty initial_orientation).

if nargin < 2 || isempty(logf_h)
    logf_h = @(varargin) [];   % no-op fallback
end

errs  = {};
warns = {};

%% ---- Geometry mode dispatch ----
geometry_mode = 'parametric';
if isfield(cfg, 'geometry_mode') && ~isempty(cfg.geometry_mode)
    geometry_mode = cfg.geometry_mode;
end

if strcmp(geometry_mode, 'custom')
    %% ---- Custom geometry validation ----
    if ~isfield(cfg, 'custom_geometry') || isempty(cfg.custom_geometry)
        errs{end+1} = 'geometry_mode is ''custom'' but custom_geometry is not provided.';
    else
        [errs, warns] = validate_custom_geometry(cfg.custom_geometry, errs, warns);
    end

    % Enforce no-mixing: parametric-specific fields must NOT be set
    if isfield(cfg, 'minor_axis') && isnumeric(cfg.minor_axis) && ~isnan(cfg.minor_axis) && cfg.minor_axis > 0
        errs{end+1} = ['geometry_mode is ''custom'' but minor_axis is set. ', ...
                        'Parametric body parameters must not be used in custom mode.'];
    end
    if isfield(cfg, 'major_axis') && isnumeric(cfg.major_axis) && ~isnan(cfg.major_axis) && cfg.major_axis > 0
        errs{end+1} = ['geometry_mode is ''custom'' but major_axis is set. ', ...
                        'Parametric body parameters must not be used in custom mode.'];
    end
    if isfield(cfg, 'N_body') && isnumeric(cfg.N_body) && ~isnan(cfg.N_body) && cfg.N_body > 0
        errs{end+1} = ['geometry_mode is ''custom'' but N_body is set. ', ...
                        'Parametric body parameters must not be used in custom mode.'];
    end
    % v7.1: stricter check — any non-empty Pili_Matrix with actual columns is rejected
    if isfield(cfg, 'Pili_Matrix') && ~isempty(cfg.Pili_Matrix) && size(cfg.Pili_Matrix, 2) > 0
        errs{end+1} = ['geometry_mode is ''custom'' but Pili_Matrix is non-empty. ', ...
                        'Parametric pili parameters must not be used in custom mode.'];
    end

else
    %% ---- Parametric mode validation (existing checks) ----

    % initial_orientation is custom-mode-only
    if isfield(cfg, 'initial_orientation') && ~isempty(cfg.initial_orientation)
        errs{end+1} = ['initial_orientation is only supported in ', ...
                       'geometry_mode = ''custom''. In parametric mode, ', ...
                       'use alpha_0 instead.'];
    end

    % Enforce no-mixing: custom_geometry must NOT be set in parametric mode
    if isfield(cfg, 'custom_geometry') && ~isempty(cfg.custom_geometry)
        errs{end+1} = ['geometry_mode is ''parametric'' but custom_geometry is set. ', ...
                        'Custom geometry parameters must not be used in parametric mode.'];
    end

    % Body geometry (guarded field access)
    if isfield(cfg, 'minor_axis') && isnumeric(cfg.minor_axis) && ~isnan(cfg.minor_axis)
        if cfg.minor_axis <= 0
            errs{end+1} = 'minor_axis must be > 0.';
        end
    end
    if isfield(cfg, 'major_axis') && isnumeric(cfg.major_axis) && ~isnan(cfg.major_axis)
        if cfg.major_axis <= 0
            errs{end+1} = 'major_axis must be > 0.';
        elseif isfield(cfg, 'minor_axis') && isnumeric(cfg.minor_axis) ...
                && ~isnan(cfg.minor_axis) && cfg.minor_axis > 0 ...
                && cfg.major_axis < cfg.minor_axis
            errs{end+1} = sprintf('major_axis (%.2e) must be >= minor_axis (%.2e).', ...
                cfg.major_axis, cfg.minor_axis);
        end
    end
    if isfield(cfg, 'N_body') && isnumeric(cfg.N_body) && ~isnan(cfg.N_body)
        if cfg.N_body < 10
            errs{end+1} = sprintf('N_body = %d is too small for capsule geometry (minimum 10).', ...
                cfg.N_body);
        end
        % P4: O(N^2) distance matrix in init_coupled_params may exhaust memory
        if cfg.N_body > 2000
            warns{end+1} = sprintf(['N_body = %d > 2000: O(N^2) distance matrix in ', ...
                'init_coupled_params may exhaust memory (~%.0f MB for the body block).'], ...
                cfg.N_body, cfg.N_body^2 * 8 / 1e6);
        end
    end

    % Pili matrix (only if non-empty)
    PM = [];
    if isfield(cfg, 'Pili_Matrix')
        PM = cfg.Pili_Matrix;
    end
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

    % Per-pili: intrinsic curvature and active moment
    M_pili = size(PM, 2);
    if isfield(cfg, 'pili_ext') && M_pili > 0
        % P1: explicit count consistency — silent fallback in init_coupled_params
        % would otherwise mask user config errors.
        if numel(cfg.pili_ext) ~= M_pili
            errs{end+1} = sprintf(['pili_ext count (%d) != Pili_Matrix column count (%d). ', ...
                'Each pilus must have a matching pili_ext entry (use parse_pili_config ', ...
                'or provide exactly %d entries).'], numel(cfg.pili_ext), M_pili, M_pili);
        else
            % Build aux struct for stall-curve diagnostics (parametric mode):
            % pili lengths sit in row 2 of Pili_Matrix in microns.
            aux.pili_lengths_m = PM(2, 1:M_pili) * 1e-6;       % [1 x M_pili] in m
            aux.r_pili         = cfg.r_pili;
            aux.mu             = cfg.mu;
            [errs, warns] = validate_pili_ext(cfg.pili_ext, M_pili, errs, warns, aux);
        end
    end

    if M_pili == 0 && isfield(cfg, 'pili_ext') && ~isempty(cfg.pili_ext)
        warns{end+1} = 'pili_ext is non-empty but no pili are defined (Pili_Matrix is empty).';
    end
end

%% ---- Per-pili extended field validation for custom mode ----
if strcmp(geometry_mode, 'custom') ...
        && isfield(cfg, 'pili_ext') && ~isempty(cfg.pili_ext) ...
        && isfield(cfg, 'custom_geometry') && isstruct(cfg.custom_geometry) ...
        && isfield(cfg.custom_geometry, 'pili') && ~isempty(cfg.custom_geometry.pili)
    M_pili_custom_expected = numel(cfg.custom_geometry.pili);
    % P1: count consistency for custom mode
    if numel(cfg.pili_ext) ~= M_pili_custom_expected
        errs{end+1} = sprintf(['pili_ext count (%d) != custom_geometry.pili count (%d). ', ...
            'Each custom pilus must have a matching pili_ext entry.'], ...
            numel(cfg.pili_ext), M_pili_custom_expected);
    else
        % Custom mode: pili lengths require reading user files; skip stall
        % diagnostics in the validator.  init_coupled_params will still log
        % each pilus length, and the stall check could be re-enabled there.
        aux_custom = struct();    % empty -> no stall check inside validator
        [errs, warns] = validate_pili_ext(cfg.pili_ext, M_pili_custom_expected, errs, warns, aux_custom);
    end
end

%% ---- Common validations (both modes) ----

% Simulation time
if isfield(cfg, 'T_obs') && isnumeric(cfg.T_obs) && ~isempty(cfg.T_obs)
    if cfg.T_obs <= 0
        errs{end+1} = 'T_obs must be > 0.';
    end
end
if isfield(cfg, 'dt') && isnumeric(cfg.dt) && ~isempty(cfg.dt)
    if cfg.dt <= 0
        errs{end+1} = 'dt (output time step) must be > 0.';
    elseif isfield(cfg, 'T_obs') && isnumeric(cfg.T_obs) && ~isempty(cfg.T_obs) && cfg.dt > cfg.T_obs
        errs{end+1} = sprintf('dt = %g > T_obs = %g: output step exceeds total time.', ...
            cfg.dt, cfg.T_obs);
    end
end

% Material
if isfield(cfg, 'E') && isnumeric(cfg.E) && ~isempty(cfg.E) && cfg.E <= 0
    errs{end+1} = 'Young''s modulus E must be > 0.';
end
if isfield(cfg, 'r_pili') && isnumeric(cfg.r_pili) && ~isempty(cfg.r_pili) && cfg.r_pili <= 0
    errs{end+1} = 'Pili radius r_pili must be > 0.';
end
if isfield(cfg, 'sigma_p') && isnumeric(cfg.sigma_p) && ~isempty(cfg.sigma_p)
    if cfg.sigma_p < 0 || cfg.sigma_p >= 0.5
        errs{end+1} = sprintf('Poisson ratio sigma_p = %g outside valid range [0, 0.5).', ...
            cfg.sigma_p);
    end
end

% Fluid
if isfield(cfg, 'mu') && isnumeric(cfg.mu) && ~isempty(cfg.mu) && cfg.mu <= 0
    errs{end+1} = 'Dynamic viscosity mu must be > 0.';
end

% Numerics
if isfield(cfg, 'tlim') && isnumeric(cfg.tlim) && ~isempty(cfg.tlim) && cfg.tlim <= 0
    errs{end+1} = 'ODE timeout tlim must be > 0.';
end
if isfield(cfg, 'animation_fps') && isnumeric(cfg.animation_fps) ...
        && ~isempty(cfg.animation_fps) && cfg.animation_fps <= 0
    errs{end+1} = 'animation_fps must be > 0.';
end

% P2: delta (Euler-angle singularity threshold) — was unchecked; extreme
% values (delta <= 0 or >= pi/4) break singularity_event in run_simulation.
if isfield(cfg, 'delta') && isnumeric(cfg.delta) && ~isempty(cfg.delta) && isfinite(cfg.delta)
    if cfg.delta <= 0
        errs{end+1} = sprintf('delta = %g must be > 0 (Euler-angle singularity margin).', cfg.delta);
    elseif cfg.delta >= pi/4
        warns{end+1} = sprintf(['delta = %g >= pi/4 (%.4f): singularity event may trigger ', ...
            'prematurely, shrinking the ODE integration step. Recommended: pi/20 ~ pi/10.'], ...
            cfg.delta, pi/4);
    end
end

% P6 / A1: alpha_0 range — moved here from init_coupled_params so all
% parameter errors are reported in one pass.
%   - parametric mode: alpha_0 is always consumed to build q_0_init.
%   - custom mode + empty initial_orientation: alpha_0 is consumed as a
%     fallback (init_coupled_params L274-276), so the range must hold.
%   - custom mode + non-empty initial_orientation: alpha_0 is ignored
%     (a warning is emitted by init_coupled_params), so no range check
%     is needed here.
alpha0_used = false;
if strcmp(geometry_mode, 'parametric')
    alpha0_used = true;
elseif strcmp(geometry_mode, 'custom') ...
        && (~isfield(cfg, 'initial_orientation') || isempty(cfg.initial_orientation))
    alpha0_used = true;
end
if alpha0_used ...
        && isfield(cfg, 'alpha_0') && isnumeric(cfg.alpha_0) ...
        && ~isempty(cfg.alpha_0) && isfinite(cfg.alpha_0)
    if cfg.alpha_0 < -pi/2 || cfg.alpha_0 > pi/2
        errs{end+1} = sprintf('alpha_0 = %.4f rad is outside valid range [-pi/2, pi/2].', ...
            cfg.alpha_0);
    end
end

% Brownian
if isfield(cfg, 'enable_brown') && logical(cfg.enable_brown)
    if isfield(cfg, 'dt_B') && isnumeric(cfg.dt_B) && ~isempty(cfg.dt_B) && cfg.dt_B <= 0
        errs{end+1} = 'Brownian time step dt_B must be > 0 when enable_brown = true.';
    end
    if ~isfield(cfg, 'Temper') || ~isnumeric(cfg.Temper) || ~isfinite(cfg.Temper) || cfg.Temper <= 0
        errs{end+1} = 'Temper must be a finite positive temperature [K] when enable_brown = true.';
    end
end

%% ---- Warnings (both modes) ----

if strcmp(geometry_mode, 'parametric') && isfield(cfg, 'N_body') ...
        && isnumeric(cfg.N_body) && ~isnan(cfg.N_body) ...
        && cfg.N_body >= 10 && cfg.N_body < 100
    warns{end+1} = sprintf('N_body = %d < 100: hydrodynamic accuracy may be poor.', cfg.N_body);
end

if isfield(cfg, 'seg_per_um') && isnumeric(cfg.seg_per_um) ...
        && ~isempty(cfg.seg_per_um) && cfg.seg_per_um < 5
    warns{end+1} = sprintf('seg_per_um = %d < 5: RSS accuracy may be severely degraded.', ...
        cfg.seg_per_um);
end

if isfield(cfg, 'rho_bac') && isnumeric(cfg.rho_bac) && ~isempty(cfg.rho_bac) ...
        && isfield(cfg, 'rho_fluid') && isnumeric(cfg.rho_fluid) && ~isempty(cfg.rho_fluid) ...
        && cfg.rho_bac <= cfg.rho_fluid
    warns{end+1} = sprintf('rho_bac (%.1f) <= rho_fluid (%.3f): bacterium will not sink.', ...
        cfg.rho_bac, cfg.rho_fluid);
end

% P9: improved wording — make it clear that 1e-6 ~ 1e-3 is the recommended
% range, so the default 1e-4 is not mistaken for "out of range".
if isfield(cfg, 'AbsTol') && isnumeric(cfg.AbsTol) && ~isempty(cfg.AbsTol)
    if cfg.AbsTol <= 0 || cfg.AbsTol > 1
        warns{end+1} = sprintf(['AbsTol = %g outside valid range (0, 1]. ', ...
            'Recommended: 1e-6 ~ 1e-3 for double-precision ODE accuracy.'], cfg.AbsTol);
    end
end
if isfield(cfg, 'RelTol') && isnumeric(cfg.RelTol) && ~isempty(cfg.RelTol)
    if cfg.RelTol <= 0 || cfg.RelTol > 1
        warns{end+1} = sprintf(['RelTol = %g outside valid range (0, 1]. ', ...
            'Recommended: 1e-6 ~ 1e-3 for double-precision ODE accuracy.'], cfg.RelTol);
    end
end

if isfield(cfg, 'enable_brown') && logical(cfg.enable_brown) ...
        && isfield(cfg, 'dt_B') && isnumeric(cfg.dt_B) && ~isempty(cfg.dt_B) ...
        && isfield(cfg, 'dt') && isnumeric(cfg.dt) && ~isempty(cfg.dt) ...
        && cfg.dt_B > 0 && cfg.dt_B >= cfg.dt
    warns{end+1} = sprintf('dt_B = %g >= dt = %g: Brownian step >= output step.', ...
        cfg.dt_B, cfg.dt);
end

% dt_B vs. Stokes-Einstein (only in parametric mode where major_axis is known)
if isfield(cfg, 'enable_brown') && logical(cfg.enable_brown) ...
        && isfield(cfg, 'dt_B') && isnumeric(cfg.dt_B) && ~isempty(cfg.dt_B) && cfg.dt_B > 0 ...
        && isfield(cfg, 'Temper') && isnumeric(cfg.Temper) && ~isempty(cfg.Temper) ...
            && isfinite(cfg.Temper) && cfg.Temper > 0 ...
        && isfield(cfg, 'mu') && isnumeric(cfg.mu) && ~isempty(cfg.mu) && cfg.mu > 0 ...
        && isfield(cfg, 'major_axis') && isnumeric(cfg.major_axis) && ~isnan(cfg.major_axis) && cfg.major_axis > 0 ...
        && isfield(cfg, 'minor_axis') && isnumeric(cfg.minor_axis) && ~isnan(cfg.minor_axis) && cfg.minor_axis > 0
    kB    = 1.380649e-23;
    a_eff = max(cfg.major_axis, cfg.minor_axis);
    D_T   = kB * cfg.Temper / (6 * pi * cfg.mu * a_eff);
    D_R   = kB * cfg.Temper / (8 * pi * cfg.mu * a_eff^3);
    tau_trans = a_eff^2 / D_T;
    tau_rot   = 1 / (2 * D_R);   % first-order orientational correlation time
    tau_min   = min(tau_trans, tau_rot);
    if isfinite(tau_min) && tau_min > 0 && cfg.dt_B > 0.1 * tau_min
        warns{end+1} = sprintf(['dt_B = %.3g s exceeds 10%% of estimated ', ...
            'min(tau_trans, tau_rot) ~ %.3g s (Stokes-Einstein with a = %.2g m). ', ...
            'Consider dt_B <= %.3g s for accurate Brownian statistics.'], ...
            cfg.dt_B, tau_min, a_eff, 0.1 * tau_min);
    end
end

% Background flow
if isfield(cfg, 'background_flow') && isstruct(cfg.background_flow)
    bg = cfg.background_flow;
    valid_bg = {'off', 'uniform', 'shear_z'};
    bg_type = 'off';
    if isfield(bg, 'type') && ~isempty(bg.type)
        bg_type = char(bg.type);
    end
    if ~any(strcmp(bg_type, valid_bg))
        errs{end+1} = sprintf('background_flow.type must be one of {off, uniform, shear_z}, got ''%s''.', bg_type);
    end
    if ~strcmp(bg_type, 'off')
        if ~isfield(bg, 'U_max') || ~isnumeric(bg.U_max) || bg.U_max <= 0
            errs{end+1} = 'background_flow.U_max must be > 0 when flow is enabled.';
        end
        if isfield(bg, 'direction') && numel(bg.direction) ~= 3
            errs{end+1} = 'background_flow.direction must be a 3-element vector.';
        end
        if isfield(bg, 'shear_axis') && numel(bg.shear_axis) ~= 3
            errs{end+1} = 'background_flow.shear_axis must be a 3-element vector.';
        end
        if strcmp(bg_type, 'shear_z')
            if ~isfield(bg, 'shear_height_m') || ~isnumeric(bg.shear_height_m) || bg.shear_height_m <= 0
                errs{end+1} = 'background_flow.shear_height_m must be > 0 for shear_z flow.';
            end
            % P5: warn when direction is nearly parallel to shear_axis.
            % In coupled_background_flow, omega_lab = (U_max/H) * cross(shear_axis, direction).
            % Parallel inputs yield zero vorticity, so pili experience no
            % shear-driven rotation even though the translational profile
            % is non-zero.  This is a common user intent error.
            if isfield(bg, 'direction') && isfield(bg, 'shear_axis') ...
                    && isnumeric(bg.direction) && isnumeric(bg.shear_axis) ...
                    && numel(bg.direction) == 3 && numel(bg.shear_axis) == 3
                n_d = norm(bg.direction(:));
                n_s = norm(bg.shear_axis(:));
                if n_d > 1e-30 && n_s > 1e-30
                    cos_ang = dot(bg.direction(:)/n_d, bg.shear_axis(:)/n_s);
                    if abs(cos_ang) > 0.99
                        warns{end+1} = sprintf(['background_flow: direction nearly parallel ', ...
                            'to shear_axis (cos(theta) = %.3f). Vorticity = shear_axis x direction ', ...
                            'will be ~0, so pili experience no shear-driven rotation despite a ', ...
                            'non-zero translational flow. Verify this is intended.'], cos_ang);
                    end
                end
            end
        end
    end
    if ~strcmp(bg_type, 'off') && isfield(cfg, 'enable_brown') && logical(cfg.enable_brown)
        warns{end+1} = ['Background flow + Brownian: flow is deterministic; ' ...
            'diffusion tensor computed without flow contribution.'];
    end
end

%% ---- Report ----

for i = 1:numel(warns)
    logf_h('[VALIDATE] WARNING: %s\n', warns{i});
    warning('validate_params:SuboptimalValue', '%s', warns{i});
end

if ~isempty(errs)
    % P10: build the aggregated message with sprintf in one pass instead of
    % cellfun + strjoin, which avoids the intermediate cell allocation when
    % the error list is long (e.g. per-pilus errors across many pili).
    msg_parts = cell(numel(errs), 1);
    for i = 1:numel(errs)
        msg_parts{i} = sprintf('  - %s', errs{i});
    end
    msg = strjoin(msg_parts, newline);
    error('validate_params:InvalidParams', ...
        'Parameter validation failed with %d error(s):\n%s', numel(errs), msg);
end

end


%% ==================== Local helper functions ====================

function [errs, warns] = validate_pili_ext(pili_ext, M_pili, errs, warns, aux)
%VALIDATE_PILI_EXT  Validate per-pilus extended fields (intrinsic_curvature,
%                   active_moment, root_coupling, motor_torque).
%
%  Used by both the parametric branch (M_pili = size(Pili_Matrix, 2))
%  and the custom branch (M_pili = numel(custom_geometry.pili)).
%
%  aux (optional) is a struct with auxiliary info used for the
%  stall-curve diagnostic on free_torsion + constant motors.  Fields:
%    .pili_lengths_m  - [1 x M_pili] pili total lengths [m]
%    .r_pili          - pili radius [m]
%    .mu              - dynamic viscosity [Pa*s]
%  An empty aux disables the stall-curve check (used by the custom-mode
%  call where centreline lengths require reading the user files).

    if nargin < 5
        aux = struct();
    end

    valid_ic = {'frozen', 'zero', 'constant', 'helix'};
    valid_am = {'off', 'traveling_wave', 'custom'};
    valid_rc = {'clamped', 'free_torsion'};
    valid_mt = {'off', 'constant'};

    aux_has_lengths = isstruct(aux) && isfield(aux, 'pili_lengths_m') && ...
                      isfield(aux, 'r_pili') && isfield(aux, 'mu') && ...
                      ~isempty(aux.pili_lengths_m);

    for kk = 1:min(numel(pili_ext), M_pili)
        pe = pili_ext(kk);

        % ---- Intrinsic curvature ----
        % P7: unified string/struct acceptance (same as motor_torque).
        if isfield(pe, 'intrinsic_curvature') && ~isempty(pe.intrinsic_curvature)
            ic_val = pe.intrinsic_curvature;
            if ischar(ic_val) || isstring(ic_val)
                ic_val = struct('type', char(ic_val));
            end
            if isstruct(ic_val) && isfield(ic_val, 'type')
                ic_type = char(ic_val.type);
            else
                ic_type = '';
            end
            if ~any(strcmp(ic_type, valid_ic))
                errs{end+1} = sprintf(['Pili %d: intrinsic_curvature type must be ', ...
                    'one of {frozen, zero, constant, helix}, got ''%s''.'], kk, ic_type);
            end
            if strcmp(ic_type, 'constant')
                if isstruct(ic_val) && isfield(ic_val, 'kappa')
                    if numel(ic_val.kappa) ~= 3
                        errs{end+1} = sprintf(['Pili %d: intrinsic_curvature.kappa must be ', ...
                            '3-element, got %d.'], kk, numel(ic_val.kappa));
                    end
                else
                    errs{end+1} = sprintf(['Pili %d: intrinsic_curvature ''constant'' ', ...
                        'requires .kappa field.'], kk);
                end
            end
            if strcmp(ic_type, 'helix')
                if isstruct(ic_val) && isfield(ic_val, 'kappa_d2') && isfield(ic_val, 'kappa_d3')
                    if ~isnumeric(ic_val.kappa_d2) || ~isscalar(ic_val.kappa_d2) || ~isfinite(ic_val.kappa_d2)
                        errs{end+1} = sprintf(['Pili %d: intrinsic_curvature.kappa_d2 must be ', ...
                            'a finite scalar.'], kk);
                    end
                    if ~isnumeric(ic_val.kappa_d3) || ~isscalar(ic_val.kappa_d3) || ~isfinite(ic_val.kappa_d3)
                        errs{end+1} = sprintf(['Pili %d: intrinsic_curvature.kappa_d3 must be ', ...
                            'a finite scalar.'], kk);
                    end
                else
                    errs{end+1} = sprintf(['Pili %d: intrinsic_curvature ''helix'' ', ...
                        'requires .kappa_d2 and .kappa_d3 fields.'], kk);
                end
            end
        end

        % ---- Active moment ----
        % P7: unified string/struct acceptance (matches intrinsic_curvature
        % and motor_torque).  parse_pili_config / parse_custom_pili_ext
        % always emit structs, but direct callers may pass strings.
        if isfield(pe, 'active_moment') && ~isempty(pe.active_moment)
            am_val = pe.active_moment;
            if ischar(am_val) || isstring(am_val)
                am_val = struct('type', char(am_val));
            end
            if isstruct(am_val)
                am_type = '';
                if isfield(am_val, 'type')
                    am_type = char(am_val.type);
                end
                if ~any(strcmp(am_type, valid_am))
                    errs{end+1} = sprintf(['Pili %d: active_moment type must be one of ', ...
                        '{off, traveling_wave, custom}, got ''%s''.'], kk, am_type);
                end
                if strcmp(am_type, 'traveling_wave')
                    if ~isfield(am_val, 'amplitude_Nm') || ~isnumeric(am_val.amplitude_Nm) ...
                            || am_val.amplitude_Nm < 0
                        errs{end+1} = sprintf('Pili %d: active_moment.amplitude_Nm must be >= 0.', kk);
                    end
                    if ~isfield(am_val, 'omega_Hz') || ~isnumeric(am_val.omega_Hz) || am_val.omega_Hz <= 0
                        errs{end+1} = sprintf('Pili %d: active_moment.omega_Hz must be > 0.', kk);
                    end
                    if isfield(am_val, 'wavenumber') && isnumeric(am_val.wavenumber) && am_val.wavenumber < 0
                        warns{end+1} = sprintf(['Pili %d: active_moment.wavenumber = %.4g is ', ...
                            'negative. Wave propagates in -s direction (from tip to root).'], ...
                            kk, am_val.wavenumber);
                    end
                    if isfield(am_val, 'd3_amp') && isnumeric(am_val.d3_amp) && am_val.d3_amp ~= 0
                        warns{end+1} = sprintf(['Pili %d: active_moment.d3_amp = %.4g (non-zero). ', ...
                            'This drives torsion rather than bending.'], kk, am_val.d3_amp);
                    end
                    if isfield(am_val, 'phase_rad')
                        if ~isnumeric(am_val.phase_rad) || ~isscalar(am_val.phase_rad) || ~isfinite(am_val.phase_rad)
                            errs{end+1} = sprintf('Pili %d: active_moment.phase_rad must be a finite scalar.', kk);
                        end
                    end
                end
            else
                errs{end+1} = sprintf('Pili %d: active_moment must be a struct or string, got %s.', ...
                    kk, class(am_val));
            end
        end

        % ---- Root coupling ----
        if isfield(pe, 'root_coupling') && ~isempty(pe.root_coupling)
            rc_val = char(pe.root_coupling);
            if ~any(strcmp(rc_val, valid_rc))
                errs{end+1} = sprintf(['Pili %d: root_coupling must be one of ', ...
                    '{clamped, free_torsion}, got ''%s''.'], kk, rc_val);
            end
        end

        % ---- Motor torque ----
        if isfield(pe, 'motor_torque') && ~isempty(pe.motor_torque)
            mt_val = pe.motor_torque;
            if ischar(mt_val) || isstring(mt_val)
                mt_val = struct('type', char(mt_val));
            end
            if ~isstruct(mt_val)
                errs{end+1} = sprintf('Pili %d: motor_torque must be a struct or string, got %s.', ...
                    kk, class(mt_val));
            else
                mt_type = '';
                if isfield(mt_val, 'type')
                    mt_type = char(mt_val.type);
                end
                if ~any(strcmp(mt_type, valid_mt))
                    errs{end+1} = sprintf(['Pili %d: motor_torque type must be one of ', ...
                        '{off, constant}, got ''%s''.'], kk, mt_type);
                end
                if strcmp(mt_type, 'constant')
                    if ~isfield(mt_val, 'torque_pN_nm') || ~isnumeric(mt_val.torque_pN_nm) ...
                            || ~isfinite(mt_val.torque_pN_nm) || mt_val.torque_pN_nm == 0
                        errs{end+1} = sprintf('Pili %d: motor_torque.torque_pN_nm must be a finite non-zero value.', kk);
                    end
                    if ~isfield(mt_val, 'max_omega_Hz') || ~isnumeric(mt_val.max_omega_Hz) ...
                            || ~isfinite(mt_val.max_omega_Hz) || mt_val.max_omega_Hz <= 0
                        errs{end+1} = sprintf('Pili %d: motor_torque.max_omega_Hz must be > 0.', kk);
                    end
                end

                rc_actual = 'clamped';
                if isfield(pe, 'root_coupling') && ~isempty(pe.root_coupling)
                    rc_actual = char(pe.root_coupling);
                end
                if ~strcmp(mt_type, 'off')
                    if strcmp(rc_actual, 'clamped')
                        errs{end+1} = sprintf(['Pili %d: motor_torque is active (''%s'') but ', ...
                            'root_coupling = clamped. Motor torque has no effect on clamped ', ...
                            'pili. Set root_coupling = ''free_torsion'' to enable torsional dynamics.'], ...
                            kk, mt_type);
                    end

                    % --- Stall-curve diagnostic (Issue #7+#9, v7.3.0) ---
                    %
                    % The motor model is tau_eff = tau_motor*(1 - |omega|/omega_max).
                    % In steady state, motor torque balances RFT torsional drag:
                    %     tau_motor*(1 - omega/omega_max) = xi_t * omega
                    % with xi_t ~ 4*pi*mu*r_pili^2*L_pilus (no-slip rotational drag
                    % of a slender rod about its own axis, RFT leading order).
                    % Solving: omega_steady = tau_motor / (xi_t + tau_motor/omega_max)
                    % => omega_steady/omega_max = tau_motor / (xi_t*omega_max + tau_motor)
                    %
                    % If this ratio is too close to 1 the motor is operating at
                    % stall, where:
                    %   (a) the speed-limiter dominates the dynamics, so the
                    %       ode15s step may shrink as it resolves the stiff
                    %       torque-balance (the stall-curve torque itself is
                    %       now solved EXACTLY in-call -- see ode_coupled_rhs.m
                    %       Step 9 -- so there is no lag-approximation error,
                    %       but the underlying dynamics are still stiff);
                    %   (b) results become highly sensitive to the empirical
                    %       linear stall curve and rarely represent real
                    %       biological motors well.
                    if strcmp(mt_type, 'constant') && strcmp(rc_actual, 'free_torsion') && ...
                       aux_has_lengths && kk <= numel(aux.pili_lengths_m) && ...
                       isfield(mt_val, 'torque_pN_nm') && isnumeric(mt_val.torque_pN_nm) ...
                           && isfinite(mt_val.torque_pN_nm) && ...
                       isfield(mt_val, 'max_omega_Hz') && isnumeric(mt_val.max_omega_Hz) ...
                           && mt_val.max_omega_Hz > 0
                        L_k         = aux.pili_lengths_m(kk);            % [m]
                        tau_motor_N = abs(mt_val.torque_pN_nm) * 1e-21;  % pN*nm -> N*m
                        omega_max   = 2*pi * mt_val.max_omega_Hz;        % Hz -> rad/s
                        xi_t        = 4*pi*aux.mu*aux.r_pili^2*L_k;      % N*m*s/rad
                        ratio       = tau_motor_N / (xi_t*omega_max + tau_motor_N);
                        if ratio >= 0.7
                            warns{end+1} = sprintf(['Pili %d: motor is operating ', ...
                                'near stall (omega_steady/omega_max ~ %.2f). The ', ...
                                'speed-limiter dominates the dynamics; although ', ...
                                'tau_eff(dpsi_1) is now solved exactly in-call, the ', ...
                                'stiff torque-balance may still shrink the ode15s step. ', ...
                                'Consider a larger max_omega_Hz, smaller torque_pN_nm, ', ...
                                'or a longer pilus.'], kk, ratio);
                        elseif ratio >= 0.4
                            warns{end+1} = sprintf(['Pili %d: motor steady-state ', ...
                                'speed reaches ~%.0f%% of max_omega_Hz under RFT drag ', ...
                                'alone. Body+pilus hydrodynamic feedback may push this ', ...
                                'higher; if step size shrinks unexpectedly, increase ', ...
                                'max_omega_Hz or check stall behaviour.'], ...
                                kk, 100*ratio);
                        end
                    end
                elseif strcmp(rc_actual, 'free_torsion')
                    warns{end+1} = sprintf(['Pili %d: root_coupling = free_torsion but ', ...
                        'motor_torque = off. Psi_1 will evolve purely under elastic torsion.'], kk);
                end
            end
        end
    end
end
