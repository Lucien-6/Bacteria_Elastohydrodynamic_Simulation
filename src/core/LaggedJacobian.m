classdef LaggedJacobian < handle
%LAGGEDJACOBIAN  Lagged (re-used) finite-difference Jacobian for ode15s
%
%  Description:
%    ode15s evaluates a dense finite-difference Jacobian (n_ode extra RHS
%    calls) every time it is started.  A Brownian run restarts ode15s at
%    every dt_B node, so this fixed cost dominates (70-85 % of the run time
%    of a 4-pilus case).  This handle class keeps the last exact
%    finite-difference Jacobian and hands it back to ode15s at the start of
%    the following segments.
%
%  Mathematical status:
%    ode15s uses the Jacobian only inside the modified-Newton iteration
%    matrix W = I - gamma*h*J of the implicit BDF/NDF step.  A stale J slows
%    the Newton iteration (contraction factor rho(I - W_old^{-1} W_true));
%    it does not change the discrete equations that the iteration solves
%    nor the local error control.  The solution therefore differs from the
%    fresh-Jacobian solution only at the level of the Newton / integration
%    tolerance.  If Newton fails with a stale J, ode15s asks for J again
%    within the same call; that second request is answered with a freshly
%    computed exact finite-difference Jacobian (refresh on demand).
%
%  Refresh policy:
%    - first request of a segment      : cached J (if any, and not too old)
%    - second or later request          : fresh exact J   (Newton failed)
%    - cache age >= max_age nodes       : fresh exact J   (the cache is
%                                         therefore used for at most
%                                         max_age - 1 following nodes)
%    - after invalidate()               : fresh exact J   (Euler-basis
%                                         rotation re-parameterises the
%                                         pilus angles)
%
%  Usage (see run_simulation.m):
%    jc = LaggedJacobian(50);
%    jc.beginSegment(@(t, z) rhs(t, z));
%    opts = odeset(opts, 'Jacobian', @(t, z) jc.evaluate(t, z));
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-09

    properties
        max_age = 50;            % maximum age of the cache [Brownian nodes]
    end

    properties (SetAccess = private)
        cache            = [];   % cached exact Jacobian
        rhs_fun          = [];   % RHS handle used for finite differences
        calls_in_segment = 0;    % Jacobian requests in the current ode15s call
        age              = 0;    % Brownian nodes since the cache was computed
        n_fresh          = 0;    % exact Jacobians computed
        n_lagged         = 0;    % requests answered from the cache
        n_on_demand      = 0;    % refreshes triggered by a second request
        n_age_refresh    = 0;    % refreshes triggered by cache age
        n_invalidated    = 0;    % explicit invalidations (basis rotation)
        n_fallback       = 0;    % segments re-run with the default Jacobian
        n_fd_rhs         = 0;    % RHS calls spent on finite differences
    end

    methods
        function obj = LaggedJacobian(max_age)
            if nargin >= 1 && ~isempty(max_age)
                if ~(isnumeric(max_age) && isscalar(max_age) && ...
                        isfinite(max_age) && max_age >= 1)
                    error('LaggedJacobian:InvalidMaxAge', ...
                        'max_age must be a finite scalar >= 1 (Brownian nodes).');
                end
                obj.max_age = double(max_age);
            end
        end

        function beginSegment(obj, rhs_fun)
        %BEGINSEGMENT  Call before each ode15s call.
            obj.rhs_fun = rhs_fun;
            obj.calls_in_segment = 0;
        end

        function tick(obj)
        %TICK  Call once at the start of every new Brownian node.
            obj.age = obj.age + 1;
        end

        function invalidate(obj)
        %INVALIDATE  Drop the cache (e.g. after an Euler-basis rotation).
            if ~isempty(obj.cache)
                obj.n_invalidated = obj.n_invalidated + 1;
            end
            obj.cache = [];
        end

        function noteFallback(obj)
        %NOTEFALLBACK  Record that a segment was re-run with default Jacobian.
            obj.n_fallback = obj.n_fallback + 1;
            obj.cache = [];
        end

        function J = evaluate(obj, t, y)
        %EVALUATE  Jacobian callback for ode15s ('Jacobian' option).
            obj.calls_in_segment = obj.calls_in_segment + 1;
            stale_request = obj.calls_in_segment > 1;
            too_old       = obj.age >= obj.max_age;
            need_fresh = isempty(obj.cache) || stale_request || too_old ...
                || size(obj.cache, 2) ~= numel(y);
            if need_fresh
                if stale_request
                    obj.n_on_demand = obj.n_on_demand + 1;
                elseif too_old && ~isempty(obj.cache)
                    obj.n_age_refresh = obj.n_age_refresh + 1;
                end
                f0 = obj.rhs_fun(t, y);
                Jn = LaggedJacobian.fd_jacobian(obj.rhs_fun, t, y, f0);
                if ~all(isfinite(Jn(:)))
                    error('LaggedJacobian:NonFinite', ...
                        'Finite-difference Jacobian contains non-finite entries.');
                end
                obj.cache = Jn;
                obj.age = 0;
                obj.n_fresh = obj.n_fresh + 1;
                obj.n_fd_rhs = obj.n_fd_rhs + 1 + numel(y);
            else
                obj.n_lagged = obj.n_lagged + 1;
            end
            J = obj.cache;
        end

        function s = stats(obj)
        %STATS  Counters as a plain struct (safe to save in results).
            s = struct('enabled', true, 'max_age', obj.max_age, ...
                'n_fresh', obj.n_fresh, 'n_lagged', obj.n_lagged, ...
                'n_on_demand', obj.n_on_demand, ...
                'n_age_refresh', obj.n_age_refresh, ...
                'n_invalidated', obj.n_invalidated, ...
                'n_fallback', obj.n_fallback, 'n_fd_rhs', obj.n_fd_rhs);
        end
    end

    methods (Static)
        function J = fd_jacobian(f, t, y, f0)
        %FD_JACOBIAN  Forward-difference Jacobian, step 1e-7*max(|y_j|,1).
            n = numel(y);
            J = zeros(numel(f0), n);
            for j = 1:n
                dy = 1e-7 * max(abs(y(j)), 1);
                yj = y;
                yj(j) = y(j) + dy;
                J(:, j) = (f(t, yj) - f0) / dy;
            end
        end
    end
end
