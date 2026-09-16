function [T_ice, f_ice, f_liq, k_eff, ok, iter] = solve_column_temperature( ...
      T_sfc, T_ice, f_ice, f_liq, dz, delz, fn, dt, settings)
   %SOLVE_COLUMN_TEMPERATURE Solve the 1-dimensional column conduction equation.
   %
   %  [T_ice, f_ice, f_liq, k_eff, ok, iter] = ...
   %     icemodel.column.solve_column_temperature(T_sfc, T_ice, f_ice, ...
   %     f_liq, dz, delz, fn, dt, settings)
   %
   % SETTINGS is the struct from icemodel.couplers.initialize_solver_settings.
   % This solver reads these fields:
   %  - tol: convergence tolerance on the node temperature change [K]
   %  - maxiter: iteration limit
   %  - alpha: relaxation factor, set to 1 when maxiter is 1
   %  - debug: true to dump a failed solve
   % The fields use_aitken and jumpmax stay in SETTINGS so the disabled Aitken
   % block below can use them.
   %
   % See also: icemodel.couplers.solve_skin_surface_column
   %
   %#codegen

   persistent cv_ice cv_liq Ls
   if isempty(cv_ice)
      [cv_ice, cv_liq, Ls] = icemodel.physicalConstant( ...
         'cv_ice', 'cv_liq', 'Ls');
   end

   % Top and bottom node indices
   N = 1;
   S = numel(T_ice);

   % Solver options
   tol = settings.tol;
   maxiter = settings.maxiter;
   alpha = settings.alpha;
   debug = settings.debug;
   if maxiter == 1
      alpha = 1;
   end

   % Thermal conductivity without vapor diffusion (skinmodel).
   k_eff = icemodel.column.bulk_thermal_conductivity(T_ice, f_ice, f_liq, 0);

   % The enthalpy budget excludes the vapor density derivative (skinmodel).
   drovdT = 0;

   % To reinstate vapor-aware conductivity and enthalpy:
   % [~, drovdT] = icemodel.vapor.saturation_vapor_density(T_ice, f_liq);
   % k_vap = icemodel.vapor.vapor_thermal_conductivity(T_ice, f_liq, drovdT);
   % k_eff = icemodel.column.bulk_thermal_conductivity(T_ice, f_ice, f_liq, k_vap);
   %
   % The iterations need the same update. See solve_column_enthalpy.

   % Initial past Picard iterates for Aitken-acceleration
   % T_1 = nan(size(T_ice));
   % T_2 = nan(size(T_ice));

   % Iterate to solve the nonlinear heat equation (p. 47)
   ok = false;
   for iter = 1:maxiter

      % Capture current T_ice iterate
      T_iter = T_ice;

      % Compute gamma at the control volume interfaces (eq. 4.9, p. 45) (JJ+1)
      g_ns = [k_eff(N); k_eff(N:S); k_eff(S)];
      gb_ns = 1.0 ./ ( (1.0 - fn) ./ g_ns(N:S+1) + fn ./ g_ns(N+1:S+2));

      % Compute the enthalpy coefficient for each c.v. for the current timestep
      aP0 = (cv_ice * f_ice + cv_liq * f_liq ...
         + Ls * (1 - f_liq - f_ice) .* drovdT) .* dz / dt;

      % Compute the aN and aS coefficients
      aN = gb_ns(N:S)     ./ delz(N:S);
      aS = gb_ns(N+1:S+1) ./ delz(N+1:S+1);

      % Account for the boundary conditions.
      bc_N = aN(N) * T_sfc;
      bc_S = 0.0;
      aS(S) = 0.0;

      % Compute the aP coefficient and solution vector b
      aP = aN(N:S) + aS(N:S) + aP0(N:S);
      b = aP0(N:S) .* T_ice(N:S);

      % Account for Dirichlet upper and Neumann lower boundary conditions
      b(N) = b(N) + bc_N;
      b(S) = b(S) + bc_S;

      % Solve the equation
      T_ice = icemodel.numerics.trisolve(-aN, aP, -aS, b);

      % Prep for next iteration
      if all(abs(T_ice - T_iter) < tol)
         ok = true;
         break
      end

      % Apply relaxation
      T_ice = alpha * T_ice + (1 - alpha) * T_iter;

      % Aitken acceleration (node-by-node) with relaxed value as fallback.
      % if settings.use_aitken
      %    T_0 = T_ice;
      %    for mm = 1:numel(T_ice)
      %       T_ice(mm) = icemodel.numerics.aitkenscalar( ...
      %          T_2(mm), T_1(mm), T_0(mm), T_ice(mm), settings.jumpmax);
      %    end
      %    T_2 = T_1;
      %    T_1 = T_0;
      % end

      % Update thermal conductivity. This also ensures T_ice-k_eff consistency
      % on the final iteration.
      k_eff = icemodel.column.bulk_thermal_conductivity(T_ice, f_ice, f_liq, 0);

      % To reinstate k_vap:
      % k_eff = icemodel.column.bulk_thermal_conductivity( ...
      %    T_ice, f_ice, f_liq, k_vap);
   end

   % Debug dump on a failed solve
   if ~ok && debug
      dumpSkinSolveFailure( ...
         T_ice, f_ice, f_liq, k_eff, dz, delz, dt, T_sfc, iter, maxiter);
   end
end

function dumpSkinSolveFailure(T_ice, f_ice, f_liq, k_eff, dz, delz, dt, ...
      T_sfc, iter, maxiter)
   %DUMPSKINSOLVEFAILURE Save skin conduction solver diagnostics on demand.

   debug_file = getenv('ICEMODEL_DEBUG_SKINSOLVE_FILE');
   if isempty(debug_file)
      return
   end

   debug_state = struct();
   debug_state.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   debug_state.T_ice = T_ice;
   debug_state.f_ice = f_ice;
   debug_state.f_liq = f_liq;
   debug_state.k_eff = k_eff;
   debug_state.dz = dz;
   debug_state.delz = delz;
   debug_state.dt = dt;
   debug_state.T_sfc = T_sfc;
   debug_state.iter = iter;
   debug_state.maxiter = maxiter;

   save(debug_file, 'debug_state');
end
