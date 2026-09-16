function [T_ice, f_ice, f_liq, k_eff, U_vap, L_vap, ok, iter, a1, err] = ...
      solve_column_enthalpy(T_sfc, T_ice, f_ice, f_liq, Fc, Fp, Sc, Sp, dz, ...
      delz, fn, dt, f_res_por, settings)
   %SOLVE_COLUMN_ENTHALPY Solve the column enthalpy balance.
   %
   %  [T_ice, f_ice, f_liq, k_eff, U_vap, L_vap, ok, iter, a1, err] = ...
   %     icemodel.column.solve_column_enthalpy(T_sfc, T_ice, f_ice, f_liq, ...
   %     Fc, Fp, Sc, Sp, dz, delz, fn, dt, f_res_por, settings)
   %
   % T_ICE, F_ICE, and F_LIQ are the accepted temperature and phase fractions.
   % K_EFF is the vapor-free node conductivity. U_VAP and L_VAP are the
   % accepted face vapor flux and donor latent heat. OK is false when the
   % melt-zone transform or iteration limit rejects the solve.
   %
   % U_VAP is positive downward [kg m-2 s-1] and is zero at both boundary
   % faces. L_VAP is Lv for a wet donor cell and Ls for a dry donor cell.
   % F_RES_POR sets the residual-liquid threshold used to select that phase.
   %
   % SETTINGS is the struct from icemodel.couplers.initialize_solver_settings.
   % This solver reads these fields:
   %  - solver: upper boundary type passed to assemble_enthalpy_system
   %  - tol: convergence tolerance on the node temperature change [K]
   %  - maxiter: iteration limit
   %  - debug: true to dump a failed solve
   % Aitken acceleration is disabled here. The fields alpha, use_aitken, and
   % jumpmax stay in SETTINGS so the disabled Aitken block below can use them.
   %
   % See also: icemodel.column.assemble_enthalpy_system,
   %  icemodel.column.vapor_transport_terms,
   %  icemodel.column.bulk_enthalpy
   %
   %#codegen

   persistent Ls
   if isempty(Ls)
      Ls = icemodel.physicalConstant('Ls');
   end

   % Read the solver controls this kernel uses.
   solver = settings.solver;
   tol = settings.tol;
   maxiter = settings.maxiter;
   debug = settings.debug;

   % Convert ice and liquid fractions to total water fraction.
   f_wat = icemodel.column.water_fraction(f_ice, f_liq);

   % Get the melt-zone liquid-fraction bounds at T_ice=TL and T_ice=TH.
   [f_liq_min, f_liq_max] = icemodel.column.meltzone_bounds(f_wat);

   % Compute vapor properties for the first iteration.
   [ro_vap, dro_vapdT] = icemodel.vapor.saturation_vapor_density(T_ice, f_liq);
   [~, De] = icemodel.vapor.vapor_thermal_conductivity(T_ice, f_liq, dro_vapdT);

   % Store the entering enthalpy [J m-3], used for the conservative solve.
   H_old = icemodel.column.bulk_enthalpy(T_ice, f_ice, f_liq, f_wat, ro_vap);

   % Store the entering temperature and liquid fraction, used for the
   % subsurface linearization error estimate.
   T_ice_old = T_ice;
   f_liq_old = f_liq;

   % Make the first convergence check fail.
   T_iter = T_ice_old + 2 * tol;

   % Initial past Picard iterates for Aitken-acceleration (disabled).
   % T_1 = nan(size(T_ice));
   % T_2 = nan(size(T_ice));

   % Initialize optional outputs.
   err = nan;
   U_vap = zeros(numel(T_ice) + 1, 1);
   L_vap = Ls * ones(numel(T_ice) + 1, 1);

   % Iterate to solve the nonlinear heat equation.
   ok = false;
   for iter = 0:maxiter-1

      % Build the vapor-free node conductivity and face vapor transport terms.
      k_eff = icemodel.column.bulk_thermal_conductivity(T_ice, f_ice, f_liq, 0);
      [k_eff_faces, ~, q_deferred_faces, U_vap, L_vap] = ...
         icemodel.column.vapor_transport_terms(T_ice, f_ice, f_liq, k_eff, ...
         ro_vap, dro_vapdT, De, delz, fn, f_res_por);

      % Update enthalpy and its temperature derivatives.
      [H, dHdT, dFdT] = icemodel.column.bulk_enthalpy( ...
         T_ice, f_ice, f_liq, f_wat, ro_vap, dro_vapdT);

      % Assemble the enthalpy equation.
      [aN, aP, aS, b, iM, a1, a2] = ...
         icemodel.column.assemble_enthalpy_system( ...
         T_ice, f_ice, f_liq, dHdT, dFdT, dro_vapdT, H - H_old, Sc, Sp, ...
         k_eff_faces, delz, dz, dt, T_sfc, Fc, Fp, solver, ...
         q_deferred_faces);

      % % Check diagonal dominance and condition number.
      % icemodel.checkdiags(aP, aN, aS)

      % Stop here so the state is updated on the final iteration.
      if all(abs(T_ice - T_iter) < tol)
         ok = true;
         break
      end

      % Save the current temperature before the next solve.
      T_iter = T_ice;

      % Solve the linearized enthalpy equation.
      T_ice = icemodel.numerics.trisolve(-aN, aP, -aS, b);

      % Convert enthalpy back to temperature and phase fractions.
      [T_ice, f_ice, f_liq, ok_mz] = icemodel.column.meltzone_transform( ...
         T_ice, T_iter, f_liq, f_wat, dFdT, f_liq_min, f_liq_max, iM, debug);

      % Return a rejected solve when the melt-zone transform fails.
      if ~ok_mz
         if debug
            dumpIceEnbalFailure("mztransform_rejected_state", ...
               T_ice, T_ice_old, T_iter, f_ice, f_liq, f_wat, k_eff, ...
               Sc, dt, T_sfc, iM, iter, maxiter, aN, aP, aS, b);
         end
         return
      end

      % Check the maximum water fraction.
      assertF(@() icemodel.column.assert_max_water(f_ice, f_liq));

      % Update vapor properties for the next iteration.
      [ro_vap, dro_vapdT] = icemodel.vapor.saturation_vapor_density( ...
         T_ice, f_liq);
      [~, De] = icemodel.vapor.vapor_thermal_conductivity( ...
         T_ice, f_liq, dro_vapdT);

      % Relaxation and Aitken (not implemented). Proper implementation requires
      % a liquid_fraction_function/meltzone_transform-consistent update of
      % T_ice, f_ice, and f_liq.
      % T_ice = min(Tf, settings.alpha * T_ice + (1 - settings.alpha) * T_iter);
      % if settings.use_aitken
      %    T_0 = T_ice;
      %    for mm = 1:numel(T_ice)
      %       T_ice(mm) = icemodel.numerics.aitkenscalar( ...
      %          T_2(mm), T_1(mm), T_0(mm), T_ice(mm), settings.jumpmax);
      %    end
      %    T_2 = T_1;
      %    T_1 = T_0;
      % end
   end

   % Return the iteration limit when the loop does not converge.
   if ~ok
      iter = maxiter;
   end

   % Debug dump on a failed solve.
   if ~ok && debug
      dumpIceEnbalFailure("maxiter_nonconvergence", ...
         T_ice, T_ice_old, T_iter, f_ice, f_liq, f_wat, k_eff, Sc, ...
         dt, T_sfc, iM, iter, maxiter, aN, aP, aS, b);
   end

   % Return the subsurface energy-balance linearization error [K].
   if nargout > 9
      err = icemodel.column.subsurface_linearization_error( ...
         T_ice, T_ice_old, f_ice, f_liq, f_liq_old, dro_vapdT, ...
         dHdT, Sc, q_deferred_faces, dz, dt, a2, Fc, Fp, a1);
   end

   % Surface energy-balance linearization error [K]. Currently disabled, retain
   % for reference.
   % surface_err = icemodel.column.surface_linearization_error( ...
   %    T_sfc, T_ice(1), Fc, Fp, a1); %#ok<NASGU>
end

function dumpIceEnbalFailure(reason, T_ice, T_ice_old, T_iter, f_ice, f_liq, ...
      f_wat, k_eff, Sc, dt, T_sfc, iM, iter, maxiter, aN, aP, aS, b)
   %DUMPICEENBALFAILURE Save enthalpy-solver failure diagnostics on demand.

   debug_file = getenv('ICEMODEL_DEBUG_ICEENBAL_FILE');
   if isempty(debug_file)
      return
   end

   % Pack the failed solve's inputs, iterates, and matrix rows.
   debug_state = struct();
   debug_state.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   debug_state.reason = reason;
   debug_state.iter = iter;
   debug_state.maxiter = maxiter;
   debug_state.dt = dt;
   debug_state.T_sfc = T_sfc;
   debug_state.max_abs_dT = max(abs(T_ice - T_iter));
   debug_state.iM = iM;
   debug_state.T_ice = T_ice;
   debug_state.T_ice_old = T_ice_old;
   debug_state.T_iter = T_iter;
   debug_state.f_ice = f_ice;
   debug_state.f_liq = f_liq;
   debug_state.f_wat = f_wat;
   debug_state.k_eff = k_eff;
   debug_state.Sc = Sc;
   debug_state.aN = aN;
   debug_state.aP = aP;
   debug_state.aS = aS;
   debug_state.b = b;

   save(debug_file, 'debug_state');
end
