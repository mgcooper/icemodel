function [T, f_ice, f_liq, k_eff, U_vap, L_vap, ok, iter, a1, err] = ...
      solve_column_enthalpy(T_sfc, T, f_ice, f_liq, Fc, Fp, Sc, Sp, dz, ...
      delz, fn, dt, solver, tol, maxiter, ~, ~, ~, debug, f_res_por)
   %SOLVE_COLUMN_ENTHALPY Solve the column enthalpy balance.
   %
   %  [T, f_ice, f_liq, k_eff, U_vap, L_vap, ok, iter, a1, err] = ...
   %     icemodel.column.solve_column_enthalpy(...)
   %
   % T, F_ICE, and F_LIQ are the accepted temperature and phase fractions.
   % K_EFF is the vapor-free node conductivity. U_VAP and L_VAP are the
   % accepted face vapor flux and donor latent heat. OK is false when the
   % melt-zone transform or iteration limit rejects the solve.
   %
   % U_VAP is positive downward [kg m-2 s-1] and is zero at both boundary
   % faces. L_VAP is Lv for a wet donor cell and Ls for a dry donor cell.
   % F_RES_POR sets the residual-liquid threshold used to select that phase.
   %
   % Aitken acceleration is disabled here. Keep alpha, use_aitken, and
   % jumpmax in the signature to match the thermal-solver.
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

   % Convert ice and liquid fractions to total water fraction.
   f_wat = icemodel.column.water_fraction(f_ice, f_liq);

   % Get the melt-zone liquid-fraction bounds at T=TL and T=TH.
   [f_liq_min, f_liq_max] = icemodel.column.meltzone_bounds(f_wat);

   % Compute vapor properties for the first iteration.
   [ro_vap, dro_vapdT] = icemodel.vapor.saturation_vapor_density(T, f_liq);
   [~, De] = icemodel.vapor.vapor_thermal_conductivity(T, f_liq, dro_vapdT);

   % Store the entering enthalpy [J m-3], used for the conservative solve.
   H_old = icemodel.column.bulk_enthalpy(T, f_ice, f_liq, f_wat, ro_vap);

   % Store the entering temperature and liquid fraction, used for the
   % subsurface linearization error estimate.
   T_ice_old = T;
   f_liq_old = f_liq;

   % Make the first convergence check fail.
   T_iter = T_ice_old + 2 * tol;

   % Initial past Picard iterates for Aitken-acceleration (disabled).
   % T_1 = nan(size(T));
   % T_2 = nan(size(T));

   % Initialize optional outputs.
   err = nan;
   U_vap = zeros(numel(T) + 1, 1);
   L_vap = Ls * ones(numel(T) + 1, 1);

   % Iterate to solve the nonlinear heat equation.
   ok = false;
   for iter = 0:maxiter-1

      % Build the vapor-free node conductivity and face vapor transport terms.
      k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);
      [k_eff_faces, ~, q_deferred_faces, U_vap, L_vap] = ...
         icemodel.column.vapor_transport_terms(T, f_ice, f_liq, k_eff, ...
         ro_vap, dro_vapdT, De, delz, fn, f_res_por);

      % Update enthalpy and its temperature derivatives.
      [H, dHdT, dFdT] = icemodel.column.bulk_enthalpy( ...
         T, f_ice, f_liq, f_wat, ro_vap, dro_vapdT);

      % Assemble the enthalpy equation.
      [aN, aP, aS, b, iM, a1, a2] = ...
         icemodel.column.assemble_enthalpy_system( ...
         T, f_ice, f_liq, dHdT, dFdT, dro_vapdT, H - H_old, Sc, Sp, ...
         k_eff_faces, delz, dz, dt, T_sfc, Fc, Fp, solver, ...
         q_deferred_faces);

      % % Check diagonal dominance and condition number.
      % icemodel.checkdiags(aP, aN, aS)

      % Stop here so the state is updated on the final iteration.
      if all(abs(T - T_iter) < tol)
         ok = true;
         break
      end

      % Save the current temperature before the next solve.
      T_iter = T;

      % Solve the linearized enthalpy equation.
      T = icemodel.numerics.trisolve(-aN, aP, -aS, b);

      % Convert enthalpy back to temperature and phase fractions.
      [T, f_ice, f_liq, ok_mz] = icemodel.column.meltzone_transform( ...
         T, T_iter, f_liq, f_wat, dFdT, f_liq_min, f_liq_max, iM, debug);

      % Return a rejected solve when the melt-zone transform fails.
      if ~ok_mz
         if debug
            dumpIceEnbalFailure("mztransform_rejected_state", T, ...
               T_ice_old, T_iter, f_ice, f_liq, f_wat, k_eff, Sc, dt, ...
               T_sfc, iM, iter, maxiter, aN, aP, aS, b);
         end
         return
      end

      % Check the maximum water fraction.
      assertF(@() icemodel.column.assert_max_water(f_ice, f_liq));

      % Update vapor properties for the next iteration.
      [ro_vap, dro_vapdT] = icemodel.vapor.saturation_vapor_density(T, f_liq);
      [~, De] = icemodel.vapor.vapor_thermal_conductivity(T, f_liq, dro_vapdT);

      % Relaxation and Aitken (not implemented). Proper implementation requires
      % a liquid_fraction_function/meltzone_transform-consistent update of
      % T, f_ice, and f_liq.
      % T = min(Tf, alpha * T + (1 - alpha) * T_iter);
      % if use_aitken
      %    T_0 = T;
      %    for mm = 1:numel(T)
      %       T(mm) = icemodel.numerics.aitkenscalar(T_2(mm), T_1(mm), ...
      %          T_0(mm), T(mm), ...
      %          jumpmax);
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
      dumpIceEnbalFailure("maxiter_nonconvergence", T, T_ice_old, ...
         T_iter, f_ice, f_liq, f_wat, k_eff, Sc, dt, T_sfc, iM, iter, ...
         maxiter, aN, aP, aS, b);
   end

   % Return the subsurface energy-balance linearization error [K].
   if nargout > 9
      err = icemodel.column.subsurface_linearization_error( ...
         T, T_ice_old, f_ice, f_liq, f_liq_old, dro_vapdT, ...
         dHdT, Sc, q_deferred_faces, dz, dt, a2, Fc, Fp, a1);
   end

   % Surface energy-balance linearization error [K]. Currently disabled, retain
   % for reference.
   % surface_err = icemodel.column.surface_linearization_error( ...
   %    T_sfc, T(1), Fc, Fp, a1); %#ok<NASGU>
end

function dumpIceEnbalFailure(reason, T, T_old, T_iter, f_ice, f_liq, ...
      f_wat, k_eff, Sc, dt, Ts, iM, iter, maxiter, aN, aP, aS, b)
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
   debug_state.Ts = Ts;
   debug_state.max_abs_dT = max(abs(T - T_iter));
   debug_state.iM = iM;
   debug_state.T = T;
   debug_state.T_old = T_old;
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
