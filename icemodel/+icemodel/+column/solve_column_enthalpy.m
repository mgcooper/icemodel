function [T, f_ice, f_liq, k_eff, ok, iter, a1, err] = solve_column_enthalpy( ...
      T_sfc, T, f_ice, f_liq, Fc, Fp, Sc, Sp, dz, delz, fn, dt, solver, ...
      tol, maxiter, ~, ~, ~, debug)
   %SOLVE_COLUMN_ENTHALPY Solve the column enthalpy balance.
   %
   % The alpha, use_aitken, and jumpmax inputs are kept in the function
   % signature so the older thermal-solver option surface remains stable while
   % node-wise Aitken acceleration stays disabled here.
   %
   %#codegen

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   % Update the water fraction
   f_wat = icemodel.water_fraction(f_ice, f_liq);

   % Canonical melt-zone bounds and liquid-fraction limits.
   [TL, TH, f_ell_min, f_ell_max] = icemodel.column.meltzone_bounds();

   % Update the melt-zone boundaries. These are the volumetric liquid
   % fractions at T=TL and T=TH, given the current total water fraction.
   f_liq_min = f_wat .* f_ell_min;
   f_liq_max = f_wat .* f_ell_max;

   % Compute vapor density [kg m-3]
   ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq);

   % Compute enthalpy [J m-3]
   H_old = icemodel.column.total_enthalpy(T, f_ice, f_liq, f_wat, ro_vap);

   % Store past values
   T_ice_old = T;
   f_liq_old = f_liq;

   % Initialize current values
   T_iter = T_ice_old + 2 * tol;

   % Initial past Picard iterates for Aitken-acceleration (disabled)
   % T_1 = nan(size(T));
   % T_2 = nan(size(T));

   % Iterate to solve the nonlinear heat equation
   ok = true;
   for iter = 0:maxiter-1

      % Update vapor density and derivative [kg m-3, kg m-3 K-1]
      [ro_vap, dro_vapdT] = icemodel.vapor.saturation_vapor_density( ...
         T, f_liq);

      % Update vapor thermal diffusion coefficient [W m-1 K-1]
      k_vap = icemodel.vapor.vapor_thermal_diffusion_coefficient( ...
         T, f_liq, dro_vapdT);

      % Update bulk (effective) thermal conductivity
      k_eff = icemodel.column.bulk_thermal_conductivity( ...
         T, f_ice, f_liq, k_vap);

      % Update total enthalpy and derivative wrt temperature
      [H, dHdT, dLdT] = icemodel.column.total_enthalpy( ...
         T, f_ice, f_liq, f_wat, ro_vap);

      % Update the general equation coefficients
      [aN, aP, aS, b, iM, a1, a2] = icemodel.column.assemble_enthalpy_system( ...
         T, f_ice, f_liq, dHdT, dLdT, dro_vapdT, H - H_old, Sc, Sp, ...
         k_eff, delz, fn, dz, dt, TL, TH, T_sfc, Fc, Fp, solver);

      % % Check diagonal dominance and condition number
      % icemodel.checkdiags(aP, aN, aS)

      % Exit here so the state variables are updated on the final iteration
      if all(abs(T - T_iter) < tol)
         break
      end

      % Capture past values
      T_iter = T;

      % Solve the equation (predictor step)
      T = icemodel.numerics.trisolve(-aN, aP, -aS, b);

      % Update the temperature-enthalpy relationship (corrector step)
      [T, f_ice, f_liq, ok] = icemodel.column.meltzone_transform( ...
         T, T_iter, f_liq, f_wat, dLdT, TL, TH, f_liq_min, f_liq_max, ...
         iM, ok, debug);

      % If failure, return to the main program and shorten the timestep
      if ~ok
         if debug
            dumpIceEnbalFailure("mztransform_rejected_state", T, ...
               T_ice_old, T_iter, f_ice, f_liq, f_wat, k_eff, Sc, dt, ...
               T_sfc, iM, iter, maxiter, aN, aP, aS, b);
         end
         return
      end

      % Mass conservation check
      assertF(@() all(icemodel.water_fraction(f_ice, f_liq) <= ...
         ro_ice / ro_liq + eps))

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

   ok = iter < maxiter;

   if ~ok && debug
      dumpIceEnbalFailure("maxiter_nonconvergence", T, T_ice_old, ...
         T_iter, f_ice, f_liq, f_wat, k_eff, Sc, dt, T_sfc, iM, iter, ...
         maxiter, aN, aP, aS, b);
   end

   % Surface energy balance linearization error [K].
   surface_err = icemodel.column.surface_linearization_error( ...
      T_sfc, T(1), Fc, Fp, a1); %#ok<NASGU>

   % Subsurface energy balance linearization error [K].
   err = icemodel.column.subsurface_linearization_error(T, T_ice_old, ...
      f_ice, f_liq, f_liq_old, dro_vapdT, dHdT, Sc, dz, dt, a2, Fc, Fp, a1);
end

function dumpIceEnbalFailure(reason, T, T_old, T_iter, f_ice, f_liq, ...
      f_wat, k_eff, Sc, dt, Ts, iM, iter, maxiter, aN, aP, aS, b)
   %DUMPICEENBALFAILURE Save enthalpy-solver failure diagnostics on demand.

   debug_file = getenv('ICEMODEL_DEBUG_ICEENBAL_FILE');
   if isempty(debug_file)
      return
   end

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
