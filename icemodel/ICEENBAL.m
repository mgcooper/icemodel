function [T, f_ice, f_liq, k_eff, ok, iter, a1, err] = ICEENBAL(T, f_ice, ...
      f_liq, dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, cv_liq, ro_ice, ...
      ro_liq, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, f_ell_min, f_ell_max, ...
      Fc, Fp, bc, tol, maxiter, ~, ~, ~, debug)
   %ICEENBAL Solve the ice energy balance.
   %
   % use_aitken and jumpmax are kept in the public contract so the older
   % thermal-solver option surface remains stable while node-wise Aitken
   % acceleration stays disabled here.
   %
   %#codegen

   % Update the water fraction
   f_wat = f_liq + f_ice * ro_ice / ro_liq;

   % Update the melt-zone boundaries. These are the volumetric liquid
   % fractions at T=TL and T=TH, given the current total water fraction.
   f_liq_min = f_wat .* f_ell_min;
   f_liq_max = f_wat .* f_ell_max;

   % Compute vapor density [kg m-3]
   ro_vap = VAPORDENSITY(T, f_liq);

   % Compute enthalpy [J m-3]
   H_old = TOTALHEAT(T, f_ice, f_liq, cv_ice, cv_liq, roLf, Ls * ro_vap, Tf);

   % Store past values
   T_old = T;
   f_liq_old = f_liq;

   % Initialize current values
   T_iter = T_old + 2 * tol;

   % Initial past Picard iterates for Aitken-acceleration (disabled)
   % T_1 = nan(size(T));
   % T_2 = nan(size(T));

   % Iterate to solve the nonlinear heat equation
   ok = true;
   for iter = 0:maxiter-1

      % Update vapor density and derivative [kg m-3, kg m-3 K-1]
      [ro_vap, drovdT] = VAPORDENSITY(T, f_liq);

      % Update vapor thermal diffusion coefficient [W m-1 K-1]
      k_vap = VAPORK(T, f_liq, drovdT);

      % Update bulk (effective) thermal conductivity
      k_eff = BULKTHERMALK(T, f_ice, f_liq, ro_ice, k_liq, k_vap);

      % Update total enthalpy
      H = TOTALHEAT(T, f_ice, f_liq, cv_ice, cv_liq, roLf, Ls * ro_vap, Tf);

      % Update the derivative of enthalpy wrt temperature
      dHdT = cv_ice * f_ice + cv_liq * f_liq;
      dLdT = 2.0 * fcp ^ 2.0 * (Tf - min(T, Tf)) .* f_wat ...
         ./ (1.0 + fcp ^ 2.0 * (Tf - min(T, Tf)) .^ 2.0) .^ 2.0;

      % Update the general equation coefficients
      [aN, aP, aS, b, iM, a1, a2] = GECOEFS(T, f_ice, f_liq, dHdT, dLdT, ...
         drovdT, H-H_old, Sc, k_eff, delz, fn, dz, dt, Ts, Ls, Lf, ro_liq, ...
         TL, JJ, Fc, Fp, bc);

      % % Check diagonal dominance and condition number
      % icemodel.checkdiags(aP, aN, aS)

      % Exit here so the state variables are updated on the final iteration
      if all(abs(T - T_iter) < tol)
         break
      end

      % Capture past values
      T_iter = T;

      % Solve the equation (predictor step)
      T = TRISOLVE(-aN, aP, -aS, b);

      % Update the temperature-enthalpy relationship (corrector step)
      [T, f_ice, f_liq, ok] = MZTRANSFORM(T, T_iter, f_liq, f_wat, dLdT, ...
         ro_ice, ro_liq, Tf, TL, TH, fcp, f_liq_min, f_liq_max, iM, ok, debug);

      % If failure, return to the main program and shorten the timestep
      if ~ok
         if debug
            dumpIceEnbalFailure("mztransform_rejected_state", T, T_old, ...
               T_iter, f_ice, f_liq, f_wat, k_eff, Sc, dt, Ts, iM, iter, ...
               maxiter, aN, aP, aS, b);
         end
         return
      end

      % Mass conservation check
      assertF(@() all(f_ice + f_liq * ro_liq / ro_ice <= 1 + eps))

      % Relaxation and Aitken (not implemented). Proper implementation requires
      % a MELTCURVE/MZTRANSFORM-consistent update of T, f_ice, and f_liq.
      % T = min(Tf, alpha * T + (1 - alpha) * T_iter);
      % if use_aitken
      %    T_0 = T;
      %    for mm = 1:numel(T)
      %       T(mm) = aitkenscalar(T_2(mm), T_1(mm), T_0(mm), T(mm), ...
      %          jumpmax);
      %    end
      %    T_2 = T_1;
      %    T_1 = T_0;
      % end
   end

   ok = iter < maxiter;

   if ~ok && debug
      dumpIceEnbalFailure("maxiter_nonconvergence", T, T_old, T_iter, ...
         f_ice, f_liq, f_wat, k_eff, Sc, dt, Ts, iM, iter, maxiter, ...
         aN, aP, aS, b);
   end

   % Surface energy balance linearization error [K]
   % err0 = -(Fc + Fp * Ts) / a1 - (T(1) - Ts);

   % Subsurface energy balance linearization error [K]
   err = (dt/dz(1) ...
      * (a2 * (T(2) - T(1)) ...
      + Fc + Fp * (Fc + a1 * T(1)) / (a1 - Fp) ...
      + Sc(1) * dz(1)) ...
      - ro_liq * Lf * (f_liq(1) - f_liq_old(1)) ) ...
      / (dHdT(1) + Ls * drovdT(1) * (1 - f_ice(1) - f_liq(1)))...
      - (T(1) - T_old(1));
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
