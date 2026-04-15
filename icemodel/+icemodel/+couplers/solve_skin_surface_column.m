function [Ts, T, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok, n_iters] = ...
      solve_skin_surface_column(xTs, xT, xf_ice, xf_liq, dz, delz, fn, dt, ...
      tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, chi, ro_air_Lv, ...
      br_coefs, liqflag, ~, tol, maxiter, alpha, cpl_maxiter, cpl_Ts_tol, ...
      cpl_seb_tol, cpl_alpha, cpl_aitken, cpl_jumpmax, debug, ro_sfc, snow_depth, ...
      opts)
   %SOLVE_SKIN_SURFACE_COLUMN Coupled skin-subsurface Ts-T solve.
   % Skinmodel-specific predictor-corrector coupler: iterates between surface
   % energy balance (Ts) and subsurface enthalpy (T, f_ice, f_liq) until Ts
   % and the SEB residual converge to within the specified tolerances.
   %
   %#codegen

   persistent Tf
   if isempty(Tf)
      Tf = icemodel.physicalConstant('Tf');
   end

   % Pre-coupler Ts predictor using checkpoint state.
   k_eff = icemodel.column.bulk_thermal_conductivity(xT, xf_ice, xf_liq, 0);
   [Ts, ok_seb] = icemodel.surface.solve_surface_energy_balance(xTs, tair, ...
      swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, ...
      liqflag, chi, xT, k_eff, dz, ro_sfc, snow_depth, opts);
   Ts = icemodel.surface.physical_surface_temperature(Ts);

   % To reinstate vapor-aware thermal conductivity, replace call above with:
   % k_eff = icemodel.column.bulk_thermal_conductivity(xT, xf_ice, xf_liq, k_vap);

   % Initialize past Picard iterates for Aitken-acceleration
   Ts_1 = nan;
   Ts_2 = nan;

   % Initial values for convergence checks
   ok_cpl = false;
   Ts_old = Ts;
   Ts_diag = Ts;
   seb_res = nan;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling).
   for cpliter = 1:cpl_maxiter

      % Inner subsurface solve from checkpoint state w/o physical advancement
      [T, f_ice, f_liq, k_eff, ok_ieb, n_iters] = ...
         icemodel.column.solve_column_temperature(xT, xf_ice, xf_liq, ...
         dz, delz, fn, dt, Ts, tol, maxiter, alpha, debug);

      if ~ok_ieb
         if debug
            dumpSkinEbSolveFailure("skinsolve_failed", Ts, Ts_diag, Ts_old, ...
               T, f_ice, f_liq, k_eff, dt, cpliter, cpl_maxiter, cpl_Ts_tol, ...
               cpl_seb_tol, seb_res, n_iters, ok_seb, ok_ieb, ok_cpl);
         end
         break
      end

      % Update surface density for the surface turbulent heat flux scheme.
      ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));

      % Inner surface solve (in-loop corrector using updated trial state)
      Ts_old = Ts;
      [Ts, ok_seb] = icemodel.surface.solve_surface_energy_balance( ...
         Ts, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, ...
         br_coefs, ro_air_Lv, liqflag, chi, T, k_eff, dz, ro_sfc, snow_depth, opts);
      Ts = icemodel.surface.physical_surface_temperature(Ts);
      Ts_diag = Ts;

      if not(ok_seb)
         if debug
            dumpSkinEbSolveFailure("sebsolve_failed", Ts, Ts_diag, Ts_old, T, ...
               f_ice, f_liq, k_eff, dt, cpliter, cpl_maxiter, cpl_Ts_tol, ...
               cpl_seb_tol, seb_res, n_iters, ok_seb, ok_ieb, ok_cpl);
         end
         break
      end

      % SEB residual
      seb_res = abs(icemodel.surface.surface_energy_balance_residual( ...
         Ts, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, ...
         br_coefs, ro_air_Lv, liqflag, chi, T, k_eff, dz, ...
         ro_sfc, snow_depth, opts));

      % At melt cap (Ts ~= Tf), positive residual is melt energy (Qm).
      if Ts >= Tf
         seb_res = 0.0;
      end

      % Check convergence (bypass coupler if cpl_maxiter == 1)
      if (cpl_maxiter == 1) || ...
            abs(Ts - Ts_old) < cpl_Ts_tol && seb_res < cpl_seb_tol
         ok_cpl = true;
         break
      end

      % Aitken acceleration with relaxation-fallback
      Ts_0 = Ts;
      Ts = icemodel.surface.physical_surface_temperature( ...
         icemodel.numerics.aitkenscalar(Ts_2, Ts_1, Ts_0, ...
         (1.0 - cpl_alpha) * Ts_old + cpl_alpha * Ts, ... % relaxation
         cpl_jumpmax, cpl_aitken) ...
         );
      Ts_2 = Ts_1;
      Ts_1 = Ts_0;

      if abs(Ts - Ts_old) < cpl_Ts_tol && seb_res < cpl_seb_tol
         ok_cpl = true;
         break
      end
   end

   % Hitting max coupling iterations without ok_cpl is a substep fail
   ok = ok_seb && ok_ieb && ok_cpl;

   if ~ok && debug
      dumpSkinEbSolveFailure("coupler_nonconvergence", Ts, Ts_diag, Ts_old, ...
         T, f_ice, f_liq, k_eff, dt, cpliter, cpl_maxiter, cpl_Ts_tol, ...
         cpl_seb_tol, seb_res, n_iters, ok_seb, ok_ieb, ok_cpl);
   end
end

function dumpSkinEbSolveFailure(reason, Ts, Ts_diag, Ts_old, T, f_ice, ...
      f_liq, k_eff, dt, cpliter, cpl_maxiter, cpl_Ts_tol, cpl_seb_tol, ...
      seb_res, n_iters, ok_seb, ok_ieb, ok_cpl)
   %DUMPSKINEBSOLVEFAILURE Save coupled skin-model solver diagnostics.

   debug_file = getenv('ICEMODEL_DEBUG_SKINEBSOLVE_FILE');
   if isempty(debug_file)
      return
   end

   debug_state = struct();
   debug_state.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   debug_state.reason = reason;
   debug_state.Ts = Ts;
   debug_state.Ts_diag = Ts_diag;
   debug_state.Ts_old = Ts_old;
   debug_state.T = T;
   debug_state.f_ice = f_ice;
   debug_state.f_liq = f_liq;
   debug_state.k_eff = k_eff;
   debug_state.dt = dt;
   debug_state.cpliter = cpliter;
   debug_state.cpl_maxiter = cpl_maxiter;
   debug_state.cpl_Ts_tol = cpl_Ts_tol;
   debug_state.cpl_seb_tol = cpl_seb_tol;
   debug_state.seb_res = seb_res;
   debug_state.n_iters = n_iters;
   debug_state.ok_seb = ok_seb;
   debug_state.ok_ieb = ok_ieb;
   debug_state.ok_cpl = ok_cpl;

   save(debug_file, 'debug_state');
end
