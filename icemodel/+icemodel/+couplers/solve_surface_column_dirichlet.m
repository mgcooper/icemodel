function [T_sfc, T_ice, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok, n_iters] = ...
      solve_surface_column_dirichlet(xT_sfc, xT_ice, xf_ice, xf_liq, dz, delz, ...
      fn, Sc, Sp, dt, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ...
      ea_atm, br_coefs, ro_air_Lv, liqflag, chi, solver, tol, maxiter, alpha, ...
      use_aitken, jumpmax, cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, ...
      cpl_alpha, cpl_aitken, cpl_jumpmax, ro_sfc, snow_depth, opts)
   %SOLVE_SURFACE_COLUMN_DIRICHLET Coupled icemodel Dirichlet SEB solve.
   %
   % Run an outer Ts-T Picard loop so the accepted Dirichlet surface state,
   % top-node temperature, and conductive closure are mutually consistent at
   % the end of a substep. Ts remains the internal solver boundary state here;
   % physical diagnosed fluxes apply the accepted physical-surface
   % contract downstream via icemodel.surface.physical_surface_temperature
   % handled by icemodel.surface.diagnose_surface_fluxes.
   %
   %#codegen

   % Cache zero-valued Robin linearization values, Dirichlet bc used here.
   persistent Fc Fp
   if isempty(Fc)
      Fc = 0.0;
      Fp = 0.0;
   end

   debug = opts.debug;

   % Pre-coupler Ts predictor using checkpoint state.
   k_eff = icemodel.column.bulk_thermal_conductivity(xT_ice, xf_ice, xf_liq);
   [T_sfc, ok_seb] = icemodel.surface.solve_surface_energy_balance(xT_sfc, ...
      tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ...
      ro_air_Lv, liqflag, chi, xT_ice, k_eff, dz, ro_sfc, snow_depth, opts);

   % Initial past Picard iterates for Aitken acceleration.
   Ts_1 = nan;
   Ts_2 = nan;

   % Initial values for convergence checks.
   ok_ieb = false;
   ok_cpl = false;
   Ts_old = T_sfc;
   Ts_diag = T_sfc;
   seb_res = nan;
   n_iters = 0;
   T_ice = xT_ice;
   f_ice = xf_ice;
   f_liq = xf_liq;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling).
   for cpliter = 1:cpl_maxiter

      % Inner subsurface solve from checkpoint state using the trial Ts.
      [T_ice, f_ice, f_liq, k_eff, ok_ieb, n_iters] = ...
         icemodel.column.solve_column_enthalpy(T_sfc, xT_ice, xf_ice, ...
         xf_liq, Fc, Fp, Sc, Sp, dz, delz, fn, dt, solver, tol, maxiter, ...
         alpha, use_aitken, jumpmax, debug);

      if ~ok_ieb
         if debug
            dumpIceEbSolveDirichletFailure("iceenbal_failed", T_sfc, Ts_diag, ...
               Ts_old, T_ice, f_ice, f_liq, k_eff, Sc, dt, cpliter, ...
               cpl_maxiter, cpl_Ts_tol, cpl_seb_tol, seb_res, ok_seb, ...
               ok_ieb, ok_cpl, n_iters);
         end
         break
      end

      % Update surface density for the surface turbulent heat flux scheme.
      ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));

      % Inner surface solve using the updated trial state.
      Ts_old = T_sfc;
      [T_sfc, ok_seb] = icemodel.surface.solve_surface_energy_balance(T_sfc, ...
         tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ...
         ro_air_Lv, liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts);
      Ts_diag = T_sfc;

      if ~ok_seb
         if debug
            dumpIceEbSolveDirichletFailure("sebsolve_failed", T_sfc, Ts_diag, ...
               Ts_old, T_ice, f_ice, f_liq, k_eff, Sc, dt, cpliter, ...
               cpl_maxiter, cpl_Ts_tol, cpl_seb_tol, seb_res, ok_seb, ...
               ok_ieb, ok_cpl, n_iters);
         end
         break
      end

      % Diagnose SEB residual using the updated conductive state.
      seb_res = abs( ...
         icemodel.surface.surface_energy_balance_residual(T_sfc, tair, swd, ...
         lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, ...
         liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts) ...
         );

      % Check convergence (bypass coupler if cpl_maxiter == 1).
      if (cpl_maxiter == 1) || ...
            (abs(T_sfc - Ts_old) < cpl_Ts_tol && seb_res < cpl_seb_tol)
         ok_cpl = true;
         break
      end

      % Aitken acceleration with relaxation fallback.
      Ts_0 = T_sfc;
      T_sfc = icemodel.numerics.aitkenscalar(Ts_2, Ts_1, Ts_0, ...
         (1.0 - cpl_alpha) * Ts_old + cpl_alpha * T_sfc, ... % relaxation
         cpl_jumpmax, cpl_aitken);
      Ts_2 = Ts_1;
      Ts_1 = Ts_0;

      if abs(T_sfc - Ts_old) < cpl_Ts_tol && seb_res < cpl_seb_tol
         ok_cpl = true;
         break
      end
   end

   ok = ok_seb && ok_ieb && ok_cpl;

   if ~ok && debug
      dumpIceEbSolveDirichletFailure("coupler_nonconvergence", T_sfc, Ts_diag, ...
         Ts_old, T_ice, f_ice, f_liq, k_eff, Sc, dt, cpliter, cpl_maxiter, ...
         cpl_Ts_tol, cpl_seb_tol, seb_res, ok_seb, ok_ieb, ok_cpl, ...
         n_iters);
   end
end

function dumpIceEbSolveDirichletFailure(reason, Ts, Ts_diag, Ts_old, T, ...
      f_ice, f_liq, k_eff, Sc, dt, cpliter, cpl_maxiter, cpl_Ts_tol, ...
      cpl_seb_tol, seb_res, ok_seb, ok_ieb, ok_cpl, n_iters)
   %DUMPICEEBSOLVEDIRICHLETFAILURE Save coupled Dirichlet solver diagnostics.

   debug_file = getenv('ICEMODEL_DEBUG_ICEEBSOLVE_FILE');
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
   debug_state.Sc = Sc;
   debug_state.dt = dt;
   debug_state.cpliter = cpliter;
   debug_state.cpl_maxiter = cpl_maxiter;
   debug_state.cpl_Ts_tol = cpl_Ts_tol;
   debug_state.cpl_seb_tol = cpl_seb_tol;
   debug_state.seb_res = seb_res;
   debug_state.ok_seb = ok_seb;
   debug_state.ok_ieb = ok_ieb;
   debug_state.ok_cpl = ok_cpl;
   debug_state.n_iters = n_iters;

   save(debug_file, 'debug_state');
end
