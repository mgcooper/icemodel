function [T_sfc, T_ice, f_ice, f_liq, U_vap, L_vap, k_eff, ok_seb, ...
      ok_ieb, ok_cpl, diag] = solve_surface_column_dirichlet(xT_sfc, xT_ice, ...
      xf_ice, xf_liq, Sc, Sp, dz, delz, fn, dt, tair, swd, lwd, albedo, ...
      wspd, ppt, tppt, psfc, ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, ...
      hv_atm, br_coefs, liqflag, chi, solver, tol, maxiter, alpha, ...
      use_aitken, jumpmax, cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, ...
      cpl_alpha, cpl_aitken, cpl_jumpmax, cpl_alpha_min, ro_sfc, ...
      snow_depth, f_res_por, opts)
   %SOLVE_SURFACE_COLUMN_DIRICHLET Coupled icemodel Dirichlet SEB solve.
   %
   % Run an outer Ts-T Picard loop. The loop makes the accepted Dirichlet
   % surface state, the top-node temperature, and the conductive closure
   % consistent with each other at the end of a substep. Ts is the internal
   % solver boundary state here. Downstream, the diagnosed fluxes use the
   % physical surface temperature from
   % icemodel.surface.physical_surface_temperature and
   % icemodel.surface.diagnose_surface_fluxes.
   % U_VAP and its face donor latent heat L_VAP follow the prognostic
   % phase state in the output list because the driver treats accepted
   % face transport as a first-class quantity.
   % DIAG is the fixed-schema observability struct from
   % icemodel.couplers.initialize_solver_diag; the ok flags stay plain
   % returns because the driver's control flow consumes them.
   %
   % CPL_ALPHA_MIN exists to keep both surface-column couplers on one
   % signature. This coupler has no recovery phase today: the known outer
   % limit cycle is a Robin-map behavior, and the Dirichlet solvers have
   % no recorded healthy-inner outer failure. A future Dirichlet recovery
   % would consume it the same way the Robin coupler does.
   %
   %#codegen

   % Zero-valued Robin linearization values, Dirichlet bc used here.
   Fc = 0.0;
   Fp = 0.0;

   debug = opts.debug;

   % Pre-coupler Ts predictor using checkpoint state.
   %
   % Vapor transport lives on the faces, and the surface face carries the
   % turbulent exchange. Leave vapor out of this node-wise predictor so it
   % matches solve_column_enthalpy. With solver 0, or cpl_maxiter of 1, this
   % predictor supplies the only sweep's boundary.
   k_eff = icemodel.column.bulk_thermal_conductivity( ...
      xT_ice, xf_ice, xf_liq, 0);
   [T_sfc, ok_seb] = icemodel.surface.solve_surface_energy_balance( ...
      xT_sfc, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ...
      ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, ...
      liqflag, chi, xT_ice, k_eff, dz, ro_sfc, snow_depth, opts);

   % Initial solver histories for acceleration.
   hist = icemodel.couplers.initialize_coupler_history();

   % Initial values for convergence checks and the observability record.
   diag = icemodel.couplers.initialize_solver_diag();
   diag.cpl_phase = 1;
   ok_cpl = false;
   Ts_old = T_sfc;
   Ts_diag = T_sfc;
   seb_res = nan;
   cpl_res = nan;

   % Signed outer-residual ring, newest last, shaped by the diag schema
   % so the record reads the same from both couplers.
   res_hist = diag.cpl_res_hist;

   % Hold the inner solve's accepted face flux and donor latent heat for
   % the production coupler.
   U_vap = zeros(numel(xT_ice) + 1, 1);
   L_vap = zeros(numel(xT_ice) + 1, 1);
   n_iters = nan;

   % Code generation requires assignment on every path; the loop below
   % reassigns cpliter on entry because cpl_maxiter is at least 1.
   cpliter = 0;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling).
   for cpliter = 1:cpl_maxiter

      % Inner subsurface solve from checkpoint state using the trial Ts.
      [T_ice, f_ice, f_liq, k_eff, ok_ieb, n_iters, ~, ~, U_vap, ...
         L_vap] = ...
         icemodel.column.solve_column_enthalpy(T_sfc, xT_ice, xf_ice, ...
         xf_liq, Fc, Fp, Sc, Sp, dz, delz, fn, dt, solver, tol, maxiter, ...
         alpha, use_aitken, jumpmax, debug, f_res_por);

      % Debug dump and break on subsurface solve failure.
      if ~ok_ieb
         if debug
            dumpIceEbSolveDirichletFailure( ...
               "iceenbal_failed", T_sfc, Ts_diag, Ts_old, T_ice, ...
               f_ice, f_liq, k_eff, Sc, dt, cpliter, cpl_maxiter, ...
               cpl_Ts_tol, cpl_seb_tol, seb_res, ok_seb, ok_ieb, ...
               ok_cpl, n_iters);
         end
         break
      end

      % Update surface density for the surface turbulent heat flux scheme.
      ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));

      % Inner surface solve using the updated trial state.
      Ts_old = T_sfc;
      [T_sfc, ok_seb] = icemodel.surface.solve_surface_energy_balance( ...
         T_sfc, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
         ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, ...
         chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts);
      Ts_diag = T_sfc;

      % Debug dump and break on surface solve failure.
      if ~ok_seb
         if debug
            dumpIceEbSolveDirichletFailure( ...
               "sebsolve_failed", T_sfc, Ts_diag, Ts_old, T_ice, ...
               f_ice, f_liq, k_eff, Sc, dt, cpliter, cpl_maxiter, ...
               cpl_Ts_tol, cpl_seb_tol, seb_res, ok_seb, ok_ieb, ...
               ok_cpl, n_iters);
         end
         break
      end

      % Diagnose SEB residual using the updated conductive state.
      seb_res = abs( ...
         icemodel.surface.surface_energy_balance_residual(T_sfc, ...
         tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
         ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, ...
         liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, ...
         opts));

      % Check convergence (bypass coupler if cpl_maxiter == 1).
      cpl_res = T_sfc - Ts_old;
      res_hist = [res_hist(2:end); cpl_res];
      if (cpl_maxiter == 1) || ...
            (abs(cpl_res) < cpl_Ts_tol && seb_res < cpl_seb_tol)
         ok_cpl = true;
         break
      end

      % Apply hybrid Aitken-secant step acceleration.
      [T_sfc, hist] = icemodel.couplers.accelerate_coupler_iterate( ...
         hist, Ts_old, T_sfc, cpl_alpha, cpl_jumpmax, cpl_aitken);

      % The accelerated iterate is not accepted here: its residual has not been
      % evaluated, and T_ice and k_eff belong to the sweep that produced the
      % pre-acceleration iterate. It is tested on the next sweep against its
      % own column solve.
   end

   % Dump the outer failure only when neither inner dump ran. Both inner
   % dumps write the same debug file, so the outer snapshot would overwrite
   % the inner-solver state.
   if debug && ok_seb && ok_ieb && ~ok_cpl
      dumpIceEbSolveDirichletFailure( ...
         "coupler_nonconvergence", T_sfc, Ts_diag, Ts_old, T_ice, ...
         f_ice, f_liq, k_eff, Sc, dt, cpliter, cpl_maxiter, ...
         cpl_Ts_tol, cpl_seb_tol, seb_res, ok_seb, ok_ieb, ...
         ok_cpl, n_iters);
   end

   % Assemble the observability record. This coupler has one phase, so
   % cpl_iters is the loop count and cpl_recovered stays false.
   diag.ok_seb = ok_seb;
   diag.ok_ieb = ok_ieb;
   diag.ok_cpl = ok_cpl;
   diag.n_iters = n_iters;
   diag.cpl_iters = cpliter;
   diag.cpl_res = cpl_res;
   diag.seb_res = seb_res;
   diag.cpl_res_hist = res_hist;
end

function dumpIceEbSolveDirichletFailure(reason, Ts, Ts_diag, Ts_old, ...
      T, f_ice, f_liq, k_eff, Sc, dt, cpliter, cpl_maxiter, ...
      cpl_Ts_tol, cpl_seb_tol, seb_res, ok_seb, ok_ieb, ok_cpl, ...
      n_iters)
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
