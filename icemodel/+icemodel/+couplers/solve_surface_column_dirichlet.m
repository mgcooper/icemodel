function [T_sfc, T_ice, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok_cpl, n_iters] = ...
      solve_surface_column_dirichlet(xT_sfc, xT_ice, xf_ice, ...
      xf_liq, Sc, Sp, dz, delz, fn, dt, tair, swd, lwd, albedo, wspd, ...
      ppt, tppt, psfc, ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, ...
      hv_atm, br_coefs, liqflag, chi, solver, tol, maxiter, alpha, ...
      use_aitken, jumpmax, cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, ...
      cpl_alpha, cpl_aitken, cpl_jumpmax, ro_sfc, snow_depth, ...
      f_res_por, opts)
   %SOLVE_SURFACE_COLUMN_DIRICHLET Coupled icemodel Dirichlet SEB solve.
   %
   % Run an outer Ts-T Picard loop. The loop makes the accepted Dirichlet
   % surface state, the top-node temperature, and the conductive closure
   % consistent with each other at the end of a substep. Ts is the internal
   % solver boundary state here. Downstream, the diagnosed fluxes use the
   % physical surface temperature from
   % icemodel.surface.physical_surface_temperature and
   % icemodel.surface.diagnose_surface_fluxes.
   %
   %#codegen

   % Cache zero-valued Robin linearization values, Dirichlet bc used here.
   persistent Fc Fp
   if isempty(Fc)
      Fc = 0.0;
      Fp = 0.0;
   end

   debug = opts.debug;

   % The coupled vapor mode moves the vapor term from the node conductivity
   % to the face conductance. It is opt-in and nothing sets it by default.
   % Guard the read so a caller-built opts struct without the field still
   % runs the default path.
   use_coupled_vapor = isfield(opts, 'use_coupled_vapor') ...
      && opts.use_coupled_vapor;

   % Pre-coupler Ts predictor using checkpoint state.
   %
   % Default mode: the three-argument form evaluates the vapor term inside
   % bulk_thermal_conductivity, so saturation_vapor_density runs here on the
   % same state solve_column_enthalpy evaluates it on twice more below: once
   % before its loop for H_old, and once on iteration 0. Bead icemodel-3xg
   % reduces the three to one, with the C3 retirement of the node vapor
   % term, because removing this one alone would change default results.
   %
   % Coupled mode: the vapor term lives on the faces, and the surface face
   % carries the turbulent exchange rather than a diffusive vapor term, so
   % the predictor k_eff must leave the node vapor term out. The explicit
   % zero matches what solve_column_enthalpy returns in coupled mode. With
   % solver 0, or any cpl_maxiter of 1, this predictor supplies the only
   % sweep's boundary, so mixing the two formulations here would put the
   % legacy vapor term into the accepted surface state.
   if use_coupled_vapor
      k_eff = icemodel.column.bulk_thermal_conductivity( ...
         xT_ice, xf_ice, xf_liq, 0);
   else
      k_eff = icemodel.column.bulk_thermal_conductivity( ...
         xT_ice, xf_ice, xf_liq);
   end
   [T_sfc, ok_seb] = icemodel.surface.solve_surface_energy_balance( ...
      xT_sfc, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ...
      ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, ...
      liqflag, chi, xT_ice, k_eff, dz, ro_sfc, snow_depth, opts);

   % Initial solver histories for acceleration.
   hist = icemodel.couplers.initialize_coupler_history();

   % Initial values for convergence checks.
   ok_cpl = false;
   Ts_old = T_sfc;
   Ts_diag = T_sfc;
   seb_res = nan;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling).
   for cpliter = 1:cpl_maxiter

      % Inner subsurface solve from checkpoint state using the trial Ts.
      [T_ice, f_ice, f_liq, k_eff, ok_ieb, n_iters] = ...
         icemodel.column.solve_column_enthalpy( ...
         T_sfc, xT_ice, xf_ice, xf_liq, Fc, Fp, Sc, Sp, dz, delz, ...
         fn, dt, solver, tol, maxiter, alpha, use_aitken, jumpmax, debug, ...
         use_coupled_vapor, f_res_por);

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
      [T_sfc, ok_seb] = ...
         icemodel.surface.solve_surface_energy_balance( ...
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
      if (cpl_maxiter == 1) || ...
            (abs(T_sfc - Ts_old) < cpl_Ts_tol && seb_res < cpl_seb_tol)
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
