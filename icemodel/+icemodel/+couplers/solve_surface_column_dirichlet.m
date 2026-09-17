function [T_sfc, T_ice, f_ice, f_liq, k_eff, U_vap, L_vap, diag] = ...
      solve_surface_column_dirichlet( ...
      xT_sfc, xT_ice, xf_ice, xf_liq, Sc, Sp, dz, delz, fn, dt, tair, ...
      swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ro_atm, cv_atm, ...
      nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, chi, ro_sfc, snow_depth, ...
      f_res_por, settings, opts)
   %SOLVE_SURFACE_COLUMN_DIRICHLET Coupled surface-subsurface solve for
   %Dirichlet-type surface boundary condition.
   %
   % The coupler runs one outer T_sfc-T_ice Picard loop. The loop makes the
   % accepted Dirichlet surface temperature (T_sfc), the top-node temperature,
   % and the conductive closure consistent with each other at the end of a
   % substep. solver 0 is the single-iteration special case (cpl_maxiter = 1),
   % solver 1 runs the full outer iterations. The coupler makes one attempt. If
   % it fails, icemodel.timestepping.checksubstep sets recovery mode
   % (underrelaxation and no acceleration).
   %
   % U_VAP is the accepted face vapor mass flux and L_VAP is its face donor
   % latent heat (see icemodel.column.vapor_transport_terms).
   %
   % See also: icemodel,
   %  icemodel.couplers.solve_surface_column_robin
   %
   %#codegen

   debug = settings.debug;

   % Zero-valued Robin linearization values (Dirichlet bc used here).
   Fc = 0.0;
   Fp = 0.0;

   % Pre-coupler T_sfc predictor using checkpoint state. The predictor uses the
   % vapor-free node conductivity, like solve_column_enthalpy (vapor heat is in
   % the face terms).
   k_eff = icemodel.column.bulk_thermal_conductivity( ...
      xT_ice, xf_ice, xf_liq, 0);
   [T_sfc, ok_seb] = icemodel.surface.solve_surface_energy_balance( ...
      xT_sfc, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ...
      ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, ...
      liqflag, chi, xT_ice, k_eff, dz, ro_sfc, snow_depth, opts);

   % Initial solver histories for acceleration.
   hist = icemodel.couplers.initialize_coupler_history();

   % Initial values for convergence checks and the diag record.
   [~, diag] = icemodel.couplers.initialize_solver_diag();
   ok_cpl = false;
   seb_res = nan;
   cpl_res = nan;
   T_sfc_old = T_sfc;
   T_sfc_diag = T_sfc;

   % Signed outer-residual ring, newest last, sized to diag.cpl_res_hist. The
   % ring distinguishes bounded oscillation (alternating signs, stable
   % magnitude) from accelerator overshoot (growing magnitude).
   res_hist = diag.cpl_res_hist;

   % Run outer T_sfc-T_ice convergence loop (iterative block/Picard coupling).
   for cpliter = 1:settings.cpl_maxiter

      % Inner subsurface solve using the trial T_sfc and checkpoint state.
      [T_ice, f_ice, f_liq, k_eff, U_vap, L_vap, ok_ieb, n_iters] = ...
         icemodel.column.solve_column_enthalpy(T_sfc, xT_ice, xf_ice, ...
         xf_liq, Fc, Fp, Sc, Sp, dz, delz, fn, dt, f_res_por, settings);

      % Debug dump and break on subsurface solve failure.
      if ~ok_ieb
         if debug
            dumpIceEbSolveDirichletFailure("iceenbal_failed", ...
               T_sfc, T_sfc_diag, T_sfc_old, T_ice, f_ice, f_liq, k_eff, ...
               Sc, dt, cpliter, settings, seb_res, ok_seb, ok_ieb, ok_cpl, ...
               n_iters, res_hist);
         end
         break
      end

      % Update surface density for the surface turbulent heat flux scheme.
      ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));

      % Inner surface solve using the updated trial state.
      T_sfc_old = T_sfc;
      [T_sfc, ok_seb] = icemodel.surface.solve_surface_energy_balance( ...
         T_sfc, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
         ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, ...
         chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts);
      T_sfc_diag = T_sfc;

      % Debug dump and break on surface solve failure.
      if ~ok_seb
         if debug
            dumpIceEbSolveDirichletFailure("sebsolve_failed", ...
               T_sfc, T_sfc_diag, T_sfc_old, T_ice, f_ice, f_liq, k_eff, ...
               Sc, dt, cpliter, settings, seb_res, ok_seb, ok_ieb, ok_cpl, ...
               n_iters, res_hist);
         end
         break
      end

      % Diagnose SEB residual using the updated conductive state.
      seb_res = abs( ...
         icemodel.surface.surface_energy_balance_residual(T_sfc, ...
         tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
         ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, ...
         liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, ...
         opts) ...
         );

      % Check convergence (bypass coupler if cpl_maxiter == 1).
      cpl_res = T_sfc - T_sfc_old;
      res_hist = [res_hist(2:end); cpl_res];
      if (settings.cpl_maxiter == 1) || ...
            (abs(cpl_res) < settings.cpl_Ts_tol ...
            && seb_res < settings.cpl_seb_tol)
         ok_cpl = true;
         break
      end

      % Apply hybrid Aitken-secant step acceleration.
      [T_sfc, hist] = icemodel.couplers.accelerate_coupler_iterate( ...
         hist, T_sfc_old, T_sfc, settings.cpl_alpha, settings.cpl_jumpmax, ...
         settings.cpl_aitken);

   end

   % Dump the outer failure only when neither inner dump ran. Both inner
   % dumps write the same debug file, so the outer dump would overwrite
   % the inner-solver state.
   if debug && ok_seb && ok_ieb && ~ok_cpl
      dumpIceEbSolveDirichletFailure("coupler_nonconvergence", ...
         T_sfc, T_sfc_diag, T_sfc_old, T_ice, f_ice, f_liq, k_eff, Sc, dt, ...
         cpliter, settings, seb_res, ok_seb, ok_ieb, ok_cpl, n_iters, res_hist);
   end

   % Assemble the diag record.
   diag.ok_seb = ok_seb;
   diag.ok_ieb = ok_ieb;
   diag.ok_cpl = ok_cpl;
   diag.n_iters = n_iters;
   diag.cpl_iters = cpliter;
   diag.cpl_res = cpl_res;
   diag.seb_res = seb_res;
   diag.cpl_res_hist = res_hist;
end

function dumpIceEbSolveDirichletFailure(reason, T_sfc, T_sfc_diag, ...
      T_sfc_old, T_ice, f_ice, f_liq, k_eff, Sc, dt, cpliter, settings, ...
      seb_res, ok_seb, ok_ieb, ok_cpl, n_iters, res_hist)
   %DUMPICEEBSOLVEDIRICHLETFAILURE Save coupled Dirichlet solver diagnostics.

   debug_file = getenv('ICEMODEL_DEBUG_ICEEBSOLVE_FILE');
   if isempty(debug_file)
      return
   end

   % Number the dump files so one run keeps every failure event, to see if
   % failures cluster or how the residual trajectory evolves across events.
   persistent seq
   if isempty(seq)
      seq = 0;
   end
   seq = seq + 1;

   % Both couplers write to ICEMODEL_DEBUG_ICEEBSOLVE_FILE with their own
   % persistent counters, so name their debug files by the coupler name.
   [dump_dir, dump_name, dump_ext] = fileparts(debug_file);
   debug_file = fullfile(dump_dir, ...
      sprintf('%s_dirichlet_%04d%s', dump_name, seq, dump_ext));

   debug_state = struct();
   debug_state.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   debug_state.reason = reason;
   debug_state.T_sfc = T_sfc;
   debug_state.T_sfc_diag = T_sfc_diag;
   debug_state.T_sfc_old = T_sfc_old;
   debug_state.T_ice = T_ice;
   debug_state.f_ice = f_ice;
   debug_state.f_liq = f_liq;
   debug_state.k_eff = k_eff;
   debug_state.Sc = Sc;
   debug_state.dt = dt;
   debug_state.cpliter = cpliter;
   debug_state.settings = settings;
   debug_state.seb_res = seb_res;
   debug_state.ok_seb = ok_seb;
   debug_state.ok_ieb = ok_ieb;
   debug_state.ok_cpl = ok_cpl;
   debug_state.n_iters = n_iters;
   debug_state.res_hist = res_hist;

   save(debug_file, 'debug_state');
end
