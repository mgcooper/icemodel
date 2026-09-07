function [T_sfc, T_ice, f_ice, f_liq, k_eff, diag] = ...
      solve_skin_surface_column(xT_sfc, xT_ice, xf_ice, xf_liq, dz, delz, ...
      fn, dt, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
      ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, ...
      chi, ro_sfc, snow_depth, settings, opts)
   %SOLVE_SKIN_SURFACE_COLUMN Coupled skin-subsurface Ts-T solve.
   %
   % Skinmodel-specific coupler. The coupler runs one outer Ts-T Picard loop.
   % The loop makes the accepted Dirichlet surface temperature (Ts), the
   % top-node temperature, and the conductive closure consistent with each other
   % at the end of a substep. solver 0 is the single-iteration special case
   % (cpl_maxiter = 1) solver 1 runs the full outer iterations. The coupler
   % makes one attempt. If it fails, icemodel.timestepping.checksubstep sets
   % recovery mode (underrelaxation and no acceleration).
   %
   % See also: icemodel.couplers.solve_surface_column_dirichlet,
   %  icemodel.couplers.solve_surface_column_robin,
   %  icemodel.couplers.initialize_solver_diag
   %
   %#codegen

   persistent Tf
   if isempty(Tf)
      Tf = icemodel.physicalConstant('Tf');
   end

   debug = settings.debug;

   % Pre-coupler Ts predictor using checkpoint state.
   k_eff = icemodel.column.bulk_thermal_conductivity(xT_ice, xf_ice, xf_liq, 0);

   % To reinstate vapor-aware thermal conductivity, replace call above with:
   % k_eff = ...
   %    icemodel.column.bulk_thermal_conductivity(xT, xf_ice, xf_liq, k_vap);

   [T_sfc, ok_seb] = icemodel.surface.solve_surface_energy_balance( ...
      xT_sfc, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
      ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, ...
      chi, xT_ice, k_eff, dz, ro_sfc, snow_depth, opts);
   T_sfc = icemodel.surface.physical_surface_temperature(T_sfc);

   % Initial solver histories for acceleration.
   hist = icemodel.couplers.initialize_coupler_history();

   % Load the default solve-attempt diagnostics.
   step_diag = icemodel.couplers.initialize_solver_diag();
   diag = step_diag.substep;
   ok_cpl = false;
   res_hist = diag.cpl_res_hist;

   % nan marks "no evaluated outer residual": an inner-solve failure
   % breaks out of the loop before the first residual evaluation below.
   cpl_res = nan;
   Ts_old = T_sfc;
   Ts_diag = T_sfc;
   seb_res = nan;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling).
   for cpliter = 1:settings.cpl_maxiter

      % Inner subsurface solve from checkpoint state w/o physical advancement.
      [T_ice, f_ice, f_liq, k_eff, ok_ieb, n_iters] = ...
         icemodel.column.solve_column_temperature(T_sfc, xT_ice, xf_ice, ...
         xf_liq, dz, delz, fn, dt, settings.tol, settings.maxiter, ...
         settings.alpha, debug);

      if ~ok_ieb
         if debug
            dumpSkinEbSolveFailure("skinsolve_failed", T_sfc, Ts_diag, ...
               Ts_old, T_ice, f_ice, f_liq, k_eff, dt, cpliter, ...
               settings, seb_res, n_iters, ok_seb, ok_ieb, ok_cpl, res_hist);
         end
         break
      end

      % Update surface density for the surface turbulent heat flux scheme.
      ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));

      % Inner surface solve (in-loop corrector using updated trial state).
      Ts_old = T_sfc;
      [T_sfc, ok_seb] = ...
         icemodel.surface.solve_surface_energy_balance(T_sfc, tair, ...
         swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ro_atm, ...
         cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, ...
         chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts);
      T_sfc = icemodel.surface.physical_surface_temperature(T_sfc);
      Ts_diag = T_sfc;

      % Debug dump and break on surface solve failure.
      if not(ok_seb)
         if debug
            dumpSkinEbSolveFailure("sebsolve_failed", T_sfc, Ts_diag, ...
               Ts_old, T_ice, f_ice, f_liq, k_eff, dt, cpliter, ...
               settings, seb_res, n_iters, ok_seb, ok_ieb, ok_cpl, res_hist);
         end
         break
      end

      % SEB residual
      seb_res = abs( ...
         icemodel.surface.surface_energy_balance_residual(T_sfc, ...
         tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
         ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, ...
         liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts) ...
         );

      % At melt cap (Ts ~= Tf), positive residual is melt energy (Qm).
      if T_sfc >= Tf
         seb_res = 0.0;
      end

      % Record the evaluated residual before acceleration replaces Ts
      % with an unevaluated trial, so a nonconvergent exit reports the
      % last evaluated iterate.
      cpl_res = T_sfc - Ts_old;
      res_hist = [res_hist(2:end); cpl_res];

      % Check convergence (bypass coupler if cpl_maxiter == 1).
      if (settings.cpl_maxiter == 1) || ...
            abs(cpl_res) < settings.cpl_Ts_tol ...
            && seb_res < settings.cpl_seb_tol
         ok_cpl = true;
         break
      end

      % Apply hybrid Aitken-secant step acceleration.
      [Ts_accel, hist] = icemodel.couplers.accelerate_coupler_iterate( ...
         hist, Ts_old, T_sfc, settings.cpl_alpha, settings.cpl_jumpmax, ...
         settings.cpl_aitken);

      % Apply the physical surface temperature.
      T_sfc = icemodel.surface.physical_surface_temperature(Ts_accel);

   end

   % Dump the outer failure only when neither inner dump ran. Both inner dumps
   % write the same debug file. Without this guard, the outer snapshot
   % overwrites the inner-solver state.
   if debug && ok_seb && ok_ieb && ~ok_cpl
      dumpSkinEbSolveFailure("coupler_nonconvergence", T_sfc, Ts_diag, ...
         Ts_old, T_ice, f_ice, f_liq, k_eff, dt, cpliter, settings, ...
         seb_res, n_iters, ok_seb, ok_ieb, ok_cpl, res_hist);
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

function dumpSkinEbSolveFailure(reason, Ts, Ts_diag, Ts_old, T, f_ice, ...
      f_liq, k_eff, dt, cpliter, settings, seb_res, n_iters, ok_seb, ...
      ok_ieb, ok_cpl, res_hist)
   %DUMPSKINEBSOLVEFAILURE Save coupled skin-model solver diagnostics.

   debug_file = getenv('ICEMODEL_DEBUG_SKINEBSOLVE_FILE');
   if isempty(debug_file)
      return
   end

   % Sequence-number the dump files so one run keeps every failure event.
   % A single overwritten file keeps only the last event, which cannot
   % show whether failures cluster or how the residual trajectory evolves
   % across events.
   persistent seq
   if isempty(seq)
      seq = 0;
   end
   seq = seq + 1;
   [dump_dir, dump_name, dump_ext] = fileparts(debug_file);
   debug_file = fullfile(dump_dir, ...
      sprintf('%s_%04d%s', dump_name, seq, dump_ext));

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
   debug_state.settings = settings;
   debug_state.seb_res = seb_res;
   debug_state.n_iters = n_iters;
   debug_state.ok_seb = ok_seb;
   debug_state.ok_ieb = ok_ieb;
   debug_state.ok_cpl = ok_cpl;
   debug_state.res_hist = res_hist;

   save(debug_file, 'debug_state');
end
