function [T_sfc, T_ice, f_ice, f_liq, U_vap, L_vap, k_eff, ok_seb, ...
      ok_ieb, ok_cpl, diag] = solve_surface_column_robin(T_sfc, xT_ice, xf_ice, ...
      xf_liq, Sc, Sp, dz, delz, fn, dt, tair, swd, lwd, albedo, wspd, ppt, ...
      tppt, psfc, ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, ...
      br_coefs, liqflag, chi, solver, tol, maxiter, alpha, use_aitken, ...
      jumpmax, cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, cpl_alpha, ...
      cpl_aitken, cpl_jumpmax, cpl_alpha_min, ro_sfc, snow_depth, ...
      f_res_por, opts)
   %SOLVE_SURFACE_COLUMN_ROBIN Coupled icemodel Robin SEB solve.
   %
   % solver = 2 is the single-sweep special case of this Robin coupler
   % (cpl_maxiter = 1). solver = 3 runs the full outer Ts-T iterations.
   % U_VAP and its face donor latent heat L_VAP follow the prognostic
   % phase state in the output list, because the driver treats
   % accepted face transport as a first-class quantity. DIAG is the fixed-schema observability struct from
   % icemodel.couplers.initialize_solver_diag; the ok flags stay plain
   % returns because the driver's control flow consumes them.
   %
   % Recovery is default behavior, owned here. Phase 1 runs the caller's
   % (CPL_ALPHA, CPL_AITKEN). When the inner solve stays healthy but the
   % outer loop exhausts its iterations (ok_ieb && ~ok_cpl), phase 2 reruns
   % the same outer loop from the same entry state with the conservative
   % pair (CPL_ALPHA_MIN, false). The known failure is a bounded period-4
   % limit cycle of the accelerated map (NUK_U 2013 evidence); damped
   % relaxation converges it. Phase 2 self-suppresses when its policy
   % would repeat phase 1's, and never runs after an inner failure, so it
   % cannot mask one. Both phases restart the inner solve from the
   % checkpoint inputs (xT_ice, xf_ice, xf_liq), which this coupler always
   % receives, so no caller-side checkpoint re-entry exists. A failed
   % phase 2 returns ~ok_cpl and the driver's ordinary substep handling
   % (dt reduction, forced advance) proceeds unchanged.
   %
   %#codegen

   % The entry surface temperature doubles as the phase-2 restart state:
   % at every call site it equals the accepted checkpoint xT_sfc.
   Ts_entry = T_sfc;

   % Phase 1: the caller's primary policy.
   [T_sfc, T_ice, f_ice, f_liq, U_vap, L_vap, k_eff, ok_ieb, ok_cpl, ...
      diag] = robinOuterLoop(Ts_entry, xT_ice, xf_ice, xf_liq, Sc, Sp, dz, delz, ...
      fn, dt, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
      ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, chi, ...
      solver, tol, maxiter, alpha, use_aitken, jumpmax, cpl_Ts_tol, ...
      cpl_seb_tol, cpl_maxiter, cpl_alpha, cpl_aitken, cpl_jumpmax, ...
      ro_sfc, snow_depth, f_res_por, opts, 1);

   % Return ok_seb true to align coupler function signatures (no Ts solve here).
   ok_seb = true;

   % Phase 2: one conservative rerun of the same outer loop from the same
   % entry state. Suppressed when phase 1's policy already was the
   % conservative pair, because the rerun would repeat it bit for bit.
   if ok_ieb && ~ok_cpl && (cpl_aitken || cpl_alpha ~= cpl_alpha_min)
      diag_primary = diag;
      [T_sfc, T_ice, f_ice, f_liq, U_vap, L_vap, k_eff, ok_ieb, ok_cpl, ...
         diag] = robinOuterLoop(Ts_entry, xT_ice, xf_ice, xf_liq, Sc, Sp, dz, ...
         delz, fn, dt, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ...
         ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, ...
         liqflag, chi, solver, tol, maxiter, alpha, use_aitken, jumpmax, ...
         cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, cpl_alpha_min, false, ...
         cpl_jumpmax, ro_sfc, snow_depth, f_res_por, opts, 2);
      diag.cpl_iters = diag.cpl_iters + diag_primary.cpl_iters;
      diag.cpl_recovered = ok_ieb && ok_cpl;
   end
   diag.ok_seb = ok_seb;
end

function [T_sfc, T_ice, f_ice, f_liq, U_vap, L_vap, k_eff, ok_ieb, ...
      ok_cpl, diag] = robinOuterLoop(T_sfc, xT_ice, xf_ice, xf_liq, Sc, Sp, dz, ...
      delz, fn, dt, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ...
      ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, ...
      liqflag, chi, solver, tol, maxiter, alpha, use_aitken, jumpmax, ...
      cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, cpl_alpha, cpl_aitken, ...
      cpl_jumpmax, ro_sfc, snow_depth, f_res_por, opts, phase)
   %ROBINOUTERLOOP Run one outer Ts-T Picard loop under one policy.
   %
   % One phase of the two-phase solve above: the loop structure is the
   % production Robin coupling, unchanged; PHASE only labels the diag
   % record and the debug dumps.

   debug = opts.debug;

   % Initial values for SEB linearization coefficients Fc, Fp.
   [Fc, Fp] = icemodel.surface.surface_flux_linearization( ...
      T_sfc, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
      ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, ...
      chi, ro_sfc, snow_depth, opts);

   % Initial solver histories for acceleration.
   hist = icemodel.couplers.initialize_coupler_history();

   % Initial values for convergence checks and the observability record.
   diag = icemodel.couplers.initialize_solver_diag();
   diag.cpl_phase = phase;
   ok_cpl = false;
   Ts_old = T_sfc;
   Ts_diag = T_sfc;
   seb_res = nan;
   cpl_res = nan;

   % Signed outer-residual ring, newest last, shaped by the diag schema.
   % Sixteen entries cover the tail of any outer loop at the production
   % cpl_maxiter. The ring separates bounded oscillation (alternating
   % signs, stable magnitude) from accelerator overshoot (growing
   % magnitude), which one final iterate cannot. It rides the diag record
   % and the failure dumps.
   res_hist = diag.cpl_res_hist;

   % Hold the inner solve's accepted face flux and donor latent heat for
   % the production coupler.
   U_vap = zeros(numel(xT_ice) + 1, 1);
   L_vap = zeros(numel(xT_ice) + 1, 1);
   [T_ice, f_ice, f_liq] = deal(xT_ice, xf_ice, xf_liq);
   k_eff = zeros(numel(xT_ice), 1);
   n_iters = nan;

   % Code generation requires assignment on every path; the loop below
   % reassigns cpliter on entry because cpl_maxiter is at least 1.
   cpliter = 0;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling).
   for cpliter = 1:cpl_maxiter
      Ts_old = T_sfc;

      % Inner subsurface solve with updated Ts, Fc, Fp and checkpoint state.
      [T_ice, f_ice, f_liq, k_eff, ok_ieb, n_iters, a1, ~, U_vap, ...
         L_vap] = ...
         icemodel.column.solve_column_enthalpy(T_sfc, xT_ice, xf_ice, ...
         xf_liq, Fc, Fp, Sc, Sp, dz, delz, fn, dt, solver, tol, maxiter, ...
         alpha, use_aitken, jumpmax, debug, f_res_por);

      % Debug dump and break on subsurface solve failure.
      if ~ok_ieb
         if debug
            dumpIceEbSolveFailure("iceenbal_failed", T_sfc, ...
               Ts_diag, Ts_old, T_ice, f_ice, f_liq, k_eff, Sc, ...
               dt, Fc, Fp, cpliter, cpl_maxiter, cpl_Ts_tol, ...
               cpl_seb_tol, seb_res, ok_ieb, ok_cpl, n_iters, ...
               res_hist, phase);
         end
         break
      end

      % Update surface density for the surface turbulent heat flux scheme.
      ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));

      % Diagnose Ts from frozen coefficients and updated (a1, T1).
      T_sfc = (Fc + a1 * T_ice(1)) / (a1 - Fp);
      Ts_diag = T_sfc;

      % Diagnose SEB residual. Take the magnitude here, as the Dirichlet and
      % skin couplers do.
      seb_res = abs( ...
         icemodel.surface.surface_energy_balance_residual( ...
         T_sfc, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ...
         ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, ...
         br_coefs, liqflag, chi, T_ice, k_eff, dz, ro_sfc, ...
         snow_depth, opts));

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

      % Use the new T_sfc solution to update Fc, Fp for the next column solve.
      [Fc, Fp] = icemodel.surface.surface_flux_linearization( ...
         T_sfc, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ...
         ea_atm, ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, ...
         br_coefs, liqflag, chi, ro_sfc, snow_depth, opts);
   end

   % Dump outer-loop failure only when the inner-failure dump did not run.
   if ~(ok_ieb && ok_cpl) && debug && ok_ieb
      dumpIceEbSolveFailure("coupler_nonconvergence", T_sfc, ...
         Ts_diag, Ts_old, T_ice, f_ice, f_liq, k_eff, Sc, dt, ...
         Fc, Fp, cpliter, cpl_maxiter, cpl_Ts_tol, cpl_seb_tol, ...
         seb_res, ok_ieb, ok_cpl, n_iters, res_hist, phase);
   end

   % Assemble the observability record for this phase.
   diag.ok_ieb = ok_ieb;
   diag.ok_cpl = ok_cpl;
   diag.n_iters = n_iters;
   diag.cpl_iters = cpliter;
   diag.cpl_res = cpl_res;
   diag.seb_res = seb_res;
   diag.cpl_res_hist = res_hist;
end

function dumpIceEbSolveFailure(reason, Ts, Ts_diag, Ts_old, T, f_ice, ...
      f_liq, k_eff, Sc, dt, Fc, Fp, cpliter, cpl_maxiter, cpl_Ts_tol, ...
      cpl_seb_tol, seb_res, ok_ieb, ok_cpl, n_iters, res_hist, phase)
   %DUMPICEEBSOLVEFAILURE Save coupled ice-SEB solver diagnostics on demand.

   debug_file = getenv('ICEMODEL_DEBUG_ICEEBSOLVE_FILE');
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
   debug_state.phase = phase;
   debug_state.Ts = Ts;
   debug_state.Ts_diag = Ts_diag;
   debug_state.Ts_old = Ts_old;
   debug_state.T = T;
   debug_state.f_ice = f_ice;
   debug_state.f_liq = f_liq;
   debug_state.k_eff = k_eff;
   debug_state.Sc = Sc;
   debug_state.dt = dt;
   debug_state.Fc = Fc;
   debug_state.Fp = Fp;
   debug_state.cpliter = cpliter;
   debug_state.cpl_maxiter = cpl_maxiter;
   debug_state.cpl_Ts_tol = cpl_Ts_tol;
   debug_state.cpl_seb_tol = cpl_seb_tol;
   debug_state.seb_res = seb_res;
   debug_state.ok_ieb = ok_ieb;
   debug_state.ok_cpl = ok_cpl;
   debug_state.n_iters = n_iters;
   debug_state.res_hist = res_hist;

   save(debug_file, 'debug_state');
end
