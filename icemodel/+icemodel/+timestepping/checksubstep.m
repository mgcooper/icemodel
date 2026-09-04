function [Ts, T, f_ice, f_liq, k_eff, substep, dt, all_ok, forced_advance, ...
      force_advance_streak_dt, settings, diag] = checksubstep( ...
      Ts, T, f_ice, f_liq, k_eff, xTs, xT, xf_ice, xf_liq, xk_eff, ...
      dt_sum, dt, timestep, numsteps, substep, ...
      force_advance_streak_dt, model_name, settings, settings0, diag)
   %CHECKSUBSTEP Accept the substep or decide what the next attempt looks like.
   %
   %  [Ts, T, f_ice, f_liq, k_eff, substep, dt, ok, forced_advance, ...
   %     force_advance_streak_dt, settings, diag] = ...
   %     icemodel.timestepping.checksubstep(...)
   %
   % This function decides the next substep attempt after every solve. DIAG is
   % the forcing step's solver-diagnostics record; diag.substep holds the
   % attempt the coupler just returned. The substep is accepted when the
   % attempt's ok_seb, ok_ieb, and ok_cpl are all true.
   %
   % On an accepted substep, assert the cv max water check and reset the
   % forced-advance streak. icemodel.timestepping.acceptsubstep restores
   % the default solver settings after an accepted solve.
   %
   % On a rejected substep, restore the x-prefixed checkpoint and pick the
   % retry, in this order:
   %   1. Settings retry (coupler fail): Retry at the same dt using "recovery
   %      mode" outer-loop relaxation settings (acceleration off and
   %      underrelaxation via cpl_recovery_alpha). Recovery dampens an
   %      outer-loop (coupler) residual oscillation that a shorter dt cannot,
   %      whereas a failed SEB or inner solve goes directly to dt shortening
   %      (outer relaxation cannot repair those). Use this when ok_seb and
   %      ok_ieb are true, but ok_cpl is false. SETTINGS must not already be in
   %      recovery mode (repeating a failed solve in recovery mode with the same
   %      dt yields the same failed solve). Note: the recovery settings must
   %      differ from the primary ones for "not in recovery mode" check to work.
   %   2. Shorten dt (icemodel.timestepping.resetsubstep) for every other
   %      failure, and when the settings retry already ran. The primary
   %      SETTINGS0 is restored here: every point in the main control flow that
   %      resets dt restores the primary settings. Each dt shortening event adds
   %      +1 to diag.n_failed_substeps; the settings retry is not incremented.
   %   3. After settings.maxsubstep failures, force-advance at the
   %      checkpoint (OK returns true with FORCED_ADVANCE true so the caller
   %      advances time only) and add +1 to diag.n_forced_advances.
   %
   % FORCE_ADVANCE_STREAK_DT tracks consecutive forced-advance time across
   % forcing steps. One full step (settings.dt_full_step) of consecutive
   % forced time is the limit; above it this function errors, because a
   % run that solves nothing must stop rather than integrate checkpoints
   % forward. MODEL_NAME identifies the driver in the error dump.
   %
   % See also: icemodel, skinmodel
   %
   %#codegen

   % Compute the minimum allowable dt and initialize forced_advance false.
   dt_min = settings.dt_full_step / settings.maxsubstep;
   forced_advance = false;

   % Get the surface solve, subsurface solve, and coupler success flags.
   [ok_seb, ok_ieb, ok_cpl] = deal( ...
      diag.substep.ok_seb, diag.substep.ok_ieb, diag.substep.ok_cpl);

   % all_ok means the surface, subsurface, and coupler solves succeeded.
   % inner_ok means the surface and subsurface solves succeeded but the
   % coupler failed.
   all_ok = ok_seb && ok_ieb && ok_cpl;
   inner_ok = ok_seb && ok_ieb && ~ok_cpl;

   % Allow recovery mode when it wasn't just tried and acceleration is active or
   % the coupler relaxation factor exceeds the configured recovery mode value.
   allow_recovery_mode = ~settings.cpl_recovery_active ...
      && (settings.cpl_aitken ...
      || settings.cpl_alpha > settings.cpl_recovery_alpha);

   if all_ok
      % Mass conservation / control volume check
      assertF(@() icemodel.column.assert_max_water(f_ice, f_liq));

   elseif inner_ok && allow_recovery_mode
      % Restore the checkpoint and retry the same dt with under-relaxation and
      % acceleration disabled. Set cpl_recovery_active so this retry runs once.
      % A dt reset restores settings0 and clears the flag. Skip this retry when
      % the primary settings are already in recovery mode (repeating a failed
      % solve in recovery mode with the same dt yields the same failed solve).
      [Ts, T, f_ice, f_liq, k_eff] = deal(xTs, xT, xf_ice, xf_liq, xk_eff);
      settings.cpl_alpha = settings.cpl_recovery_alpha;
      settings.cpl_aitken = false;
      settings.cpl_recovery_active = true;
   else
      % Solver fail. Reset the prior accepted state, shorten the timestep, and
      % retry the substep.

      % Restore settings0 before a dt reset. The shorter retry then uses
      % the primary settings.
      settings = settings0;

      % Keep the failed state for the optional debug dump.
      [Ts_fail, T_fail, f_ice_fail, f_liq_fail, k_eff_fail] = ...
         deal(Ts, T, f_ice, f_liq, k_eff);

      % Reset state, shorten dt, and add +1 to the n_subfail failure count.
      [Ts, T, f_ice, f_liq, k_eff, n_subfail, substep, dt] ...
         = icemodel.timestepping.resetsubstep(xTs, xT, xf_ice, xf_liq, ...
         xk_eff, settings.dt_full_step, substep, settings.maxsubstep, ...
         diag.n_failed_substeps, dt_sum);
      diag.n_failed_substeps = n_subfail;

      % Print an optional debug message.
      if settings.debug && dt <= dt_min + eps(dt_min) ...
            && n_subfail < settings.maxsubstep
         fprintf('timestep = %d (%.2f%%), dt = %.0f (dt_min), ok = %s\n', ...
            timestep, 100*timestep/numsteps, dt, mat2str(all_ok))
      end

      % On maxsubstep failures, force advance using the last accepted state.
      if n_subfail == settings.maxsubstep
         all_ok = true;
         forced_advance = true;
         diag.n_forced_advances = diag.n_forced_advances + 1;
         fprintf('timestep = %d, n_subfail == maxsubstep\n', timestep)
      end

      % Dump the failing state when debug mode is enabled.
      if settings.debug
         if forced_advance
            projected_force_advance_dt = force_advance_streak_dt + dt;
         else
            projected_force_advance_dt = force_advance_streak_dt;
         end
         dumpMaxsubstepDebugState(Ts_fail, T_fail, f_ice_fail, f_liq_fail, ...
            k_eff_fail, xTs, xT, xf_ice, xf_liq, xk_eff, dt_sum, dt, ...
            settings.dt_full_step, timestep, numsteps, substep, ...
            settings.maxsubstep, n_subfail, forced_advance, ...
            projected_force_advance_dt, settings.dt_full_step);
      end
   end

   % Track consecutive forced-advance time across forcing steps. A genuine
   % acceptance resets the streak; a forced advance extends it and errors
   % when consecutive forced time exceeds one full step.
   if all_ok
      if forced_advance
         force_advance_streak_dt = force_advance_streak_dt + dt;
         if force_advance_streak_dt > ...
               settings.dt_full_step + eps(settings.dt_full_step)
            error('icemodel:ForceAdvanceStreakExceeded', ...
               ['%s repeated checksubstep force-advance exceeded the ', ...
               'allowed streak at timestep %d/%d (streak_dt = %.0f s, ', ...
               'limit = %.0f s).'], model_name, timestep, numsteps, ...
               force_advance_streak_dt, settings.dt_full_step);
         end
      else
         force_advance_streak_dt = 0.0;
      end
   end
end

function dumpMaxsubstepDebugState(Ts_fail, T_fail, f_ice_fail, f_liq_fail, ...
      k_eff_fail, xTs, xT, xf_ice, xf_liq, xk_eff, dt_sum, dt, ...
      dt_FULL_STEP, timestep, numsteps, substep, maxsubstep, n_subfail, ...
      forced_advance, force_advance_streak_dt, force_advance_limit_dt)
   %DUMPMAXSUBSTEPDEBUGSTATE Save the failed and checkpoint states.

   debug_file = getenv('ICEMODEL_DEBUG_MAXSUBSTEP_FILE');
   if isempty(debug_file)
      return
   end

   debug_state = struct();
   debug_state.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   debug_state.timestep = timestep;
   debug_state.numsteps = numsteps;
   debug_state.dt_sum = dt_sum;
   debug_state.dt = dt;
   debug_state.dt_full_step = dt_FULL_STEP;
   debug_state.dt_min = dt_FULL_STEP / maxsubstep;
   debug_state.substep = substep;
   debug_state.maxsubstep = maxsubstep;
   debug_state.n_subfail = n_subfail;
   debug_state.forced_advance = forced_advance;
   debug_state.force_advance_streak_dt = force_advance_streak_dt;
   debug_state.force_advance_limit_dt = force_advance_limit_dt;
   debug_state.failed = struct('Ts', Ts_fail, 'T', T_fail, ...
      'f_ice', f_ice_fail, 'f_liq', f_liq_fail, 'k_eff', k_eff_fail);
   debug_state.checkpoint = struct('Ts', xTs, 'T', xT, ...
      'f_ice', xf_ice, 'f_liq', xf_liq, 'k_eff', xk_eff);

   save(debug_file, 'debug_state');
end
