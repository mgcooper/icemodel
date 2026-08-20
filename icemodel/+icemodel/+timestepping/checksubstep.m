function [Ts, T, f_ice, f_liq, k_eff, n_subfail, substep, dt, ok, ...
      forced_advance, force_advance_streak_dt] = checksubstep(Ts, T, ...
      f_ice, f_liq, k_eff, xTs, xT, xf_ice, xf_liq, xk_eff, dt_sum, dt, ...
      dt_FULL_STEP, timestep, numsteps, substep, maxsubstep, n_subfail, ...
      debug, ok, force_advance_streak_dt, model_name)
   %CHECKSUBSTEP Accept, retry, or force-advance the current substep.
   %
   %  [Ts, T, f_ice, f_liq, k_eff, n_subfail, substep, dt, ok, ...
   %     forced_advance, force_advance_streak_dt] = ...
   %     icemodel.timestepping.checksubstep(...)
   %
   % On an accepted substep (OK true), assert the water bounds and reset
   % the forced-advance streak. On a rejected substep, restore the
   % x-prefixed checkpoint, shorten dt, and, after MAXSUBSTEP failures,
   % force-advance at the checkpoint (OK returns true with
   % FORCED_ADVANCE true so the caller advances time only).
   %
   % FORCE_ADVANCE_STREAK_DT tracks consecutive forced-advance time
   % across forcing steps. One full step (DT_FULL_STEP) of consecutive
   % forced time is the limit; above it this function errors, because a
   % run that solves nothing must stop rather than integrate checkpoints
   % forward. MODEL_NAME names the driver in that error.
   %
   %#codegen

   dt_min = dt_FULL_STEP / maxsubstep;
   forced_advance = false;

   if ok
      % Mass conservation / control volume check
      assertF(@() icemodel.column.assert_max_water(f_ice, f_liq));
   else
      % Preserve the rejected state for optional debug dumps after the
      % timestep reset has updated the accepted retry state.
      [Ts_fail, T_fail, f_ice_fail, f_liq_fail, k_eff_fail] = ...
         deal(Ts, T, f_ice, f_liq, k_eff);

      % Adapt the timestep (shorten dt and restart the substep on failure).
      [Ts, T, f_ice, f_liq, k_eff, n_subfail, substep, dt] ...
         = icemodel.timestepping.resetsubstep(xTs, xT, xf_ice, xf_liq, ...
         xk_eff, dt_FULL_STEP, substep, maxsubstep, n_subfail, dt_sum);

      if debug == true && dt <= dt_min + eps(dt_min) && n_subfail < maxsubstep
         fprintf('timestep = %d (%.2f%%), dt = %.0f (dt_min), ok = %s\n', ...
            timestep, 100*timestep/numsteps, dt, mat2str(ok))
      end

      % On maxsubstep, force advance using last accepted state
      if n_subfail == maxsubstep
         ok = true;
         forced_advance = true;
         fprintf('timestep = %d, n_subfail == maxsubstep\n', timestep)
      end

      % Dump the failing state when debug mode is enabled.
      if debug
         if forced_advance
            projected_force_advance_dt = force_advance_streak_dt + dt;
         else
            projected_force_advance_dt = force_advance_streak_dt;
         end
         dumpMaxsubstepDebugState(Ts_fail, T_fail, f_ice_fail, f_liq_fail, ...
            k_eff_fail, xTs, xT, xf_ice, xf_liq, xk_eff, dt_sum, dt, ...
            dt_FULL_STEP, timestep, numsteps, substep, maxsubstep, ...
            n_subfail, forced_advance, projected_force_advance_dt, ...
            dt_FULL_STEP);
      end
   end

   % Track consecutive forced-advance time across forcing steps. A genuine
   % acceptance resets the streak; a forced advance extends it and errors
   % above one full step of consecutive forced time.
   if ok
      if forced_advance
         force_advance_streak_dt = force_advance_streak_dt + dt;
         if force_advance_streak_dt > dt_FULL_STEP + eps(dt_FULL_STEP)
            error('icemodel:ForceAdvanceStreakExceeded', ...
               ['%s repeated checksubstep force-advance exceeded the ', ...
               'allowed streak at timestep %d/%d (streak_dt = %.0f s, ', ...
               'limit = %.0f s).'], model_name, timestep, numsteps, ...
               force_advance_streak_dt, dt_FULL_STEP);
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
   %DUMPMAXSUBSTEPDEBUGSTATE Save the failed and checkpoint states on demand.

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
