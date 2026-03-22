function [Ts, T, f_ice, f_liq, n_subfail, substep, dt, ok] = CHECKSUBSTEP( ...
      Ts, T, f_ice, f_liq, xTs, xT, xf_ice, xf_liq, ro_ice, ro_liq, dt_sum, ...
      dt, dt_FULL_STEP, timestep, numsteps, substep, maxsubstep, n_subfail, ...
      debug, eps, ok)

   dt_min = dt_FULL_STEP / maxsubstep;

   if ok
      if debug == true
         % Mass conservation / control volume check
         assertF(@() all(f_ice + f_liq * ro_liq / ro_ice <= 1 + eps))
      end
   else
      % Dump the failing state before reset if this failure will saturate
      % the substep budget or if an explicit debug dump was requested.
      if (n_subfail + 1 >= maxsubstep) || ...
            ~isempty(getenv('ICEMODEL_DEBUG_MAXSUBSTEP_FILE'))
         dumpMaxsubstepDebugState(Ts, T, f_ice, f_liq, ...
            xTs, xT, xf_ice, xf_liq, dt_sum, dt, dt_FULL_STEP, timestep, ...
            numsteps, substep, maxsubstep, n_subfail + 1);
      end

      % Adapt the timestep (shorten dt and restart the substep on failure).
      [Ts, T, f_ice, f_liq, n_subfail, substep, dt] ...
         = RESETSUBSTEP(xTs, xT, xf_ice, xf_liq, dt_FULL_STEP, ...
         substep, maxsubstep, n_subfail, dt_sum);

      if debug == true && dt <= dt_min + eps && n_subfail < maxsubstep
         fprintf('timestep = %d (%.2f%%), dt = %.0f (dt_min), ok = %s\n', ...
            timestep, 100*timestep/numsteps, dt, mat2str(ok))
      end

      % On maxsubstep, force advance using last accepted state
      if n_subfail == maxsubstep
         ok = true;
         fprintf('timestep = %d, n_subfail == maxsubstep\n', timestep)
      end
   end
end

function dumpMaxsubstepDebugState(Ts_fail, T_fail, f_ice_fail, f_liq_fail, ...
      xTs, xT, xf_ice, xf_liq, dt_sum, dt, dt_FULL_STEP, timestep, ...
      numsteps, substep, maxsubstep, n_subfail)
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
   debug_state.failed = struct('Ts', Ts_fail, 'T', T_fail, ...
      'f_ice', f_ice_fail, 'f_liq', f_liq_fail);
   debug_state.checkpoint = struct('Ts', xTs, 'T', xT, ...
      'f_ice', xf_ice, 'f_liq', xf_liq);

   save(debug_file, 'debug_state');
end
