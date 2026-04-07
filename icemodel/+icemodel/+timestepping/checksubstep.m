function [Ts, T, f_ice, f_liq, n_subfail, substep, dt, ok, varargout] = ...
      checksubstep(Ts, T, f_ice, f_liq, xTs, xT, xf_ice, xf_liq, ro_ice, ...
      ro_liq, dt_sum, dt, dt_FULL_STEP, timestep, numsteps, substep, ...
      maxsubstep, n_subfail, debug, eps, ok, varargin)
   % Check accepted state, retry rejected substeps, and force advance on cap.

   dt_min = dt_FULL_STEP / maxsubstep;
   forced_advance = false;
   force_advance_streak_dt = NaN;
   force_advance_limit_dt = NaN;
   if ~isempty(varargin)
      force_advance_streak_dt = varargin{1};
   end
   if numel(varargin) >= 2
      force_advance_limit_dt = varargin{2};
   end

   if ok
      % Mass conservation / control volume check
      assertF(@() all(f_ice + f_liq * ro_liq / ro_ice <= 1 + eps))
   else
      % Preserve the rejected state for optional debug dumps after the
      % timestep reset has updated the accepted retry state.
      [Ts_fail, T_fail, f_ice_fail, f_liq_fail] = deal(Ts, T, f_ice, f_liq);

      % Adapt the timestep (shorten dt and restart the substep on failure).
      [Ts, T, f_ice, f_liq, n_subfail, substep, dt] ...
         = icemodel.timestepping.resetsubstep(xTs, xT, xf_ice, xf_liq, ...
         dt_FULL_STEP, substep, maxsubstep, n_subfail, dt_sum);

      if debug == true && dt <= dt_min + eps && n_subfail < maxsubstep
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
            xTs, xT, xf_ice, xf_liq, dt_sum, dt, dt_FULL_STEP, timestep, ...
            numsteps, substep, maxsubstep, n_subfail, forced_advance, ...
            projected_force_advance_dt, force_advance_limit_dt);
      end
   end

   if nargout >= 9
      varargout{1} = forced_advance;
   end
end

function dumpMaxsubstepDebugState(Ts_fail, T_fail, f_ice_fail, f_liq_fail, ...
      xTs, xT, xf_ice, xf_liq, dt_sum, dt, dt_FULL_STEP, timestep, ...
      numsteps, substep, maxsubstep, n_subfail, forced_advance, ...
      force_advance_streak_dt, force_advance_limit_dt)
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
      'f_ice', f_ice_fail, 'f_liq', f_liq_fail);
   debug_state.checkpoint = struct('Ts', xTs, 'T', xT, ...
      'f_ice', xf_ice, 'f_liq', xf_liq);

   save(debug_file, 'debug_state');
end
