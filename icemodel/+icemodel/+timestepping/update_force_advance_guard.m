function force_advance_streak_dt = update_force_advance_guard( ...
      force_advance_streak_dt, forced_advance, dt, dt_limit, timestep, ...
      numsteps, model_name)
   % Track persistent maxsubstep force-advance streaks across full steps.
   %
   %  streak_dt = icemodel.timestepping.update_force_advance_guard( ...
   %     streak_dt, forced_advance, dt, dt_limit, timestep, numsteps, ...
   %     model_name)
   %
   % Consecutive forced advances are allowed to span up to one full forcing
   % step. Beyond that, treat the run as broken and fail fast.
   %
   %#codegen

   if forced_advance
      force_advance_streak_dt = force_advance_streak_dt + dt;
      if force_advance_streak_dt > dt_limit + eps(dt_limit)
         error('icemodel:ForceAdvanceStreakExceeded', ...
            ['%s repeated checksubstep force-advance exceeded the allowed ', ...
            'streak at timestep %d/%d (streak_dt = %.0f s, limit = %.0f s).'], ...
            model_name, timestep, numsteps, force_advance_streak_dt, dt_limit);
      end
   else
      force_advance_streak_dt = 0.0;
   end
end
