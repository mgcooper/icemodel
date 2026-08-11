function force_advance_streak_dt = update_force_advance_guard( ...
      force_advance_streak_dt, forced_advance, dt, dt_limit, timestep, ...
      numsteps, model_name)
   %UPDATE_FORCE_ADVANCE_GUARD Track persistent maxsubstep force-advance
   % streaks across full steps.
   %
   %  streak_dt = icemodel.timestepping.update_force_advance_guard( ...
   %     streak_dt, forced_advance, dt, dt_limit, timestep, numsteps, ...
   %     model_name)
   %
   % Consecutive forced advances can span up to one full forcing step. Above
   % that, this function raises an error and stops the run.
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
