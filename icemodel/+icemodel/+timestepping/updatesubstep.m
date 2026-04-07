function [Ts, T, f_ice, f_liq, dt_sum, dt_new, liqflag, roL, ro_sno, cp_sno] ...
      = updatesubstep(Ts, T, f_ice, f_liq, dt_FULL_STEP, dt_sum, dt_new, TINY, ...
      ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs)
   % Checkpoint the accepted state and allocate time within the full step.
   %
   %#codegen

   persistent f_liq_phase_switch_threshold
   if isempty(f_liq_phase_switch_threshold)
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   % Allocate this substep to the timestep
   dt_sum = dt_sum + dt_new;

   % Adjust dt to exactly complete the full step without going over.
   % The first condition is true if the full step is incomplete, the second
   % is true if the next substep will exceed the full step.
   if (dt_FULL_STEP - dt_sum) > TINY && (dt_sum + dt_new - dt_FULL_STEP) > TINY
      dt_new = dt_FULL_STEP - dt_sum; % dt_new = max(dt - dt_sum, dt_min);
   end

   if nargout > 5

      % Top node contains enough liquid water to use the liquid-phase curve.
      liqflag = f_liq(1) > f_liq_phase_switch_threshold;

      % roL is for surface equations - subsurface node "wetness" varies
      if liqflag
         roL = roLv;  % ro_air * Lv
      else
         roL = roLs;  % ro_air * Ls
      end
   end

   if nargout > 7

      % UPDATE DENSITY, HEAT CAPACITY, DIFFUSION LENGTH SCALE
      ro_sno = f_ice * ro_ice + f_liq * ro_liq + (1.0 - f_liq - f_ice) * ro_air;
      cp_sno = (cv_ice * f_ice + cv_liq * f_liq) ./ ro_sno;

      % zD = sqrt(k_eff(1) * dt / (ro_sno(1) * cp_sno(1)));
      % if zD > dz(1)
      %  % placeholder
      % end
   end

   % an older option used to interpolate between timesteps
   % itime = itime + seconds(dt_new);
   % this was initialized in initialize_timesteps
   % itime = Time(1);
end
