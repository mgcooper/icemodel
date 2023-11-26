function [T, Tsfc, f_ice, f_liq, dt_sum, dt_new, dt_flag, liqflag, roL, ...
      ro_sno, cp_sno] = UPDATESUBSTEP(T, Tsfc, f_ice, f_liq, ro_ice, ro_liq, ...
      ro_air, cv_ice, cv_liq, dt, dt_sum, dt_new, roLv, roLs, dt_min, TINY)
   %UPDATESUBSTEP Return past values, allocate timestep, and update state

   % Allocate this substep to the timestep
   dt_sum = dt_sum + dt_new;

   % First is true if step incomplete, second if overage will occur
   dt_flag = (dt - dt_sum) > TINY && (dt_sum + dt_new - dt) > TINY;

   if dt_flag
      dt_new = max(dt - dt_sum, dt_min);
   end

   if nargout > 7
      
      % Top node contains >2% liquid water
      liqflag = f_liq(1) > 0.02;
      if liqflag
         roL = roLv;  % ro_air * Lv
      else
         roL = roLs;  % ro_air * Ls
      end
   end
      
   if nargout > 9

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
   % this was initialized in INITTIMESTEPS
   % itime = Time(1);
end
