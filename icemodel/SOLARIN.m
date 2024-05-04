function Qsi = SOLARIN(J_day_start, iter, dt, xlat, cloud_frac, ...
      slope_az, terrain_slope, ihrs_day, transmiss, ihour)
   %SOLARIN Compute the incoming solar radiation.
   %
   % This assumes a daily time step.
   %
   %#codegen

   J_day = iter + J_day_start - 1;
   Qsi_sum = 0.0;
   if (dt == 86400.0)
      for ihour = 1:ihrs_day
         Qsi = SOLARRAD(J_day, xlat, cloud_frac, ihour, slope_az, ...
            terrain_slope, transmiss);
         Qsi_sum = Qsi_sum + Qsi_tmp;
      end
      Qsi = Qsi_sum / ihrs_day;
   else
      J_day = iter / 24 + J_day_start - 1;
      Qsi = SOLARRAD(J_day, xlat, cloud_frac, ihour, slope_az, ...
         terrain_slope, transmiss);
   end
end
