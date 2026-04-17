function Qsi = incoming_shortwave_radiation(J_day_start, step, dt, xlat, ...
      cloud_frac, slope_az, terrain_slope, ihrs_day, transmiss, ihour)
   %incoming_shortwave_radiation Estimate downwelling shortwave radiation.
   %
   %  Qsi = icemodel.surface.incoming_shortwave_radiation( ...
   %     J_day_start, step, dt, xlat, cloud_frac, slope_az, ...
   %     terrain_slope, ihrs_day, transmiss, ihour)
   %
   % Inputs:
   %  J_day_start   - day-of-year for the first forcing step [1-366]
   %  step          - 1-based forcing-step index from J_day_start [1]
   %  dt            - forcing timestep duration [s]
   %  xlat          - site latitude [degrees north]
   %  cloud_frac    - cloud fraction [1]
   %  slope_az      - terrain aspect, north-zero clockwise [degrees]
   %  terrain_slope - terrain slope angle [degrees]
   %  ihrs_day      - hourly samples per day for daily averaging [1]
   %  transmiss     - bulk atmospheric transmissivity [1]
   %  ihour         - hour of day for subdaily evaluation [0-24]
   %
   % Output:
   %  Qsi - incoming shortwave radiation [W m^-2]
   %
   % For daily forcing (`dt == 86400`), this helper returns the arithmetic
   % mean of hourly terrain-adjusted shortwave over one day. For subdaily
   % forcing, it evaluates the instantaneous terrain-adjusted flux at the
   % provided `ihour` and uses a continuous step-based Julian-day offset so
   % the seasonal declination term evolves smoothly for arbitrary `dt`.
   %
   % This is a forcing fallback helper, not the canonical spectral or SEB
   % shortwave machinery used once `swd` is already known.
   %
   %#codegen

   if dt == 86400.0
      % Daily forcing uses the mean of hourly clear-sky / cloud-adjusted
      % terrain radiation samples across the requested day.
      J_day = J_day_start + step - 1;
      Qsi_sum = 0.0;
      for hour_sample = 1:ihrs_day
         Qsi_hour = icemodel.surface.terrain_adjusted_shortwave_radiation( ...
            J_day, xlat, cloud_frac, hour_sample, slope_az, terrain_slope, ...
            transmiss);
         Qsi_sum = Qsi_sum + Qsi_hour;
      end
      Qsi = Qsi_sum / ihrs_day;
   else
      % Preserve the historical 1-based step interpretation while
      % generalizing the day-fraction update beyond the old hourly-only
      % `step / 24` expression.
      J_day = J_day_start - 1 + step * dt / 86400.0;
      Qsi = icemodel.surface.terrain_adjusted_shortwave_radiation( ...
         J_day, xlat, cloud_frac, ihour, slope_az, terrain_slope, ...
         transmiss);
   end
end
