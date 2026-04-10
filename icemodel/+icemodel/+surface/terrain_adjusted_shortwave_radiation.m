function Qsi = terrain_adjusted_shortwave_radiation(Jday, xlat, ...
      cloudfrac, xhour, slopeaz, sfcslope, transmiss)
   %terrain_adjusted_shortwave_radiation Estimate terrain-adjusted shortwave.
   %
   %  Qsi = icemodel.surface.terrain_adjusted_shortwave_radiation( ...
   %     Jday, xlat, cloudfrac, xhour, slopeaz, sfcslope, transmiss)
   %
   % Inputs:
   %  Jday      - Julian day / day-of-year index [1-366]
   %  xlat      - latitude [degrees north]
   %  cloudfrac - cloud fraction [1]
   %  xhour     - hour of day [0-24]
   %  slopeaz   - terrain aspect, north-zero clockwise [degrees]
   %  sfcslope  - terrain slope angle [degrees]
   %  transmiss - bulk atmospheric transmissivity [1]
   %
   % Output:
   %  Qsi - terrain-adjusted incoming shortwave radiation [W m^-2]
   %
   % This is a legacy forcing fallback used when shortwave forcing must be
   % estimated from site geometry and cloud information. It applies the
   % original direct-plus-diffuse transmissivity parameterization and the
   % topographic slope/aspect correction on the instantaneous solar geometry.
   %
   %#codegen

   % Required constants from the translated legacy parameterization.
   solar_const = 1370.0;
   days_yr = 365.25;
   trop_can = 0.41;
   solstice = 173.0;
   pi_value = 2.0 * acos(0.0);
   deg2rad = pi_value / 180.0;

   % Compute the solar declination and hour-angle geometry in radians.
   sol_dec = trop_can * cos(2.0 * pi_value * (real(Jday) - solstice) / days_yr);
   hr_angl = (xhour * 15.0 - 180.0) * deg2rad;

   % The sine of the solar elevation angle equals cos_Z, the cosine of the
   % solar zenith angle. Clamp below the horizon to zero.
   cos_Z = sin(sol_dec) * sin(xlat * deg2rad) + ...
      cos(sol_dec) * cos(xlat * deg2rad) * cos(hr_angl);
   cos_Z = max(0.0, cos_Z);

   % Apply the legacy cloud / transmissivity parameterization to the direct
   % and diffuse atmospheric transmission terms.
   trans_direct = transmiss * (0.6 + 0.2 * cos_Z) * (1.0 - cloudfrac);
   trans_diffuse = transmiss * (0.3 + 0.1 * cos_Z) * cloudfrac;
   Qsi_trans_dir = solar_const * trans_direct;
   Qsi_trans_dif = solar_const * trans_diffuse;

   % Reconstruct the sun azimuth using the legacy south-zero convention.
   sin_Z = sqrt(max(0.0, 1.0 - cos_Z * cos_Z));
   if sin_Z <= 0.0 || cos_Z <= 0.0
      sun_azimuth = 0.0;
   else
      sin_sun_azimuth = cos(sol_dec) * sin(hr_angl) / sin_Z;
      sun_azimuth = asin(max(-1.0, min(1.0, sin_sun_azimuth)));

      % Keep the azimuth referenced from the slope normal even when the sun
      % moves below the local slope horizon.
      if hr_angl < 0.0
         if hr_angl < sun_azimuth
            sun_azimuth = -pi_value - sun_azimuth;
         end
      elseif hr_angl > 0.0
         if hr_angl > sun_azimuth
            sun_azimuth = pi_value - sun_azimuth;
         end
      end
   end

   % Convert the input aspect from north-zero clockwise to the south-zero
   % convention used by the translated slope-correction formula.
   if slopeaz >= 180.0
      slope_az_S0 = slopeaz - 180.0;
   else
      slope_az_S0 = slopeaz + 180.0;
   end

   % Project the direct beam onto the local slope normal.
   cos_i = cos(sfcslope * deg2rad) * cos_Z + ...
      sin(sfcslope * deg2rad) * sin_Z * ...
      cos(sun_azimuth - slope_az_S0 * deg2rad);
   cos_i = max(0.0, cos_i);
   if cos_Z <= 0.0
      cos_i = 0.0;
   end

   % Combine direct and diffuse terrain-adjusted components.
   Qsi_direct = cos_i * Qsi_trans_dir;
   Qsi_diffuse = cos_Z * Qsi_trans_dif;
   Qsi = Qsi_direct + Qsi_diffuse;
end
