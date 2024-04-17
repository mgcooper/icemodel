function Qsi = SOLARRAD(Jday,xlat,cloudfrac,xhour,slopeaz,sfcslope,transmiss)
   %SOLAR_RAD Compute the incoming solar radiation at the model timestep

   % Required constants.
   solar_const = 1370.;
   days_yr = 365.25;
   Trop_Can = 0.41;
   solstice = 173.;
   pi = 2.0 * acos(0.0);
   deg2rad = pi / 180.0;

   % COMPUTE THE BASIC SOLAR RADIATION PARAMETERS.

   % Compute the solar declination angle (radians).
   sol_dec = Trop_Can * cos(2.*pi * (real(Jday) - solstice)/days_yr);

   % Compute the sun's hour angle (radians).
   hr_angl = (xhour * 15.0 - 180.0) * deg2rad;

   % Compute cos_Z.  Note that the sin of the solar elevation angle,
   %   sin_alfa, is equal to the cosine of the solar zenith angle, cos_Z.
   cos_Z = sin(sol_dec) * sin(xlat * deg2rad) + ...
      cos(sol_dec) * cos(xlat * deg2rad) * cos(hr_angl);
   cos_Z = max(0.0,cos_Z);

   % Account for clouds, water vapor, pollution, etc.
   trans_direct = transmiss * (0.6 + 0.2 * cos_Z) * (1.0-cloudfrac);
   trans_diffuse = transmiss * (0.3 + 0.1 * cos_Z) * cloudfrac;

   % Compute the solar radiation transmitted through the atmosphere.
   Qsi_trans_dir = solar_const * trans_direct;
   Qsi_trans_dif = solar_const * trans_diffuse;

   % COMPUTE THE CORRECTIONS TO ALLOW FOR TOPOGRAPHIC SLOPE AND ASPECT.

   % The sine of the solar zenith angle.
   sin_Z = sqrt(1.0 - cos_Z*cos_Z);

   % Azimuth of the sun, with south having zero azimuth.
   sun_azimuth = asin(max(-1.0,min(1.0,cos(sol_dec)*sin(hr_angl)/sin_Z)));

   % Make the corrections so that the angles below the local horizon
   %   are still measured from the normal to the slope.
   if (hr_angl<0.0)
      if (hr_angl<sun_azimuth)
         sun_azimuth = - pi - sun_azimuth;
      elseif (hr_angl>0.0)
         if (hr_angl>sun_azimuth)
            sun_azimuth = pi - sun_azimuth;
         end
      end
   end

   % Build, from the variable with north having zero azimuth, a
   %   slope_azimuth value with south having zero azimuth.
   if (slopeaz>=180.0)
      slope_az_S0 = slopeaz - 180.0;
   else
      slope_az_S0 = slopeaz + 180.0;
   end
   % Compute the angle between the normal to the slope and the angle
   %   at which the direct solar radiation impinges on the sloping
   %   terrain (radians).
   cos_i = cos(sfcslope * deg2rad) * cos_Z + ...
      sin(sfcslope * deg2rad) * sin_Z * ...
      cos(sun_azimuth - slope_az_S0 * deg2rad);
   % Adjust the topographic correction due to local slope so that
   %   the correction is zero if the sun is below the local horizon
   %   (i.e., the slope is in the shade) or if the sun is below the
   %   global horizon.
   if (cos_i<0.0)
      cos_i = 0.0;
   end
   if (cos_Z<=0.0)
      cos_i = 0.0;
   end
   % Adjust the solar radiation for slope, etc.
   Qsi_direct = cos_i * Qsi_trans_dir;
   Qsi_diffuse = cos_Z * Qsi_trans_dif;

   % Combine the direct and diffuse solar components.
   Qsi = Qsi_direct + Qsi_diffuse;
end
