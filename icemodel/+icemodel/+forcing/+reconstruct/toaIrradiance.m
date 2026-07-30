function toa = toaIrradiance(times, latitude, longitude)
   %TOAIRRADIANCE Top-of-atmosphere irradiance on a horizontal surface.
   %
   %  toa = icemodel.forcing.reconstruct.toaIrradiance(times, 67.0, -48.8)
   %
   % Role
   %  Single source of the reconstruct namespace's clear-sky reference:
   %  every consumer of "meaningful sun" as an irradiance (the tier-1 CSI
   %  interpolation in fillShortGaps, the census daylight cut, the
   %  synthetic-draw placement) calls this one function. The swd darkness
   %  zero-fill (reconstructSeries) keys on the civil-twilight solar
   %  elevation instead (D-28) but derives it from the same
   %  icemodel.forcing.helpers.solarElevation geometry this function
   %  scales, so the solar math stays single-sourced. Solar-constant
   %  scaling by Sun elevation; the annual eccentricity correction (~3%)
   %  is second-order for CSI RATIOS because it cancels between numerator
   %  samples hours apart, so the plain form stays.
   %
   % See also: icemodel.forcing.reconstruct.fillShortGaps,
   %  icemodel.forcing.reconstruct.reconstructSeries

   elevation = icemodel.forcing.helpers.solarElevation( ...
      times, latitude, longitude);
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   toa = bands.solar_constant_wm2 * max(0, sind(elevation(:)));
end
