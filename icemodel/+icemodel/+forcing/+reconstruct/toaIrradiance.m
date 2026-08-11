function toa = toaIrradiance(times, latitude, longitude)
   %TOAIRRADIANCE Top-of-atmosphere irradiance on a horizontal surface.
   %
   %  toa = icemodel.forcing.reconstruct.toaIrradiance(times, 67.0, -48.8)
   %
   % Role
   %  This function computes the clear-sky reference for the reconstruct
   %  namespace. Every consumer of "meaningful sun" as an irradiance calls
   %  it: the tier-1 CSI interpolation in fillShortGaps, the census
   %  daylight cut, and the synthetic-draw placement. The swd darkness
   %  zero-fill (reconstructSeries) keys on the civil-twilight solar
   %  elevation instead (D-28). It takes that elevation from the same
   %  icemodel.forcing.helpers.solarElevation geometry that this function
   %  scales, so both use the same solar math. This function scales the
   %  solar constant by the Sun elevation. The annual eccentricity
   %  correction (~3%) is second-order for CSI RATIOS because it cancels
   %  between numerator samples hours apart, so the plain form stays.
   %
   % See also: icemodel.forcing.reconstruct.fillShortGaps,
   %  icemodel.forcing.reconstruct.reconstructSeries

   elevation = icemodel.forcing.helpers.solarElevation( ...
      times, latitude, longitude);
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   toa = bands.solar_constant_wm2 * max(0, sind(elevation(:)));
end
