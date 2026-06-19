function [wspd, wdir] = windFromComponents(u, v)
   %WINDFROMCOMPONENTS Wind speed and direction from zonal/meridional components.
   %
   %  [wspd, wdir] = icemodel.forcing.helpers.windFromComponents(u, v)
   %
   % Inputs
   %  u - zonal (eastward) wind component [m s-1]
   %  v - meridional (northward) wind component [m s-1]
   %
   % Outputs
   %  wspd - wind speed [m s-1]
   %  wdir - meteorological wind direction [degrees]: the direction the
   %         wind blows FROM, measured clockwise from north, in (0, 360].
   %
   % See also: icemodel.forcing.helpers.metchecks

   wspd = hypot(u, v);
   wdir = atan2d(u, v) + 180;
end
