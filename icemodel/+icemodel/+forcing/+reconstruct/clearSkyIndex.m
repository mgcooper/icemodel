function [csi, toa] = clearSkyIndex(times, values, location, kwargs)
   %CLEARSKYINDEX Normalize shortwave flux by station-specific TOA irradiance.
   %
   %  [csi, toa] = icemodel.forcing.reconstruct.clearSkyIndex( ...
   %     times, values, location)
   %
   % CSI is undefined below the policy's meaningful-sun threshold. Keeping
   % those samples missing prevents nighttime geometry from entering donor
   % fits, lag searches, or daylight interpolation.

   arguments
      times datetime
      values (:, 1) double
      location (1, 1) struct
      kwargs.toa_dark_wm2 (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().toa_dark_wm2
   end

   if numel(values) ~= numel(times)
      error('icemodel:reconstruct:clearSkyIndex:sizeMismatch', ...
         'values and times must share one sample axis');
   end
   required = ["lat_wgs84", "lon_wgs84"];
   if ~all(isfield(location, required)) ...
         || ~isscalar(location.lat_wgs84) || ~isscalar(location.lon_wgs84) ...
         || ~isfinite(location.lat_wgs84) || ~isfinite(location.lon_wgs84)
      error('icemodel:reconstruct:clearSkyIndex:missingSolarGeometry', ...
         'CSI requires finite scalar latitude and longitude');
   end

   toa = icemodel.forcing.reconstruct.toaIrradiance(times, ...
      location.lat_wgs84, location.lon_wgs84);
   csi = values ./ toa;
   csi(toa < kwargs.toa_dark_wm2) = NaN;
end
