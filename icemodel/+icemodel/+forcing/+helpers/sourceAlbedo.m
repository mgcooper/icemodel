function [albedo, qc_counts] = sourceAlbedo(swd, swu, kwargs)
   %SOURCEALBEDO Derive broadband albedo where shortwave input is valid.
   %
   % Ratios below 10 W m-2 downwelling are dominated by low-sun sensor noise.
   % Nonpositive reflected shortwave is also not a usable albedo observation.
   % When timestamp and location are supplied together, radiometer ratios below
   % 20 degrees solar elevation are rejected as low-angle measurements. Callers
   % may override that angle or impose a source-specific physical minimum.
   % Leave rejected samples missing so forcing builders can apply their explicit
   % gap-fill policy afterward; the input radiation is never modified.
   arguments
      swd
      swu
      kwargs.minimum (1, 1) double {mustBeNonnegative} = 0
      kwargs.Time datetime = datetime.empty(0, 1)
      kwargs.latitude (1, 1) double = NaN
      kwargs.longitude (1, 1) double = NaN
      kwargs.minimum_solar_elevation (1, 1) double {mustBeFinite} = 20
   end

   % Solar screening is optional, but partial geometry would silently produce a
   % different QC policy. Require the complete timestamp/location tuple.
   has_solar_geometry = [~isempty(kwargs.Time), ...
      isfinite(kwargs.latitude), isfinite(kwargs.longitude)];
   if any(has_solar_geometry) && ~all(has_solar_geometry)
      error('icemodel:forcing:helpers:sourceAlbedo:incompleteSolarGeometry', ...
         'Time, latitude, and longitude must be supplied together')
   end

   swdown_floor = 10;
   albedo = swu ./ swd;
   finite_ratio = isfinite(albedo);
   low_light = finite_ratio & swd < swdown_floor;
   low_solar_elevation = false(size(albedo));
   if all(has_solar_geometry)
      if numel(kwargs.Time) ~= numel(albedo)
         error('icemodel:forcing:helpers:sourceAlbedo:timeSizeMismatch', ...
            'Time must have one timestamp per shortwave sample')
      end
      solar_elevation = icemodel.forcing.helpers.solarElevation( ...
         kwargs.Time, kwargs.latitude, kwargs.longitude);
      solar_elevation = reshape(solar_elevation, size(albedo));
      low_solar_elevation = finite_ratio ...
         & solar_elevation <= kwargs.minimum_solar_elevation;
   end
   nonpositive_swu = finite_ratio & swd >= swdown_floor & swu <= 0;
   below_minimum = finite_ratio & swd >= swdown_floor & swu > 0 ...
      & albedo < kwargs.minimum;
   invalid = ~finite_ratio | low_light | low_solar_elevation ...
      | nonpositive_swu | below_minimum;
   albedo(invalid) = NaN;
   qc_counts = struct( ...
      'low_light', nnz(low_light), ...
      'low_solar_elevation', nnz(low_solar_elevation), ...
      'nonpositive_swu', nnz(nonpositive_swu), ...
      'below_minimum', nnz(below_minimum), ...
      'total', nnz(low_light | low_solar_elevation ...
      | nonpositive_swu | below_minimum));
end
