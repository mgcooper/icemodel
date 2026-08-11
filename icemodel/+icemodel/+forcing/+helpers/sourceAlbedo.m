function [albedo, qc_counts] = sourceAlbedo(swd, swu, kwargs)
   %SOURCEALBEDO Derive broadband albedo where shortwave input is valid.
   %
   % Ratios below 10 W m-2 downwelling are dominated by low-sun sensor noise.
   % Nonpositive reflected shortwave is also not a usable albedo observation.
   % If you supply timestamp and location together, the function rejects
   % radiometer ratios below 20 degrees solar elevation as low-angle
   % measurements. Callers can override that angle or set a source-specific
   % physical minimum. The function leaves rejected samples missing, so forcing
   % builders can then apply their own gap-fill policy. The function does not
   % change the input radiation.
   arguments
      swd
      swu
      kwargs.minimum (1, 1) double {mustBeNonnegative} = 0
      kwargs.Time datetime = datetime.empty(0, 1)
      kwargs.latitude (1, 1) double = NaN
      kwargs.longitude (1, 1) double = NaN
      kwargs.minimum_solar_elevation (1, 1) double {mustBeFinite} = 20
   end

   % Solar screening is optional. Partial geometry would apply a different QC
   % policy, so require the complete timestamp, latitude, and longitude set.
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
