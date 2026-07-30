function Data = attachLocationMetadata(Data, location)
   %ATTACHLOCATIONMETADATA Add location CustomProperties to a Data timetable.
   %
   %  Data = icemodel.forcing.helpers.attachLocationMetadata(Data, location)
   %
   % LOCATION uses the canonical WGS84/projected location fields accepted by
   % projectLocation plus elev_m. An optional slope field supplies the native
   % surface slope in m/m; sources without slope metadata receive NaN.

   % Project only locations that do not already carry finite EPSG:3413 values.
   location = icemodel.forcing.helpers.projectLocation(location);

   % Preserve a source-provided surface slope while keeping the common
   % no-source-metadata representation used by the other forcing families.
   slope = NaN;
   if isfield(location, 'slope')
      slope = location.slope;
   end

   % Attach one canonical CustomProperties schema across every Data builder.
   Data = addprop(Data, ...
      {'X', 'Y', 'Lat', 'Lon', 'Elev', 'Slope', 'ScalarUnits'}, ...
      repmat({'table'}, 1, 7));
   Data.Properties.CustomProperties.X = location.x_epsg3413;
   Data.Properties.CustomProperties.Y = location.y_epsg3413;
   Data.Properties.CustomProperties.Lat = location.lat_wgs84;
   Data.Properties.CustomProperties.Lon = location.lon_wgs84;
   Data.Properties.CustomProperties.Elev = location.elev_m;
   Data.Properties.CustomProperties.Slope = slope;
   Data.Properties.CustomProperties.ScalarUnits = ...
      ["m", "m", "degrees N", "degrees W", "m asl", "m/m"];
end
