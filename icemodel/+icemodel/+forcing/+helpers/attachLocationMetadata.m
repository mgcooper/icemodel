function Data = attachLocationMetadata(Data, location)
   %ATTACHLOCATIONMETADATA Add location CustomProperties to a Data timetable.
   %
   %  Data = icemodel.forcing.helpers.attachLocationMetadata(Data, location)
   %
   % LOCATION uses the WGS84 or projected location fields that projectLocation
   % accepts, plus elev_m. An optional slope field gives the native surface
   % slope in m/m. A source with no slope metadata gets NaN.

   % Project only locations that do not already carry finite EPSG:3413 values.
   location = icemodel.forcing.helpers.projectLocation(location);

   % Keep the surface slope that the source gives. Use NaN when the source has
   % no slope, as the other forcing families do.
   slope = NaN;
   if isfield(location, 'slope')
      slope = location.slope;
   end

   % Attach the same CustomProperties schema that every Data builder uses.
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
