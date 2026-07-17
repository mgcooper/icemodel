function metadata = artifactMetadata(value)
   %ARTIFACTMETADATA Build a source-light top-level artifact metadata record.
   %
   %  metadata = icemodel.forcing.helpers.artifactMetadata(value)
   %
   % VALUE may be a metadata struct or a table/timetable. Table UserData is
   % preserved, and legacy Lat/Lon CustomProperties fill missing
   % lat_wgs84/lon_wgs84 identity fields. Writers save this record beside the
   % payload so reuse checks need not load a large timetable.

   metadata = struct();
   if isstruct(value)
      metadata = value;
      return
   end
   if ~istable(value) && ~istimetable(value)
      return
   end
   if isstruct(value.Properties.UserData)
      metadata = value.Properties.UserData;
   end

   % Derive the actual saved cadence rather than trusting a filename or caller
   % marker. This top-level copy lets later reuse/prune checks remain source-light.
   if istimetable(value) && height(value) >= 2
      steps = seconds(diff(value.Time));
      cadence = median(steps, 'omitnan');
      if isfinite(cadence) && cadence > 0 ...
            && all(isfinite(steps)) && all(abs(steps - cadence) < 1e-6)
         metadata.artifact_cadence_seconds = cadence;
      elseif isfield(metadata, 'artifact_cadence_seconds')
         metadata = rmfield(metadata, 'artifact_cadence_seconds');
      end
   end

   % Custom location properties predate direct point fields in UserData. Fill
   % only absent facts so explicit source metadata remains authoritative.
   custom = value.Properties.CustomProperties;
   names = string(fieldnames(custom));
   if ismember("Lat", names) && ~isfield(metadata, 'lat_wgs84')
      metadata.lat_wgs84 = custom.Lat;
   end
   if ismember("Lon", names) && ~isfield(metadata, 'lon_wgs84')
      metadata.lon_wgs84 = custom.Lon;
   end
end
