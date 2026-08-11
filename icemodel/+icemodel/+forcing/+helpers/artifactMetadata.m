function metadata = artifactMetadata(value)
   %ARTIFACTMETADATA Build a source-light top-level artifact metadata record.
   %
   %  metadata = icemodel.forcing.helpers.artifactMetadata(value)
   %
   % VALUE may be a metadata struct, a table, or a timetable. The function keeps
   % table UserData. It fills missing lat_wgs84/lon_wgs84 identity fields from
   % Lat/Lon CustomProperties. Writers save this record beside the payload, so
   % reuse checks do not load a large timetable.

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
   % marker. This top-level copy keeps later reuse and prune checks
   % source-light.
   if istimetable(value)
      cadence = icemodel.forcing.helpers.uniformCadenceSeconds(value);
      if isfinite(cadence)
         metadata.artifact_cadence_seconds = cadence;
      elseif isfield(metadata, 'artifact_cadence_seconds')
         metadata = rmfield(metadata, 'artifact_cadence_seconds');
      end
   end

   % Some tables carry the location in CustomProperties instead of direct point
   % fields in UserData. Fill only absent fields, so explicit source metadata
   % stays authoritative.
   custom = value.Properties.CustomProperties;
   names = string(fieldnames(custom));
   if ismember("Lat", names) && ~isfield(metadata, 'lat_wgs84')
      metadata.lat_wgs84 = custom.Lat;
   end
   if ismember("Lon", names) && ~isfield(metadata, 'lon_wgs84')
      metadata.lon_wgs84 = custom.Lon;
   end
end
