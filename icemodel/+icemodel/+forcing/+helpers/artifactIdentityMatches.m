function tf = artifactIdentityMatches(filename, expected, variable_name)
   %ARTIFACTIDENTITYMATCHES Reject reuse across concrete provenance conflicts.
   %
   %  tf = icemodel.forcing.helpers.artifactIdentityMatches( ...
   %     filename, expected, variable_name)
   %
   % EXPECTED is a timetable/table or metadata struct. The candidate artifact is
   % read from FILENAME using its top-level artifact_metadata and, when supplied,
   % VARIABLE_NAME (for example "met" or "Data"). Missing legacy metadata is
   % treated as unknown and remains reusable; conflicting known source/product,
   % schema, sampling-method, or point metadata returns false. A malformed
   % one-coordinate candidate is rejected when EXPECTED has a concrete point.

   arguments
      filename (1, 1) string
      expected
      variable_name (1, 1) string = ""
   end

   % A missing or unreadable file cannot be a compatible reuse candidate.
   tf = false;
   if ~isfile(filename)
      return
   end
   try
      inventory = whos('-file', filename);
   catch
      return
   end
   names = string({inventory.name});
   saved = struct();
   if ismember("artifact_metadata", names)
      saved = load(filename, 'artifact_metadata');
   end

   % Prefer the durable top-level copy, then fill only absent identity values
   % from the saved table and its location CustomProperties.
   candidate = struct();
   if isfield(saved, 'artifact_metadata') ...
         && isstruct(saved.artifact_metadata)
      candidate = saved.artifact_metadata;
   end
   expected = icemodel.forcing.helpers.artifactMetadata(expected);
   [~, candidate_point_state] = identityPoint(candidate);
   [~, expected_point_state] = identityPoint(expected);
   needs_legacy_point = expected_point_state == "known" ...
      && candidate_point_state == "unknown";
   if variable_name ~= "" ...
         && (isempty(fieldnames(candidate)) || needs_legacy_point) ...
         && ismember(variable_name, names)
      % Legacy artifacts without the top-level copy require one timetable read;
      % current artifacts remain source-light even when their payload is large.
      try
         saved_table = load(filename, char(variable_name));
      catch
         saved_table = struct();
      end
      if isfield(saved_table, char(variable_name))
         candidate = fillIdentity(candidate, ...
            icemodel.forcing.helpers.artifactMetadata( ...
            saved_table.(char(variable_name))));
      end
   end

   % Share scalar and documented-alias comparison with manifest merging while
   % preserving this file-reuse boundary's directional malformed-point rule.
   if ~icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
         candidate, expected)
      return
   end

   % Point mismatches are concrete conflicts. Unknown legacy points remain
   % eligible, while half-populated candidate coordinates are malformed.
   [candidate_point, candidate_state] = identityPoint(candidate);
   [expected_point, expected_state] = identityPoint(expected);
   if expected_state == "known"
      if candidate_state == "partial"
         return
      end
      if candidate_state == "known" ...
            && any(abs(candidate_point - expected_point) > 1e-8)
         return
      end
   end
   tf = true;
end

function target = fillIdentity(target, source)
   %FILLIDENTITY Add only absent candidate facts from a saved table.
   fields = fieldnames(source);
   for k = 1:numel(fields)
      name = fields{k};
      if ~isfield(target, name)
         target.(name) = source.(name);
      end
   end
end

function [point, state] = identityPoint(metadata)
   %IDENTITYPOINT Return a known, partial, or unknown [lat lon] identity.
   % Current artifacts use direct coordinate fields; some older producers use
   % a point vector or nested location record.
   [point, state] = coordinateFields(metadata);
   if state ~= "unknown"
      return
   end
   if isfield(metadata, 'point') && isnumeric(metadata.point) ...
         && numel(metadata.point) == 2
      point = reshape(double(metadata.point), 1, 2);
      state = pointState(point);
      return
   end
   for name = ["location", "site_location"]
      if isfield(metadata, char(name)) && isstruct(metadata.(char(name)))
         [point, state] = coordinateFields(metadata.(char(name)));
         if state ~= "unknown"
            return
         end
      end
   end
end

function [point, state] = coordinateFields(metadata)
   %COORDINATEFIELDS Read a direct lat_wgs84/lon_wgs84 pair.
   point = [NaN, NaN];
   has_lat = isfield(metadata, 'lat_wgs84');
   has_lon = isfield(metadata, 'lon_wgs84');
   if ~has_lat && ~has_lon
      state = "unknown";
      return
   end
   if has_lat && has_lon && isscalar(metadata.lat_wgs84) ...
         && isscalar(metadata.lon_wgs84)
      point = [double(metadata.lat_wgs84), double(metadata.lon_wgs84)];
      state = pointState(point);
   else
      state = "partial";
   end
end

function state = pointState(point)
   %POINTSTATE Classify one two-coordinate candidate.
   if all(isfinite(point))
      state = "known";
   elseif any(isfinite(point))
      state = "partial";
   else
      state = "unknown";
   end
end
