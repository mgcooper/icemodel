function writeManifest(pathname, manifest)
   %WRITEMANIFEST Write one verification family manifest JSON file.
   %
   %  icemodel.verification.setup.writeManifest(pathname, manifest)
   %
   % Inputs
   %  pathname   Destination manifest JSON path.
   %  manifest   Family manifest struct.
   %
   % Role
   %  Setup helper used by import/update tooling. Normal verification workflow
   %  functions read manifests but do not write them.

   % Write pretty JSON so committed manifest diffs stay reviewable. MATLAB
   % jsonencode serializes a scalar struct as an object, so convert top-level
   % case/skip struct arrays to cells at the file boundary to preserve JSON
   % array shape for one-case manifests.
   manifest_for_json = normalizeManifestJsonTypes(manifest);
   if isfield(manifest_for_json, 'cases')
      manifest_for_json.cases = num2cell(manifest_for_json.cases);
   end
   if isfield(manifest_for_json, 'skipped')
      manifest_for_json.skipped = num2cell(manifest_for_json.skipped);
   end

   % Repeated additive imports must be byte- and mtime-stable. Build the exact
   % payload first and leave an identical or field-order-only equivalent
   % existing manifest untouched.
   payload = jsonencode(manifest_for_json, PrettyPrint=true);
   if isfile(pathname)
      existing_payload = fileread(pathname);
      if strcmp(existing_payload, payload) ...
            || manifestPayloadsEquivalent(existing_payload, manifest_for_json)
         return
      end
   end

   fid = fopen(pathname, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, payload, 'char');
end

function tf = manifestPayloadsEquivalent(existing_payload, manifest_for_json)
   %MANIFESTPAYLOADSEQUIVALENT Compare JSON objects independent of field order.
   tf = false;
   % jsondecode maps both scalar null and [] to numeric [], so the semantic
   % comparison below cannot see this schema type drift. Force one canonical
   % rewrite when an existing flat site_location object contains [].
   if hasNullableSiteLocationArray(existing_payload)
      return
   end
   try
      existing = jsondecode(existing_payload);
      % Decode the already-valid incoming JSON as well so both comparison
      % operands have identical MATLAB container/type inference. Comparing a
      % decoded file directly with the caller's heterogeneous cell/struct graph
      % can otherwise report a false semantic difference.
      incoming = jsondecode(jsonencode(manifest_for_json));
   catch
      % A malformed prior file is not equivalent and must be replaced by the
      % valid manifest payload supplied by the caller.
      return
   end

   % Restore schema-defined singleton object arrays that jsondecode collapses,
   % then compare canonical object-key order. Cell/struct-array element order
   % is never sorted, so case ordering and every JSON value/type stay material.
   existing = normalizeManifestJsonTypes(existing);
   incoming = normalizeManifestJsonTypes(incoming);
   if isfield(existing, 'cases')
      existing.cases = num2cell(existing.cases);
   end
   if isfield(incoming, 'cases')
      incoming.cases = num2cell(incoming.cases);
   end
   if isfield(existing, 'skipped')
      existing.skipped = num2cell(existing.skipped);
   end
   if isfield(incoming, 'skipped')
      incoming.skipped = num2cell(incoming.skipped);
   end
   % A cache-reuse diagnostic is emitted at runtime and must not turn an
   % otherwise exact additive replay into a manifest mutation. Preserve real
   % provenance notes; remove only the two exact full-coverage messages owned
   % by stageRcmForcing before the semantic comparison.
   existing = stripTransientReuseNotes(existing);
   incoming = stripTransientReuseNotes(incoming);
   existing = canonicalizeObjectFields(existing);
   incoming = canonicalizeObjectFields(incoming);
   tf = strcmp(jsonencode(existing), jsonencode(incoming));
end

function tf = hasNullableSiteLocationArray(payload)
   %HASNULLABLESITELOCATIONARRAY Detect noncanonical arrays in scalar fields.
   fields = 'lat_wgs84|lon_wgs84|x_epsg3413|y_epsg3413|elev_m';
   pattern = ['"site_location"\s*:\s*\{[^{}]*"(?:' fields ...
      ')"\s*:\s*\[\]'];
   tf = ~isempty(regexp(payload, pattern, 'once'));
end

function value = stripTransientReuseNotes(value)
   %STRIPTRANSIENTREUSENOTES Ignore exact full-coverage replay diagnostics.
   if iscell(value)
      % Cells represent JSON arrays. Normalize each member without changing
      % element order or any non-diagnostic value.
      for k = 1:numel(value)
         value{k} = stripTransientReuseNotes(value{k});
      end
      return
   end
   if ~isstruct(value)
      return
   end

   % Manifests convert object arrays to cells before this comparison, so a
   % scalar nested source leg can safely drop its one transient note field.
   transient = [ ...
      "Existing staged met/Data files fully cover requested window."
      "Existing staged Data file fully covers requested window."];
   if isscalar(value) && isfield(value, 'note')
      note = string(value.note);
      if isscalar(note) && ismember(note, transient)
         value = rmfield(value, 'note');
      end
   end

   % Walk every remaining object field so nested colocation legs receive the
   % same narrow normalization while case order and scientific notes survive.
   fields = fieldnames(value);
   for k = 1:numel(value)
      for n = 1:numel(fields)
         name = fields{n};
         value(k).(name) = stripTransientReuseNotes(value(k).(name));
      end
   end
end

function value = canonicalizeObjectFields(value)
   %CANONICALIZEOBJECTFIELDS Sort object keys recursively, never array values.
   if iscell(value)
      % Cells represent JSON arrays, whose element order is semantically
      % significant even when their members are objects.
      for k = 1:numel(value)
         value{k} = canonicalizeObjectFields(value{k});
      end
      return
   end
   if ~isstruct(value)
      return
   end

   % Struct fields represent JSON object keys and may be sorted safely. Walk
   % each array element in place so arrays of objects retain their order.
   fields = sort(fieldnames(value));
   value = orderfields(value, fields);
   for k = 1:numel(value)
      for n = 1:numel(fields)
         name = fields{n};
         value(k).(name) = canonicalizeObjectFields(value(k).(name));
      end
   end
end

function value = normalizeManifestJsonTypes(value)
   %NORMALIZEMANIFESTJSONTYPES Preserve schema-defined JSON types on rewrites.

   if iscell(value)
      % Raw-decoded manifests can contain cells at several schema levels. Walk
      % every member without changing the surrounding cell-array contract.
      for k = 1:numel(value)
         value{k} = normalizeManifestJsonTypes(value{k});
      end
      return
   end
   if ~isstruct(value)
      return
   end

   % Restore the two schema-defined types that jsondecode cannot distinguish
   % from ordinary MATLAB values: nullable coordinate scalars and nested
   % forcing-window JSON arrays.
   fields = string(fieldnames(value));
   for k = 1:numel(value)
      for field = reshape(fields, 1, [])
         name = char(field);
         child = value(k).(name);
         if field == "site_location" && isstruct(child)
            % jsondecode maps scalar JSON null coordinates to numeric []. Put
            % those schema-defined scalars back to NaN so jsonencode preserves
            % null instead of changing untouched cases to JSON arrays.
            value(k).(name) = normalizeSiteLocationNulls(child);
         elseif field == "forcing_complete_windows"
            if isempty(child)
               value(k).(name) = cell(0, 1);
            elseif isstruct(child)
               value(k).(name) = num2cell(child(:));
            elseif iscell(child)
               value(k).(name) = reshape(child, [], 1);
            else
               error('icemodel:verification:writeManifest:badForcingWindows', ...
                  'forcing_complete_windows must be empty, struct, or cell')
            end
         else
            value(k).(name) = normalizeManifestJsonTypes(child);
         end
      end
   end
end

function location = normalizeSiteLocationNulls(location)
   %NORMALIZESITELOCATIONNULLS Preserve nullable scalar coordinate JSON types.
   nullable = ["lat_wgs84", "lon_wgs84", "x_epsg3413", ...
      "y_epsg3413", "elev_m"];
   present = intersect(nullable, string(fieldnames(location)), 'stable');
   for k = 1:numel(location)
      for field = reshape(present, 1, [])
         name = char(field);
         if isnumeric(location(k).(name)) && isempty(location(k).(name))
            location(k).(name) = NaN;
         end
      end
   end
end
