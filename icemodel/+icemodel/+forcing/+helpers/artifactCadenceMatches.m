function tf = artifactCadenceMatches(filename, variable_name, expected_seconds)
   %ARTIFACTCADENCEMATCHES Prove a saved timetable has the requested cadence.
   %
   %  tf = icemodel.forcing.helpers.artifactCadenceMatches( ...
   %     filename, variable_name, expected_seconds)
   %
   % Current writers save artifact_cadence_seconds beside the payload, keeping
   % checks source-light. Legacy artifacts without that field remain reusable only
   % when the saved timetable has at least two rows and a uniform matching axis.

   arguments
      filename (1, 1) string
      variable_name (1, 1) string
      expected_seconds (1, 1) double
   end

   tf = false;
   if ~isfile(filename) || ~isfinite(expected_seconds) ...
         || expected_seconds <= 0
      return
   end
   try
      inventory = whos('-file', filename);
   catch
      return
   end
   names = string({inventory.name});

   % Trust the writer-derived top-level cadence when current metadata provides it.
   if ismember("artifact_metadata", names)
      saved = load(filename, 'artifact_metadata');
      if isstruct(saved.artifact_metadata) ...
            && isfield(saved.artifact_metadata, 'artifact_cadence_seconds')
         candidate = double(saved.artifact_metadata.artifact_cadence_seconds);
         tf = isscalar(candidate) && isfinite(candidate) ...
            && abs(candidate - expected_seconds) < 1e-6;
         return
      end
   end

   % Legacy files require the actual table axis; a label alone is not proof.
   if ~ismember(variable_name, names)
      return
   end
   saved = load(filename, char(variable_name));
   value = saved.(char(variable_name));
   if ~istimetable(value) || height(value) < 2
      return
   end
   steps = seconds(diff(value.Time));
   candidate = median(steps, 'omitnan');
   tf = isfinite(candidate) && candidate > 0 ...
      && all(isfinite(steps)) && all(abs(steps - candidate) < 1e-6) ...
      && abs(candidate - expected_seconds) < 1e-6;
end
