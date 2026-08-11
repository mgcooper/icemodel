function tf = artifactCadenceMatches(filename, variable_name, expected_seconds)
   %ARTIFACTCADENCEMATCHES Check a saved timetable for the requested cadence.
   %
   %  tf = icemodel.forcing.helpers.artifactCadenceMatches( ...
   %     filename, variable_name, expected_seconds)
   %
   % Writers save artifact_cadence_seconds beside the payload, and this
   % function reads that field when the file has it. An older file without the
   % field is reusable only when the saved timetable has at least two rows and
   % a uniform time axis that matches expected_seconds.

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

   % Use the top-level cadence the writer stored when the metadata has it.
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

   % A file without that metadata needs the actual table axis. The variable
   % name alone does not show the cadence.
   if ~ismember(variable_name, names)
      return
   end
   saved = load(filename, char(variable_name));
   value = saved.(char(variable_name));
   if ~istimetable(value)
      return
   end
   candidate = icemodel.forcing.helpers.uniformCadenceSeconds(value);
   tf = isfinite(candidate) && abs(candidate - expected_seconds) < 1e-6;
end
