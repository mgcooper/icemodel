function metadata = readRcmArtifactMetadata(filename)
   %READRCMARTIFACTMETADATA Read saved RCM provenance without payload arrays.
   %
   %  metadata = ...readRcmArtifactMetadata(filename)
   %
   % Returns the scalar artifact_metadata struct saved beside a staged met or
   % Data payload. Legacy, missing, unreadable, and malformed metadata all
   % return an empty struct so callers can apply their own unknown-provenance
   % policy without loading the potentially large timetable.

   arguments
      filename (1, 1) string
   end

   metadata = struct();
   try
      info = whos('-file', filename);
   catch
      return
   end
   if ~any(string({info.name}) == "artifact_metadata")
      return
   end

   % Load only the small provenance variable; RCM timetable payloads are large.
   saved = load(filename, 'artifact_metadata', '-mat');
   if isstruct(saved.artifact_metadata) && isscalar(saved.artifact_metadata)
      metadata = saved.artifact_metadata;
   end
end
