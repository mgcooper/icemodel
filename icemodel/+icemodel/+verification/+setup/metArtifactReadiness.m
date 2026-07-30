function [tf, reason, complete_windows] = metArtifactReadiness(met_file)
   %METARTIFACTREADINESS Diagnose the exact saved scalar-window met artifact.
   %
   %  [tf, reason, complete_windows] = ...
   %     icemodel.verification.setup.metArtifactReadiness(met_file)
   %
   % Importers call this helper with the path returned by writemet, which may be
   % a newly written file, an exact no-overwrite reuse, or a broader enclosing
   % reuse. Readiness must describe those referenced bytes rather than the
   % request timetable that happened to precede the writer call.
   %
   % COMPLETE_WINDOWS is a JSON-portable struct column with UTC ISO-8601
   % start_time/end_time strings and numeric sample_count. The caller's existing
   % scalar met_files field is the artifact link; no absolute path is duplicated
   % in each window record.
   %
   % See also: icemodel.forcing.helpers.writemet,
   %  icemodel.verification.setup.metForcingReady

   arguments
      met_file string
   end

   % The current importer contract is intentionally one naming="window" file.
   % Reject vectors explicitly instead of inventing unused cross-file semantics.
   if ~isscalar(met_file)
      error('icemodel:verification:metArtifactReadiness:scalarWindowRequired', ...
         'met_file must name one scalar-window artifact')
   end

   % The three supported importers use naming="window", so an empty or missing
   % scalar path is a broken staging result rather than an unready artifact.
   if strlength(met_file) == 0 || ~isfile(met_file)
      error('icemodel:verification:metArtifactReadiness:fileMissing', ...
         'referenced saved met artifact does not exist: %s', met_file)
   end

   % Load only the canonical payload. Unrelated MAT variables stay out of memory
   % and cannot be mistaken for the runtime artifact under inspection.
   try
      saved = load(met_file, 'met');
   catch err
      error('icemodel:verification:metArtifactReadiness:loadFailed', ...
         'could not load referenced saved met artifact %s: %s', ...
         met_file, err.message)
   end
   if ~isfield(saved, 'met') || ~istimetable(saved.met)
      error('icemodel:verification:metArtifactReadiness:badPayload', ...
         'referenced saved met artifact %s lacks a timetable named met', ...
         met_file)
   end

   % A file selected by writemet must still satisfy the structural model-met
   % contract. Fail clearly on a malformed legacy/colliding file rather than
   % publishing advisory readiness for bytes runtime cannot load safely.
   try
      icemodel.forcing.helpers.validatemet(saved.met)
      [tf, reason, windows] = ...
         icemodel.verification.setup.metForcingReady(saved.met);
   catch err
      error('icemodel:verification:metArtifactReadiness:badPayload', ...
         'referenced saved met artifact %s is invalid: %s', ...
         met_file, err.message)
   end

   % Convert MATLAB table/datetime values at this boundary so nested manifest
   % records have one explicit, portable schema across all three importers.
   template = struct('start_time', "", 'end_time', "", 'sample_count', 0);
   complete_windows = repmat(template, height(windows), 1);
   for k = 1:height(windows)
      complete_windows(k).start_time = isoUtc(windows.start_time(k));
      complete_windows(k).end_time = isoUtc(windows.end_time(k));
      complete_windows(k).sample_count = double(windows.sample_count(k));
   end
end

function text = isoUtc(value)
   %ISOUTC Format one finite manifest diagnostic timestamp explicitly in UTC.

   value = icemodel.verification.setup.ensureUtc(value);
   % ensureUtc attaches UTC to naive values but intentionally preserves an
   % existing zone. Convert that zoned instant before appending the literal Z.
   value.TimeZone = 'UTC';
   text = string(value, "yyyy-MM-dd'T'HH:mm:ss'Z'");
end
