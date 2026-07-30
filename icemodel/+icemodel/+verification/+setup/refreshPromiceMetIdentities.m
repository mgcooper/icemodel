function manifest = refreshPromiceMetIdentities(kwargs)
   %REFRESHPROMICEMETIDENTITIES Pin existing native met bytes in the manifest.
   %
   %  manifest = ...
   %     icemodel.verification.setup.refreshPromiceMetIdentities()
   %  manifest = ...
   %     icemodel.verification.setup.refreshPromiceMetIdentities( ...
   %     manifest_file="/path/to/data/eval/promice/manifest.json", ...
   %     met_root="/path/to/data/input/met")
   %
   % This is the POLICY A1/D-22 metadata-only refresh. It hashes only the
   % native met files already declared by the PROMICE producer manifest and
   % writes their relative path, byte size, and SHA-256 back to that manifest.
   % It never reads raw PROMICE sources or rebuilds staged artifacts.

   arguments
      kwargs.manifest_file (1, 1) string = ""
      kwargs.met_root (1, 1) string = ""
   end

   % Resolve the canonical manifest first; a custom manifest implies its paired
   % data/input/met root unless the caller supplies one explicitly.
   manifest_file = kwargs.manifest_file;
   if manifest_file == ""
      manifest_file = string(icemodel.internal.fullpath( ...
         'data', 'eval', 'promice', 'manifest.json'));
   end
   if ~isfile(manifest_file)
      error('icemodel:verification:refreshPromiceMetIdentities:missingManifest', ...
         'PROMICE producer manifest does not exist: %s', manifest_file)
   end
   met_root = kwargs.met_root;
   if met_root == ""
      data_root = fileparts(fileparts(fileparts(manifest_file)));
      met_root = string(fullfile(data_root, 'input', 'met'));
   end
   if ~isfolder(met_root)
      error('icemodel:verification:refreshPromiceMetIdentities:missingMetRoot', ...
         'Staged met root does not exist: %s', met_root)
   end

   % Decode and validate the one structural field this targeted repair owns.
   manifest = jsondecode(fileread(manifest_file));
   if ~isstruct(manifest) || ~isfield(manifest, 'cases') ...
         || ~isstruct(manifest.cases) || isempty(manifest.cases)
      error('icemodel:verification:refreshPromiceMetIdentities:invalidManifest', ...
         'PROMICE producer manifest has no case records: %s', manifest_file)
   end

   % Build every identity in memory before writing, so one missing, duplicate,
   % or escaping path leaves the existing manifest untouched.
   cases = manifest.cases;
   for k = 1:numel(cases)
      if ~isfield(cases(k), 'colocation') ...
            || ~isstruct(cases(k).colocation) ...
            || ~isfield(cases(k).colocation, 'promice') ...
            || ~isstruct(cases(k).colocation.promice)
         continue
      end
      leg = cases(k).colocation.promice;
      if ~isfield(leg, 'met_files') || isempty(leg.met_files)
         leg.met_file_identities = emptyIdentities();
         cases(k).colocation.promice = leg;
         continue
      end

      declared = reshape(string(leg.met_files), [], 1);
      declared = replace(declared, "\", "/");
      if any(strlength(declared) == 0) ...
            || numel(unique(declared)) ~= numel(declared)
         error(['icemodel:verification:refreshPromiceMetIdentities:' ...
            'ambiguousMetFiles'], ...
            'PROMICE case %s has empty or duplicate staged met paths.', ...
            caseId(cases(k)))
      end

      identities = repmat(struct('file', "", 'size_bytes', NaN, ...
         'sha256', ""), numel(declared), 1);
      for n = 1:numel(declared)
         filename = string(fullfile(met_root, declared(n)));
         if ~icemodel.internal.isPathInside(filename, met_root)
            error(['icemodel:verification:refreshPromiceMetIdentities:' ...
               'artifactOutsideRoot'], ...
               'PROMICE staged met path escapes %s: %s', met_root, declared(n))
         end
         if ~isfile(filename)
            error(['icemodel:verification:refreshPromiceMetIdentities:' ...
               'missingMetFile'], ...
               'PROMICE staged met file does not exist: %s', filename)
         end
         info = dir(filename);
         identities(n).file = declared(n);
         identities(n).size_bytes = info.bytes;
         identities(n).sha256 = ...
            icemodel.verification.setup.fileSha256(filename);
      end
      leg.met_file_identities = identities;
      cases(k).colocation.promice = leg;
   end

   % Use the canonical manifest writer so equivalent repeats remain byte-stable.
   manifest.cases = cases;
   icemodel.verification.setup.writeManifest(manifest_file, manifest);
end

function identities = emptyIdentities()
   %EMPTYIDENTITIES Return the schema-stable empty identity array.
   identities = struct('file', {}, 'size_bytes', {}, 'sha256', {});
end

function id = caseId(entry)
   %CASEID Return a useful case label for a manifest validation error.
   id = "<unknown>";
   if isfield(entry, 'case_id')
      id = string(entry.case_id);
   end
end
