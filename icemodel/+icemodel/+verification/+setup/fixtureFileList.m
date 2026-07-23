function [files, selection] = fixtureFileList(kwargs)
   %FIXTUREFILELIST Return manifest paths for selected release capabilities.
   %
   %  files = icemodel.verification.setup.fixtureFileList()
   %  files = icemodel.verification.setup.fixtureFileList( ...
   %     capabilities="forcing-integration", root="/tmp/data", ...
   %     manifest="/tmp/icemodel-v1.1-data-manifest.json")
   %  [files, selection] = icemodel.verification.setup.fixtureFileList(...)
   %
   % The tracked release-data manifest is the single source of truth shared by
   % packFixtures and fetchFixtures. Returned paths are relative POSIX paths
   % below root. The optional second output contains the selected manifest rows.
   %
   % Name-value
   %  capabilities  Capability names to select. The two required v1.1
   %                capabilities are selected by default.
   %  root          Installation root associated with the returned paths.
   %  manifest      Release-data manifest JSON file.
   %
   % See also: icemodel.verification.setup.packFixtures,
   %  icemodel.verification.setup.fetchFixtures

   arguments
      kwargs.capabilities string = defaultCapabilities()
      kwargs.root (1, 1) string = defaultTestDataRoot()
      kwargs.manifest (1, 1) string = defaultManifestFile()
   end

   % Load and validate the shared manifest before filtering so malformed rows
   % cannot be hidden by a narrower capability selection.
   manifest = readManifest(kwargs.manifest);
   archives = structRows(manifest.archives, "archives");
   entries = structRows(manifest.files, "files");
   validateArchiveRows(archives);
   validateFileRows(entries);
   validateCapabilityMembership(entries, archives);

   % Normalize and validate the requested set against archive declarations;
   % archives define capabilities even before optional file rows are staged.
   requested = unique(reshape(string(kwargs.capabilities), 1, []), 'stable');
   requested = requested(strlength(requested) > 0);
   if isempty(requested)
      error('icemodel:verification:fixtureFileList:noCapabilities', ...
         'Select at least one release-data capability.')
   end
   known = string({archives.capability});
   unknown = requested(~ismember(requested, known));
   if ~isempty(unknown)
      error('icemodel:verification:fixtureFileList:unknownCapability', ...
         'Unknown release-data capability: %s', strjoin(unknown, ', '))
   end

   % Canonical capability and path order makes filtered manifests and archive
   % mapping independent of caller request order. The manifest declaration
   % order is authoritative rather than lexical sorting.
   selected_mask = ismember(known, requested);
   capabilities = known(selected_mask);
   selected_archives = archives(selected_mask);
   if isempty(entries)
      selected_entries = entries;
   else
      selected_entries = entries(ismember(string({entries.capability}), ...
         capabilities));
   end
   files = reshape(string({selected_entries.path}), [], 1);

   % Return the parsed selection so pack/fetch reuse this exact schema gate.
   selection = struct( ...
      'schema_version', manifest.schema_version, ...
      'version', string(manifest.version), ...
      'capabilities', reshape(capabilities, [], 1), ...
      'root', kwargs.root, ...
      'manifest_file', kwargs.manifest, ...
      'archives', selected_archives, ...
      'files', selected_entries);
end

%% Local helpers
function capabilities = defaultCapabilities()
   %DEFAULTCAPABILITIES Required v1.1 capabilities installed together.
   capabilities = ["formal-core", "verification-showcase"];
end

function pathname = defaultTestDataRoot()
   %DEFAULTTESTDATAROOT Canonical release-provisioned formal-data root.
   pathname = string(icemodel.internal.fullpath('test', 'data'));
end

function pathname = defaultManifestFile()
   %DEFAULTMANIFESTFILE Tracked authoritative release-data manifest.
   pathname = string(icemodel.internal.fullpath('test', 'assets', ...
      'icemodel-v1.1-data-manifest.json'));
end

function manifest = readManifest(pathname)
   %READMANIFEST Decode the release-data manifest with a stable error contract.
   if ~isfile(pathname)
      error('icemodel:verification:fixtureFileList:noManifest', ...
         'Release-data manifest not found: %s', pathname)
   end
   try
      manifest = jsondecode(fileread(pathname));
   catch cause
      errorWithCause('icemodel:verification:fixtureFileList:invalidManifest', ...
         cause, 'Cannot decode release-data manifest: %s', pathname)
   end
   required = ["schema_version", "version", "archives", "files"];
   if ~isstruct(manifest) || ~all(isfield(manifest, required))
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Manifest must define schema_version, version, archives, and files.')
   end
   version = string(manifest.version);
   valid_version = ~isempty(regexp(char(version), ...
      '^v[0-9][A-Za-z0-9._-]*$', 'once'));
   if manifest.schema_version ~= 1 || ~valid_version
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Manifest schema_version must be 1 and version must be a safe v-prefixed tag.')
   end
end

function rows = structRows(value, field_name)
   %STRUCTROWS Normalize decoded JSON object arrays to struct columns.
   if isempty(value)
      rows = repmat(struct(), 0, 1);
      return
   end
   if ~isstruct(value)
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Manifest field %s must be an array of objects.', field_name)
   end
   rows = reshape(value, [], 1);
end

function validateArchiveRows(rows)
   %VALIDATEARCHIVEROWS Validate capability and archive declarations.
   required = ["capability", "name", "required", "bytes", "sha256"];
   if isempty(rows) || ~all(isfield(rows, required))
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Manifest archives must declare capability, name, required, bytes, and sha256.')
   end
   capabilities = string({rows.capability});
   names = string({rows.name});
   valid_capabilities = arrayfun(@(value) ~isempty(regexp(char(value), ...
      '^[a-z0-9][a-z0-9-]*$', 'once')), capabilities);
   valid_names = arrayfun(@isSafeArchiveName, names);
   if any(~valid_capabilities) || any(~valid_names) ...
         || numel(unique(capabilities)) ~= numel(capabilities) ...
         || numel(unique(names)) ~= numel(names)
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Archive capability and filename values must be safe and unique.')
   end
   for k = 1:numel(rows)
      if ~isscalar(rows(k).required) || ~islogical(rows(k).required) ...
            || ~isscalar(rows(k).bytes) || ~isnumeric(rows(k).bytes) ...
            || rows(k).bytes < 0
         error('icemodel:verification:fixtureFileList:invalidManifest', ...
            'Archive %s has invalid required or bytes metadata.', names(k))
      end
      validateOptionalDigest(string(rows(k).sha256), names(k));
   end
end

function tf = isSafeArchiveName(name)
   %ISSAFEARCHIVENAME Require one relative tar.gz filename without controls.
   has_controls = ~isempty(regexp(char(name), '[\x00-\x1F\x7F]', 'once'));
   backslash = string(char(92));
   tf = strlength(name) > 7 && endsWith(name, ".tar.gz") ...
      && ~contains(name, "/") && ~contains(name, backslash) && ~has_controls;
end

function validateFileRows(rows)
   %VALIDATEFILEROWS Validate install paths and immutable file metadata.
   if isempty(rows)
      return
   end
   required = ["capability", "path", "required", "bytes", "sha256"];
   if ~all(isfield(rows, required))
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Manifest files must declare capability, path, required, bytes, and sha256.')
   end
   paths = string({rows.path});
   if numel(unique(paths)) ~= numel(paths)
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Manifest file paths must be unique across capabilities.')
   end
   if ~isequal(paths, sort(paths))
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Manifest file paths must use canonical ascending order.')
   end
   for k = 1:numel(rows)
      validateRelativePath(paths(k));
      if strlength(string(rows(k).capability)) == 0 ...
            || ~isscalar(rows(k).required) || ~islogical(rows(k).required) ...
            || ~isscalar(rows(k).bytes) || ~isnumeric(rows(k).bytes) ...
            || rows(k).bytes < 0
         error('icemodel:verification:fixtureFileList:invalidManifest', ...
            'Manifest file %s has invalid metadata.', paths(k))
      end
      validateDigest(string(rows(k).sha256), paths(k));
   end
end

function validateCapabilityMembership(entries, archives)
   %VALIDATECAPABILITYMEMBERSHIP Tie every file row to one archive policy.
   if isempty(entries)
      return
   end
   archive_capabilities = string({archives.capability});
   file_capabilities = string({entries.capability});
   if any(~ismember(file_capabilities, archive_capabilities))
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Every file capability must have one matching archive declaration.')
   end
   for k = 1:numel(entries)
      archive = archives(archive_capabilities == file_capabilities(k));
      if logical(entries(k).required) ~= logical(archive.required)
         error('icemodel:verification:fixtureFileList:invalidManifest', ...
            'File %s required status must match its capability archive.', ...
            entries(k).path)
      end
   end
end

function validateRelativePath(pathname)
   %VALIDATERELATIVEPATH Reject paths that could escape the selected root.
   parts = split(pathname, "/");
   is_drive = ~isempty(regexp(char(pathname), '^[A-Za-z]:', 'once'));
   has_controls = ~isempty(regexp(char(pathname), ...
      '[\x00-\x1F\x7F]', 'once'));
   backslash = string(char(92));
   if strlength(pathname) == 0 || startsWith(pathname, "/") ...
         || contains(pathname, backslash) || is_drive || has_controls ...
         || any(parts == ".." | parts == "." | strlength(parts) == 0)
      error('icemodel:verification:fixtureFileList:unsafePath', ...
         'Manifest path must be a safe relative POSIX path: %s', pathname)
   end
end

function validateDigest(digest, label)
   %VALIDATEDIGEST Require a complete SHA-256 for every declared data file.
   if isempty(regexp(char(digest), '^[0-9A-Fa-f]{64}$', 'once'))
      error('icemodel:verification:fixtureFileList:invalidManifest', ...
         'Manifest file %s must declare a 64-character SHA-256.', label)
   end
end

function validateOptionalDigest(digest, label)
   %VALIDATEOPTIONALDIGEST Allow unstaged archives but validate populated hashes.
   if strlength(digest) > 0
      validateDigest(digest, label);
   end
end

function errorWithCause(identifier, cause, message, varargin)
   %ERRORWITHCAUSE Throw a stable public error while retaining root-cause detail.
   exception = MException(identifier, message, varargin{:});
   exception = addCause(exception, cause);
   throw(exception)
end
