function result = packFixtures(version, kwargs)
   %PACKFIXTURES Pack selected release-data capabilities into separate archives.
   %
   %  result = icemodel.verification.setup.packFixtures("v1.1")
   %  result = icemodel.verification.setup.packFixtures("v1.1", ...
   %     capabilities="forcing-integration", root="/tmp/test-data", ...
   %     manifest="/tmp/source-manifest.json", staging_dir="/tmp/stage")
   %
   % Each selected capability becomes the archive named by the release-data
   % manifest. Source bytes must already match their declared size and SHA-256.
   % The generated release manifest adds archive size and SHA-256 metadata and
   % is deterministic apart from the compressed archive bytes themselves.
   %
   % Name-value
   %  capabilities  Capabilities to pack; required v1.1 capabilities by default.
   %  root          Data root containing manifest-declared relative paths.
   %  manifest      Source release-data manifest.
   %  staging_dir   Output directory for archives and generated manifest.
   %  overwrite     Replace existing selected artifacts when true.
   %  silent        Suppress the size summary.
   %
   % See also: icemodel.verification.setup.fixtureFileList,
   %  icemodel.verification.setup.fetchFixtures

   arguments
      version (1, 1) string
      kwargs.capabilities string = defaultCapabilities()
      kwargs.root (1, 1) string = defaultTestDataRoot()
      kwargs.manifest (1, 1) string = defaultManifestFile()
      kwargs.staging_dir (1, 1) string = defaultStagingDir()
      kwargs.overwrite (1, 1) logical = false
      kwargs.silent (1, 1) logical = false
   end

   % Parse selection through the shared manifest gate so list, pack, and fetch
   % cannot disagree about capability membership or install paths.
   [~, selection] = icemodel.verification.setup.fixtureFileList( ...
      capabilities=kwargs.capabilities, root=kwargs.root, ...
      manifest=kwargs.manifest);
   if selection.version ~= version
      error('icemodel:verification:packFixtures:versionMismatch', ...
         'Requested version %s does not match manifest version %s.', ...
         version, selection.version)
   end

   % Refuse empty selected capabilities before touching the staging directory.
   for k = 1:numel(selection.capabilities)
      capability = selection.capabilities(k);
      if ~any(string({selection.files.capability}) == capability)
         error('icemodel:verification:packFixtures:noFiles', ...
            'Capability %s has no declared files to pack.', capability)
      end
   end

   % Verify source bytes against the authoritative file rows before producing
   % any artifact; packaging must never silently bless source drift.
   [missing, mismatched] = verifyFiles(kwargs.root, selection.files);
   if ~isempty(missing) || ~isempty(mismatched)
      error('icemodel:verification:packFixtures:sourceMismatch', ...
         ['Source data do not match the manifest. Missing: %d; ' ...
         'size/hash mismatches: %d.'], numel(missing), numel(mismatched))
   end

   % Resolve all final output paths and clobber policy before writing anything.
   staging_link = icemodel.verification.setup.fixtureCallerSymlink( ...
      kwargs.staging_dir);
   if strlength(staging_link) > 0
      error('icemodel:verification:packFixtures:symlinkDestination', ...
         'Release staging path traverses a symbolic link: %s', staging_link)
   end
   manifest_file = fullfile(kwargs.staging_dir, ...
      "icemodel-" + version + "-data-manifest.json");
   archive_files = strings(numel(selection.archives), 1);
   for k = 1:numel(selection.archives)
      archive_files(k) = fullfile(kwargs.staging_dir, ...
         string(selection.archives(k).name));
   end
   outputs = [archive_files; manifest_file];
   directory_conflicts = outputs(isfolder(outputs));
   if ~isempty(directory_conflicts)
      error('icemodel:verification:packFixtures:destinationTypeConflict', ...
         'Release artifact output is an existing directory: %s', ...
         directory_conflicts(1))
   end
   if ~kwargs.overwrite && any(isfile(outputs))
      error('icemodel:verification:packFixtures:exists', ...
         ['One or more v1.1 data artifacts already exist in %s. ' ...
         'Pass overwrite=true to replace them.'], kwargs.staging_dir)
   end

   % Build every artifact outside the final staging tree. A later archive or
   % manifest failure therefore cannot leave partial outputs that poison retry.
   work_dir = makeWorkDir(kwargs.staging_dir);
   cleaner = onCleanup(@() removeTree(work_dir));
   work_archive_files = fullfile(work_dir, ...
      reshape(string({selection.archives.name}), [], 1));
   work_manifest_file = fullfile(work_dir, ...
      "icemodel-" + version + "-data-manifest.json");

   % Package each capability independently so the optional forcing archive
   % never burdens users who need only formal and showcase data.
   archives = repmat(archiveTemplate(), numel(selection.archives), 1);
   source_bytes = zeros(numel(selection.archives), 1);
   archive_bytes = zeros(numel(selection.archives), 1);
   for k = 1:numel(selection.archives)
      capability = string(selection.archives(k).capability);
      mask = string({selection.files.capability}) == capability;
      entries = selection.files(mask);
      paths = reshape(string({entries.path}), [], 1);
      createArchive(work_archive_files(k), paths, kwargs.root);

      % Record immutable archive metadata only after MATLAB closes the archive.
      info = dir(char(work_archive_files(k)));
      archive_bytes(k) = info.bytes;
      source_bytes(k) = sum(double([entries.bytes]));
      archives(k) = struct( ...
         'capability', capability, ...
         'name', string(selection.archives(k).name), ...
         'required', logical(selection.archives(k).required), ...
         'bytes', info.bytes, ...
         'sha256', icemodel.verification.setup.fileSha256( ...
         work_archive_files(k)));
   end

   % Emit the filtered release manifest beside the selected archives. No
   % timestamp is included, keeping identical inputs byte-stable for review.
   manifest = struct( ...
      'schema_version', selection.schema_version, ...
      'version', version, ...
      'archives', archives, ...
      'files', selection.files);
   writeJson(work_manifest_file, manifest);

   % Promote the complete artifact set as one rollback-protected transaction.
   promoteArtifacts([work_archive_files; work_manifest_file], outputs, ...
      kwargs.staging_dir, work_dir);
   clear cleaner

   % Preserve the historical scalar fields for one-capability callers while
   % exposing vectors for the normal multi-capability v1.1 pack.
   result = struct( ...
      'version', version, ...
      'capabilities', selection.capabilities, ...
      'archive_file', archive_files, ...
      'archive_files', archive_files, ...
      'manifest_file', string(manifest_file), ...
      'manifest', manifest, ...
      'num_files', numel(selection.files), ...
      'source_bytes', sum(source_bytes), ...
      'archive_bytes', archive_bytes, ...
      'saving_bytes', sum(source_bytes) - sum(archive_bytes));

   if ~kwargs.silent
      reportSaving(result);
   end
end

%% Local helpers
function capabilities = defaultCapabilities()
   %DEFAULTCAPABILITIES Required v1.1 capability set.
   capabilities = ["formal-core", "verification-showcase"];
end

function pathname = defaultTestDataRoot()
   %DEFAULTTESTDATAROOT Canonical formal and public-verification data root.
   pathname = string(icemodel.internal.fullpath('test', 'data'));
end

function pathname = defaultManifestFile()
   %DEFAULTMANIFESTFILE Tracked authoritative release-data manifest.
   pathname = string(icemodel.internal.fullpath('test', 'assets', ...
      'icemodel-v1.1-data-manifest.json'));
end

function pathname = defaultStagingDir()
   %DEFAULTSTAGINGDIR Canonical gitignored release artifact directory.
   pathname = string(icemodel.internal.fullpath('release-staging'));
end

function row = archiveTemplate()
   %ARCHIVETEMPLATE Stable field order for generated archive rows.
   row = struct('capability', "", 'name', "", 'required', false, ...
      'bytes', 0, 'sha256', "");
end

function work_dir = makeWorkDir(staging_dir)
   %MAKEWORKDIR Create transaction storage without creating the staging tree.
   parent = string(fileparts(staging_dir));
   if strlength(parent) == 0
      parent = string(pwd);
   end
   while ~isfolder(parent)
      next = string(fileparts(parent));
      if strlength(next) == 0 || next == parent
         parent = string(pwd);
         break
      end
      parent = next;
   end
   work_dir = string(tempname(parent));
   mkdir(work_dir)
end

function promoteArtifacts(sources, destinations, staging_dir, work_dir)
   %PROMOTEARTIFACTS Commit all packed outputs or restore their prior state.
   backup_root = fullfile(work_dir, "backup");
   mkdir(backup_root)
   max_outputs = numel(destinations);
   backups = repmat( ...
      struct('destination', "", 'backup', ""), max_outputs, 1);
   promoted = strings(max_outputs, 1);
   created_dirs = strings(count(staging_dir, filesep) + 1, 1);
   backup_count = 0;
   promoted_count = 0;
   created_count = 0;
   try
      [created_dirs, created_count] = ensureParentTracked( ...
         staging_dir, created_dirs, created_count);
      for k = 1:numel(destinations)
         destination = destinations(k);
         backup = fullfile(backup_root, "output-" + string(k));

         % Preserve every replaceable pre-call output until the full commit.
          if isfolder(destination)
             error( ...
                'icemodel:verification:packFixtures:destinationTypeConflict', ...
                'Release artifact output is an existing directory: %s', ...
                destination)
          elseif isfile(destination)
             moveChecked(destination, backup);
            backup_count = backup_count + 1;
            backups(backup_count) = struct( ...
               'destination', destination, 'backup', backup);
         end
         moveChecked(sources(k), destination);
         promoted_count = promoted_count + 1;
         promoted(promoted_count) = destination;
      end
   catch cause
      try
         rollbackArtifacts(promoted(1:promoted_count), ...
            backups(1:backup_count), created_dirs(1:created_count));
      catch rollback_cause
         exception = MException( ...
            'icemodel:verification:packFixtures:rollbackFailed', ...
            'Packing failed and rollback could not restore prior artifacts.');
         exception = addCause(exception, cause);
         exception = addCause(exception, rollback_cause);
         throw(exception)
      end
      rethrow(cause)
   end
end

function [created_dirs, created_count] = ensureParentTracked( ...
      parent, created_dirs, created_count)
   %ENSUREPARENTTRACKED Create a directory chain and journal prior absence.
   missing = strings(count(parent, filesep) + 1, 1);
   missing_count = 0;
   cursor = string(parent);
   while ~isfolder(cursor)
      if isfile(cursor)
         error('icemodel:verification:packFixtures:promotionFailed', ...
            'Artifact parent is a file: %s', cursor)
      end
      missing_count = missing_count + 1;
      missing(missing_count) = cursor;
      next = string(fileparts(cursor));
      if strlength(next) == 0 || next == cursor
         break
      end
      cursor = next;
   end
   for k = missing_count:-1:1
      [ok, message, message_id] = mkdir(missing(k));
      if ~ok
         error('icemodel:verification:packFixtures:promotionFailed', ...
            'Cannot create %s (%s: %s).', ...
            missing(k), message_id, message)
      end
      created_count = created_count + 1;
      created_dirs(created_count) = missing(k);
   end
end

function rollbackArtifacts(promoted, backups, created_dirs)
   %ROLLBACKARTIFACTS Remove this call's outputs and restore all prior entries.
   for k = numel(promoted):-1:1
      removePath(promoted(k));
   end
   for k = numel(backups):-1:1
      moveChecked(backups(k).backup, backups(k).destination);
   end
   for k = numel(created_dirs):-1:1
      if isfolder(created_dirs(k))
         contents = dir(char(created_dirs(k)));
         names = string({contents.name});
         contents = contents(~ismember(names, [".", ".."]));
         if ~isempty(contents)
            error('icemodel:verification:packFixtures:rollbackFailed', ...
               'New artifact directory is not empty: %s', created_dirs(k))
         end
         [ok, message, message_id] = rmdir(created_dirs(k));
         if ~ok
            error('icemodel:verification:packFixtures:rollbackFailed', ...
               'Cannot remove %s (%s: %s).', ...
               created_dirs(k), message_id, message)
         end
      end
   end
end

function moveChecked(source, destination)
   %MOVECHECKED Move one artifact and retain the platform failure details.
   [ok, message, message_id] = movefile(source, destination, 'f');
   if ~ok
      error('icemodel:verification:packFixtures:promotionFailed', ...
         'Cannot move %s to %s (%s: %s).', ...
         source, destination, message_id, message)
   end
end

function removePath(pathname)
   %REMOVEPATH Delete one transaction-created file or directory.
   if isfile(pathname)
      delete(pathname)
   elseif isfolder(pathname)
      rmdir(pathname, 's')
   end
end

function removeTree(pathname)
   %REMOVETREE Delete isolated transaction storage after success or failure.
   if isfolder(pathname)
      rmdir(pathname, 's')
   end
end

function [missing, mismatched] = verifyFiles(root, entries)
   %VERIFYFILES Verify declared source files without accepting extra files.
   paths = reshape(string({entries.path}), [], 1);
   missing_mask = false(numel(paths), 1);
   mismatched_mask = false(numel(paths), 1);
   for k = 1:numel(entries)
      relpath = paths(k);
      assertNoSourceSymlinks(root, relpath);
      pathname = underRoot(root, relpath);
      if ~isfile(pathname)
         missing_mask(k) = true;
         continue
      end
      info = dir(char(pathname));
      digest = icemodel.verification.setup.fileSha256(pathname);
      if info.bytes ~= entries(k).bytes ...
            || ~strcmpi(digest, string(entries(k).sha256))
         mismatched_mask(k) = true;
      end
   end
   missing = paths(missing_mask);
   mismatched = paths(mismatched_mask);
end

function assertNoSourceSymlinks(root, relpath)
   %ASSERTNOSOURCESYMLINKS Reject links before reading manifest-declared bytes.
   ancestor_link = icemodel.verification.setup.fixtureCallerSymlink( ...
      string(root));
   if strlength(ancestor_link) > 0
      error('icemodel:verification:packFixtures:symlinkSource', ...
         'Manifest source root traverses a symbolic link: %s', ancestor_link)
   end
   cursor = icemodel.verification.setup.fixtureCanonicalRoot(string(root));

   % Inspect every manifest component inside the canonical source root.
   parts = split(relpath, "/");
   for k = 1:numel(parts)
      cursor = fullfile(cursor, parts(k));
      assertNotSymbolicLink(cursor, relpath);
   end
end

function assertNotSymbolicLink(pathname, relpath)
   %ASSERTNOTSYMBOLICLINK Keep packing inside the selected source tree.
   path_object = java.nio.file.Paths.get(char(pathname), ...
      javaArray('java.lang.String', 0));
   if java.nio.file.Files.isSymbolicLink(path_object)
      error('icemodel:verification:packFixtures:symlinkSource', ...
         'Manifest source path traverses a symbolic link: %s', relpath)
   end
end

function createArchive(archive_file, paths, root)
   %CREATEARCHIVE Create archives without undeclared macOS metadata members.
   if ismac
      % MATLAB's tar writer serializes com.apple.* attributes as AppleDouble
      % files. Native bsdtar also normalizes ownership and omits the gzip clock
      % while preserving source bytes, modes, mtimes, and declared path order.
      members = [ancestorDirectories(paths); paths];
      quoted_members = strings(numel(members), 1);
      for k = 1:numel(members)
         quoted_members(k) = shellQuote(members(k));
      end
      command = "COPYFILE_DISABLE=1 /usr/bin/tar --format ustar " ...
         + "--uid 0 --gid 0 --uname root --gname root " ...
         + "--no-recursion --options gzip:!timestamp -czf " ...
         + shellQuote(archive_file) + " -C " + shellQuote(root) + " -- " ...
         + strjoin(quoted_members, " ");
      [status, message] = system(char(command));
      if status ~= 0
         error('icemodel:verification:packFixtures:archiveFailed', ...
            'Cannot create release-data archive: %s', strtrim(message))
      end
      return
   end

   % MATLAB tar is portable and does not synthesize AppleDouble files away
   % from macOS, so it remains the dependency-free implementation elsewhere.
   tar(char(archive_file), cellstr(paths), char(root));
end

function directories = ancestorDirectories(files)
   %ANCESTORDIRECTORIES Derive canonical directory entries for archive headers.
   directories = strings(sum(count(files, "/")), 1);
   cursor = 0;
   for k = 1:numel(files)
      parts = split(files(k), "/");
      for n = 1:numel(parts) - 1
         cursor = cursor + 1;
         directories(cursor) = strjoin(parts(1:n), "/");
      end
   end
   directories = unique(directories);
end

function quoted = shellQuote(pathname)
   %SHELLQUOTE Quote one path for the fixed native macOS tar command.
   quote = char(39);
   escaped = strrep(char(pathname), quote, [quote '"' quote '"' quote]);
   quoted = string([quote escaped quote]);
end

function pathname = underRoot(root, relpath)
   %UNDERROOT Convert a validated relative POSIX path to a local path.
   parts = cellstr(split(relpath, "/"));
   pathname = string(fullfile(root, parts{:}));
end

function writeJson(pathname, data)
   %WRITEJSON Write deterministic pretty JSON with array shape preserved.
   payload = data;
   payload.archives = num2cell(payload.archives);
   payload.files = num2cell(payload.files);
   fid = fopen(pathname, 'w');
   if fid < 0
      error('icemodel:verification:packFixtures:cannotWrite', ...
         'Cannot write release-data manifest: %s', pathname)
   end
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, jsonencode(payload, PrettyPrint=true), 'char');
end

function reportSaving(result)
   %REPORTSAVING Print one compact release-data packing summary.
   mb = @(bytes) double(bytes) / 1024 / 1024;
   fprintf('\n=== icemodel v1.1 release data packed ===\n');
   fprintf('Version:       %s\n', result.version);
   fprintf('Capabilities:  %s\n', strjoin(result.capabilities, ', '));
   fprintf('Files:         %d\n', result.num_files);
   fprintf('Source:        %.1f MB\n', mb(result.source_bytes));
   fprintf('Archives:      %.1f MB\n', mb(sum(result.archive_bytes)));
   fprintf('Manifest:      %s\n\n', result.manifest_file);
end
