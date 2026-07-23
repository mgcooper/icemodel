function result = fetchFixtures(version, kwargs)
   %FETCHFIXTURES Transactionally provision or verify v1.1 release data.
   %
   %  result = icemodel.verification.setup.fetchFixtures("v1.1")
   %  result = icemodel.verification.setup.fetchFixtures("v1.1", ...
   %     capabilities="formal-core", archive="/tmp/formal-core.tar.gz")
   %  result = icemodel.verification.setup.fetchFixtures("v1.1", ...
   %     download=false)
   %
   % Calling this provisioning API is the explicit request to download missing
   % selected archives; the default selection is the mandatory v1.1 set. Pass
   % download=false for network-free verification. Local archive and manifest
   % overrides keep offline and pre-publication workflows fully supported.
   %
   % Name-value
   %  capabilities  Capability names; required v1.1 set by default.
   %  root          Destination data root.
   %  manifest      Tracked, staged, or local release-data manifest.
   %  archive       Local archives in canonical capability order as returned by
   %                fixtureFileList or packFixtures. A scalar archive requires
   %                exactly one selected capability.
   %  download      Download missing archives from the release; true by default.
   %  release_url   Optional release-asset base URL override.
   %  extract       Provision archives when true; verify-only when false.
   %  repo          GitHub owner/repository used by the default release URL.
   %  strict        Throw when installed data are incomplete.
   %  silent        Suppress the actionable provisioning banner.
   %
   % See also: icemodel.verification.setup.fixtureFileList,
   %  icemodel.verification.setup.packFixtures

   arguments
      version (1, 1) string = "v1.1"
      kwargs.capabilities string = defaultCapabilities()
      kwargs.root (1, 1) string = defaultTestDataRoot()
      kwargs.manifest (1, 1) string = defaultManifestFile()
      kwargs.archive string = strings(1, 0)
      kwargs.download (1, 1) logical = true
      kwargs.release_url (1, 1) string = ""
      kwargs.extract (1, 1) logical = true
      kwargs.repo (1, 1) string = "mgcooper/icemodel"
      kwargs.strict (1, 1) logical = true
      kwargs.silent (1, 1) logical = false
   end

   % Resolve capability rows through the same schema gate used for packing.
   [~, selection] = icemodel.verification.setup.fixtureFileList( ...
      capabilities=kwargs.capabilities, root=kwargs.root, ...
      manifest=kwargs.manifest);
   if selection.version ~= version
      error('icemodel:verification:fetchFixtures:versionMismatch', ...
         'Requested version %s does not match manifest version %s.', ...
         version, selection.version)
   end

   % Reject destination links before even the verified fast path can inspect
   % matching bytes through a path that escapes the selected data root.
   assertNoDestinationSymlinks(kwargs.root, {selection});

   % A complete installed capability is an idempotent success. This branch
   % returns before archive resolution, download, extraction, or writes.
   [missing, mismatched] = verifyInstalled(kwargs.root, selection);
   if isempty(missing) && isempty(mismatched)
      result = resultStruct(true, "verified", kwargs.root, selection, ...
         missing, mismatched, "");
      return
   end

   % Network-free verification calls report one exact provisioning command
   % rather than mutating or reaching the network.
   local_archives = reshape(string(kwargs.archive), [], 1);
   local_archives = local_archives(strlength(local_archives) > 0);
   if ~kwargs.extract || (isempty(local_archives) && ~kwargs.download)
      result = unavailableResult(version, kwargs.root, selection, missing, ...
         mismatched, kwargs.release_url, kwargs.repo, kwargs.strict, ...
         kwargs.silent);
      return
   end

   % Local archives map one-to-one to the selected capability order. This
   % prevents one archive from being guessed across multiple capabilities.
   if ~isempty(local_archives) ...
         && numel(local_archives) ~= numel(selection.capabilities)
      error('icemodel:verification:fetchFixtures:archiveCapabilityMismatch', ...
         ['Provide exactly one local archive per selected capability. ' ...
         'A scalar archive requires exactly one capability.'])
   end

   % Resolve and validate every missing capability before the first destination
   % write. One work tree and one later rollback journal cover the whole call.
   work_dir = makeWorkDir(kwargs.root);
   cleaner = onCleanup(@() removeTree(work_dir));
   pending = cell(numel(selection.capabilities), 1);
   extracted_roots = strings(numel(selection.capabilities), 1);
   n_pending = 0;
   for k = 1:numel(selection.capabilities)
      capability = selection.capabilities(k);
      selected = capabilitySelection(selection, capability);
      [cap_missing, cap_mismatched] = verifyInstalled(kwargs.root, selected);
      if isempty(cap_missing) && isempty(cap_mismatched)
         continue
      end

      % Keep each archive/extracted tree isolated while retaining them all until
      % the invocation-wide promotion either commits or rolls back.
      n_pending = n_pending + 1;
      capability_work_dir = fullfile(work_dir, "capability-" + n_pending);
      mkdir(capability_work_dir)
      if isempty(local_archives)
         archive_file = downloadArchive(selected.archives, version, ...
            kwargs.repo, kwargs.release_url, capability_work_dir);
      else
         archive_file = local_archives(k);
         if ~isfile(archive_file)
            error('icemodel:verification:fetchFixtures:noArchive', ...
               'Archive not found: %s', archive_file)
         end
      end

      % Validation completes for all selected capabilities before promotion.
      extracted_roots(n_pending) = validateAndExtract( ...
         archive_file, selected, capability_work_dir);
      pending{n_pending} = selected;
   end

   % One shared backup journal makes all changed capabilities commit or restore
   % together, including failure in a later capability or final verification.
   pending = pending(1:n_pending);
   extracted_roots = extracted_roots(1:n_pending);
   promoteCapabilities(extracted_roots, kwargs.root, pending, ...
      selection, work_dir);
   clear cleaner
   result = resultStruct(true, "provisioned", kwargs.root, selection, ...
      strings(0, 1), strings(0, 1), "");
end

%% Selection and status helpers
function capabilities = defaultCapabilities()
   %DEFAULTCAPABILITIES Required v1.1 capabilities installed together.
   capabilities = ["formal-core", "verification-showcase"];
end

function pathname = defaultTestDataRoot()
   %DEFAULTTESTDATAROOT Canonical release-provisioned test data root.
   pathname = string(icemodel.internal.fullpath('test', 'data'));
end

function pathname = defaultManifestFile()
   %DEFAULTMANIFESTFILE Tracked authoritative release-data manifest.
   pathname = string(icemodel.internal.fullpath('test', 'assets', ...
      'icemodel-v1.1-data-manifest.json'));
end

function selected = capabilitySelection(selection, capability)
   %CAPABILITYSELECTION Filter a parsed selection to one transaction boundary.
   selected = selection;
   selected.capabilities = capability;
   selected.archives = selection.archives( ...
      string({selection.archives.capability}) == capability);
   if isempty(selection.files)
      selected.files = selection.files;
   else
      selected.files = selection.files( ...
         string({selection.files.capability}) == capability);
   end
end

function [missing, mismatched] = verifyInstalled(root, selection)
   %VERIFYINSTALLED Verify selected declared files without rejecting unrelated data.
   if isempty(selection.files)
      missing = "capability:" + selection.capabilities;
      mismatched = strings(0, 1);
      return
   end

   % Preallocate one slot per declaration, then compact after verification.
   paths = reshape(string({selection.files.path}), [], 1);
   missing_mask = false(numel(paths), 1);
   mismatched_mask = false(numel(paths), 1);
   for k = 1:numel(selection.files)
      entry = selection.files(k);
      pathname = underRoot(root, paths(k));
      if ~isfile(pathname)
         missing_mask(k) = true;
         continue
      end
      info = dir(char(pathname));
      digest = icemodel.verification.setup.fileSha256(pathname);
      if info.bytes ~= entry.bytes ...
            || ~strcmpi(digest, string(entry.sha256))
         mismatched_mask(k) = true;
      end
   end
   missing = paths(missing_mask);
   mismatched = paths(mismatched_mask);
end

function result = unavailableResult(version, root, selection, missing, ...
      mismatched, release_url, repo, strict, silent)
   %UNAVAILABLERESULT Report exact explicit provisioning without downloading.
   command = fetchCommand(version, selection.capabilities, root, ...
      selection.manifest_file, release_url, repo);
   if ~silent
      fprintf('\n=== icemodel required release data incomplete ===\n');
      fprintf('Data root: %s\n', root);
      fprintf('Missing: %d; mismatched: %d\n', ...
         numel(missing), numel(mismatched));
      fprintf('Provision explicitly with:\n  %s\n\n', command);
   end
   if strict
      error('icemodel:verification:fetchFixtures:incompleteFixtures', ...
         ['Required release data are incomplete under %s. Run exactly: ' ...
         '%s'], root, command)
   end
   result = resultStruct(false, "missing", root, selection, missing, ...
      mismatched, command);
end

function command = fetchCommand(version, capabilities, root, manifest, ...
      release_url, repo)
   %FETCHCOMMAND Build the stable copy-paste provisioning command.
   command = icemodel.verification.setup.fixtureFetchCommand( ...
      version, capabilities, root=root, manifest=manifest, ...
      release_url=release_url, repo=repo);
end

function result = resultStruct(ok, mode, root, selection, missing, ...
      mismatched, command)
   %RESULTSTRUCT Return one stable status shape for verification and provisioning.
   result = struct( ...
      'ok', ok, ...
      'mode', mode, ...
      'root', root, ...
      'manifest_file', selection.manifest_file, ...
      'capabilities', selection.capabilities, ...
      'missing', missing, ...
      'mismatched', mismatched, ...
      'command', command);
end

%% Archive validation and extraction
function work_dir = makeWorkDir(root)
   %MAKEWORKDIR Create transaction storage without creating destination parents.
   parent = string(fileparts(root));
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

function archive_file = downloadArchive(archive, version, repo, ...
      release_url, work_dir)
   %DOWNLOADARCHIVE Download one explicitly requested release asset.
   validateArchiveMetadata(archive);
   if strlength(release_url) == 0
      release_url = "https://github.com/" + repo ...
         + "/releases/download/" + version;
   end
   release_url = strip(release_url, 'right', '/');
   archive_file = fullfile(work_dir, string(archive.name));
   url = release_url + "/" + string(archive.name);
   try
      websave(char(archive_file), char(url));
   catch cause
      throwWithCause('icemodel:verification:fetchFixtures:downloadFailed', ...
         cause, 'Failed to download release-data archive: %s', url)
   end
end

function extracted_root = validateAndExtract(archive_file, selection, work_dir)
   %VALIDATEANDEXTRACT Validate hash and tar headers before isolated extraction.
   validateArchiveMetadata(selection.archives);
   info = dir(char(archive_file));
   digest = icemodel.verification.setup.fileSha256(archive_file);
   if info.bytes ~= selection.archives.bytes
      error('icemodel:verification:fetchFixtures:archiveSizeMismatch', ...
         'Archive byte size does not match the manifest: %s', archive_file)
   end
   if ~strcmpi(digest, string(selection.archives.sha256))
      error('icemodel:verification:fetchFixtures:archiveHashMismatch', ...
         'Archive SHA-256 does not match the manifest: %s', archive_file)
   end

   % Decompress only into transaction storage, then inspect every raw tar
   % header before invoking untar on untrusted member names or types.
   header_dir = fullfile(work_dir, "headers");
   mkdir(header_dir)
   try
      tar_files = gunzip(char(archive_file), char(header_dir));
   catch cause
      throwWithCause('icemodel:verification:fetchFixtures:invalidArchive', ...
         cause, 'Cannot decompress archive: %s', archive_file)
   end
   if numel(tar_files) ~= 1
      error('icemodel:verification:fetchFixtures:invalidArchive', ...
         'Expected one tar stream in archive: %s', archive_file)
   end
   [member_files, member_dirs] = tarMembers(string(tar_files{1}));
   validateMembers(member_files, member_dirs, selection.files);

   % Extract only after headers pass, then reject missing/extra files and verify
   % every extracted byte before promotion.
   extracted_root = fullfile(work_dir, "extracted");
   mkdir(extracted_root)
   try
      untar(char(archive_file), char(extracted_root));
   catch cause
      throwWithCause('icemodel:verification:fetchFixtures:invalidArchive', ...
         cause, 'Cannot extract validated archive: %s', archive_file)
   end
   actual = inventoryFiles(extracted_root);
   expected = sort(reshape(string({selection.files.path}), [], 1));
   missing = setdiff(expected, actual);
   extra = setdiff(actual, expected);
   if ~isempty(missing)
      error('icemodel:verification:fetchFixtures:missingMember', ...
         'Archive is missing declared file: %s', missing(1))
   end
   if ~isempty(extra)
      error('icemodel:verification:fetchFixtures:extraMember', ...
         'Archive contains undeclared file: %s', extra(1))
   end
   [missing, mismatched] = verifyInstalled(extracted_root, selection);
   if ~isempty(missing) || ~isempty(mismatched)
      error('icemodel:verification:fetchFixtures:fileHashMismatch', ...
         'Extracted files do not match declared size and SHA-256 metadata.')
   end
end

function validateArchiveMetadata(archive)
   %VALIDATEARCHIVEMETADATA Require immutable metadata before trusting bytes.
   digest = string(archive.sha256);
   valid_digest = ~isempty(regexp(char(digest), ...
      '^[0-9A-Fa-f]{64}$', 'once'));
   if ~isscalar(archive.bytes) || ~isnumeric(archive.bytes) ...
         || archive.bytes <= 0 || ~valid_digest
      error('icemodel:verification:fetchFixtures:incompleteManifest', ...
         'Archive %s lacks final byte-size or SHA-256 metadata.', archive.name)
   end
end

function [files, directories] = tarMembers(tar_file)
   %TARMEMBERS Read tar headers and reject non-file/non-directory entries.
   fid = fopen(tar_file, 'r');
   if fid < 0
      error('icemodel:verification:fetchFixtures:invalidArchive', ...
         'Cannot inspect tar stream: %s', tar_file)
   end
   cleaner = onCleanup(@() fclose(fid));

   % Count valid members first so the second pass can allocate exact outputs.
   file_count = 0;
   directory_count = 0;
   while true
      [done, ~, kind] = nextTarMember(fid, tar_file);
      if done
         break
      end
      file_count = file_count + (kind == 1);
      directory_count = directory_count + (kind == 2);
   end
   if fseek(fid, 0, 'bof') ~= 0
      error('icemodel:verification:fetchFixtures:invalidArchive', ...
         'Cannot rewind tar stream: %s', tar_file)
   end

   % Read the same validated stream into its exact-size inventories.
   files = strings(file_count, 1);
   directories = strings(directory_count, 1);
   file_index = 0;
   directory_index = 0;
   while true
      [done, name, kind] = nextTarMember(fid, tar_file);
      if done
         break
      end
      if kind == 1
         file_index = file_index + 1;
         files(file_index) = name;
      elseif kind == 2
         directory_index = directory_index + 1;
         directories(directory_index) = name;
      end
   end
   if numel(unique(files)) ~= numel(files) ...
         || numel(unique(directories)) ~= numel(directories)
      error('icemodel:verification:fetchFixtures:unsafeMember', ...
         'Archive contains duplicate file or directory entries.')
   end
end

function [done, name, kind] = nextTarMember(fid, tar_file)
   %NEXTTARMEMBER Read, validate, and advance past one tar member.
   header = fread(fid, 512, '*uint8')';
   done = isempty(header) || all(header == 0);
   name = "";
   kind = 0;
   if done
      return
   end
   if numel(header) ~= 512
      error('icemodel:verification:fetchFixtures:invalidArchive', ...
         'Truncated tar header in %s.', tar_file)
   end

   % Decode and normalize the member before accepting its type.
   name = tarText(header(1:100));
   prefix = tarText(header(346:500));
   if strlength(prefix) > 0
      name = prefix + "/" + name;
   end
   name = normalizeMember(name);
   type = header(157);
   if type == 0 || type == uint8('0')
      kind = 1;
   elseif type == uint8('5') && strlength(name) > 0
      kind = 2;
   elseif type ~= uint8('5')
      error('icemodel:verification:fetchFixtures:unsafeMember', ...
         'Archive member %s has forbidden tar type %d.', name, type)
   end

   % Advance over the padded payload to the next header.
   size = tarOctal(header(125:136));
   if fseek(fid, ceil(double(size) / 512) * 512, 'cof') ~= 0
      error('icemodel:verification:fetchFixtures:invalidArchive', ...
         'Cannot seek over tar member %s.', name)
   end
end

function text = tarText(bytes)
   %TARTEXT Decode one null-terminated ASCII tar header field.
   stop = find(bytes == 0, 1) - 1;
   if isempty(stop)
      stop = numel(bytes);
   end
   text = string(char(bytes(1:stop)));
end

function value = tarOctal(bytes)
   %TAROCTAL Parse a traditional ASCII octal tar size field.
   text = strtrim(char(bytes(bytes ~= 0)));
   if isempty(text)
      value = 0;
      return
   end
   if isempty(regexp(text, '^[0-7]+$', 'once'))
      error('icemodel:verification:fetchFixtures:invalidArchive', ...
         'Tar member size is not an octal integer.')
   end
   value = base2dec(text, 8);
end

function name = normalizeMember(name)
   %NORMALIZEMEMBER Normalize benign ./ prefixes and reject unsafe paths.
   while startsWith(name, "./")
      name = extractAfter(name, 2);
   end
   if endsWith(name, "/")
      name = extractBefore(name, strlength(name));
   end
   parts = split(name, "/");
   is_drive = ~isempty(regexp(char(name), '^[A-Za-z]:', 'once'));
   backslash = string(char(92));
   if startsWith(name, "/") || contains(name, backslash) ...
         || is_drive || any(parts == ".." | parts == "." ...
         | strlength(parts) == 0)
      error('icemodel:verification:fetchFixtures:unsafeMember', ...
         'Archive member has an unsafe path: %s', name)
   end
end

function validateMembers(files, directories, entries)
   %VALIDATEMEMBERS Require exactly declared files and their ancestor directories.
   expected = sort(reshape(string({entries.path}), [], 1));
   missing = setdiff(expected, files);
   extra = setdiff(files, expected);
   if ~isempty(missing)
      error('icemodel:verification:fetchFixtures:missingMember', ...
         'Archive is missing declared file: %s', missing(1))
   end
   if ~isempty(extra)
      error('icemodel:verification:fetchFixtures:extraMember', ...
         'Archive contains undeclared file: %s', extra(1))
   end
   allowed_dirs = ancestorDirectories(expected);
   missing_dirs = setdiff(allowed_dirs, directories);
   if ~isempty(missing_dirs)
      error('icemodel:verification:fetchFixtures:missingMember', ...
         'Archive is missing declared directory: %s', missing_dirs(1))
   end
   extra_dirs = setdiff(directories, allowed_dirs);
   if ~isempty(extra_dirs)
      error('icemodel:verification:fetchFixtures:extraMember', ...
         'Archive contains undeclared directory: %s', extra_dirs(1))
   end
end

function directories = ancestorDirectories(files)
   %ANCESTORDIRECTORIES Derive the only directory members an archive may carry.
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

function files = inventoryFiles(root)
   %INVENTORYFILES Return all extracted regular-file paths relative to root.
   hits = dir(char(fullfile(root, "**", "*")));
   hits = hits(~[hits.isdir]);
   files = strings(numel(hits), 1);
   for k = 1:numel(hits)
      pathname = fullfile(string(hits(k).folder), string(hits(k).name));
      files(k) = relativePosix(root, pathname);
   end
   files = sort(files);
end

%% Promotion and rollback
function promoteCapabilities(extracted_roots, root, selections, ...
      full_selection, work_dir)
   %PROMOTECAPABILITIES Atomically promote every changed capability in one call.
   assertNoDestinationSymlinks(root, selections);
   backup_root = fullfile(work_dir, "backup");
   mkdir(backup_root)

   % Every journal has a declared-file upper bound, so failed transactions can
   % roll back without dynamically growing bookkeeping arrays.
   max_files = numel(full_selection.files);
   backups = repmat( ...
      struct('destination', "", 'backup', ""), max_files, 1);
   promoted = strings(max_files, 1);
   relative_paths = reshape(string({full_selection.files.path}), [], 1);
   destination_paths = string(root) + filesep ...
      + replace(relative_paths, "/", filesep);
   max_created_dirs = sum(count(destination_paths, filesep) + 1);
   created_dirs = strings(max_created_dirs, 1);
   backup_count = 0;
   promoted_count = 0;
   created_count = 0;
   try
      for n = 1:numel(selections)
         selection = selections{n};
         for k = 1:numel(selection.files)
            relpath = string(selection.files(k).path);
            source = underRoot(extracted_roots(n), relpath);
            destination = underRoot(root, relpath);
            backup = underRoot(backup_root, relpath);

            % A declared file may replace only a prior file. Moving a directory
            % would discard undeclared children after a successful transaction.
            if isfolder(destination)
               error('icemodel:verification:fetchFixtures:destinationTypeConflict', ...
                  'Declared file destination is an existing directory: %s', ...
                  relpath)
            elseif isfile(destination)
               icemodel.helpers.ensureDirExists(string(fileparts(backup)));
               moveChecked(destination, backup);
               backup_count = backup_count + 1;
               backups(backup_count) = struct( ...
                  'destination', destination, 'backup', backup);
            end

            % Track every new destination directory so rollback restores even
            % the absence of the pre-call root and intermediate directories.
            [created_dirs, created_count] = ensureParentTracked( ...
               string(fileparts(destination)), created_dirs, created_count);
            moveChecked(source, destination);
            promoted_count = promoted_count + 1;
            promoted(promoted_count) = destination;
         end
      end

      % Verify the complete requested set before backups become disposable.
      [missing, mismatched] = verifyInstalled(root, full_selection);
      if ~isempty(missing) || ~isempty(mismatched)
         error('icemodel:verification:fetchFixtures:promotionVerificationFailed', ...
            'Promoted capabilities failed final verification.')
      end
   catch cause
      try
         rollbackPromotion(promoted(1:promoted_count), ...
            backups(1:backup_count), created_dirs(1:created_count));
      catch rollback_cause
         exception = MException( ...
            'icemodel:verification:fetchFixtures:rollbackFailed', ...
            'Promotion failed and rollback could not restore prior data.');
         exception = addCause(exception, cause);
         exception = addCause(exception, rollback_cause);
         throw(exception)
      end
      rethrow(cause)
   end
end

function assertNoDestinationSymlinks(root, selections)
   %ASSERTNODESTINATIONSYMLINKS Preflight before verification or promotion.
   lexical_root = string(root);
   ancestor_link = icemodel.verification.setup.fixtureCallerSymlink( ...
      lexical_root);
   if strlength(ancestor_link) > 0
      error('icemodel:verification:fetchFixtures:symlinkDestination', ...
         'Fixture destination traverses a symbolic link: %s', ancestor_link)
   end
   assertNotDestinationSymlink(lexical_root, ".");
   root = icemodel.verification.setup.fixtureCanonicalRoot(lexical_root);

   % Check every existing or prospective manifest component inside the root.
   for n = 1:numel(selections)
      selection = selections{n};
      for k = 1:numel(selection.files)
         relpath = string(selection.files(k).path);
         cursor = string(root);
         parts = split(relpath, "/");
         for j = 1:numel(parts)
            cursor = fullfile(cursor, parts(j));
            assertNotDestinationSymlink(cursor, relpath);
         end
      end
   end
end

function assertNotDestinationSymlink(pathname, relpath)
   %ASSERTNOTDESTINATIONSYMLINK Prevent writes through existing symbolic links.
   path_object = java.nio.file.Paths.get(char(pathname), ...
      javaArray('java.lang.String', 0));
   if java.nio.file.Files.isSymbolicLink(path_object)
      error('icemodel:verification:fetchFixtures:symlinkDestination', ...
         'Destination path traverses a symbolic link: %s', relpath)
   end
end

function [created_dirs, created_count] = ensureParentTracked( ...
      parent, created_dirs, created_count)
   %ENSUREPARENTTRACKED Create a parent chain and journal its prior absence.
   missing = strings(count(parent, filesep) + 1, 1);
   missing_count = 0;
   cursor = parent;
   while ~isfolder(cursor)
      if isfile(cursor)
         error('icemodel:verification:fetchFixtures:promotionFailed', ...
            'Destination parent is a file: %s', cursor)
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
         error('icemodel:verification:fetchFixtures:promotionFailed', ...
            'Cannot create %s (%s: %s).', ...
            missing(k), message_id, message)
      end
      created_count = created_count + 1;
      created_dirs(created_count) = missing(k);
   end
end

function rollbackPromotion(promoted, backups, created_dirs)
   %ROLLBACKPROMOTION Restore files and directories to exact pre-call state.
   for k = numel(promoted):-1:1
      if isfile(promoted(k))
         delete(promoted(k))
      elseif isfolder(promoted(k))
         rmdir(promoted(k), 's')
      end
   end
   for k = numel(backups):-1:1
      icemodel.helpers.ensureDirExists( ...
         string(fileparts(backups(k).destination)));
      moveChecked(backups(k).backup, backups(k).destination);
   end
   for k = numel(created_dirs):-1:1
      if isfolder(created_dirs(k))
         contents = dir(char(created_dirs(k)));
         names = string({contents.name});
         contents = contents(~ismember(names, [".", ".."]));
         if ~isempty(contents)
            error('icemodel:verification:fetchFixtures:rollbackFailed', ...
               'New destination directory is not empty: %s', created_dirs(k))
         end
         [ok, message, message_id] = rmdir(created_dirs(k));
         if ~ok
            error('icemodel:verification:fetchFixtures:rollbackFailed', ...
               'Cannot remove %s (%s: %s).', ...
               created_dirs(k), message_id, message)
         end
      end
   end
end

function moveChecked(source, destination)
   %MOVECHECKED Move one path and preserve the platform error on failure.
   [ok, message, message_id] = movefile(source, destination, 'f');
   if ~ok
      error('icemodel:verification:fetchFixtures:promotionFailed', ...
         'Cannot move %s to %s (%s: %s).', ...
         source, destination, message_id, message)
   end
end

%% Path and error helpers
function pathname = underRoot(root, relpath)
   %UNDERROOT Convert a validated relative POSIX path to a local path.
   parts = cellstr(split(relpath, "/"));
   pathname = string(fullfile(root, parts{:}));
end

function relpath = relativePosix(root, pathname)
   %RELATIVEPOSIX Strip a known root prefix and normalize path separators.
   root = icemodel.verification.setup.fixtureCanonicalRoot(string(root));
   pathname = icemodel.verification.setup.fixtureCanonicalRoot(string(pathname));
   relpath = icemodel.verification.setup.fixtureRelativePosix( ...
      root, pathname);
end

function removeTree(pathname)
   %REMOVETREE Delete transaction storage after success or failure.
   if isfolder(pathname)
      rmdir(pathname, 's')
   end
end

function throwWithCause(identifier, cause, message, varargin)
   %THROWWITHCAUSE Throw a stable public error while retaining root cause.
   exception = MException(identifier, message, varargin{:});
   exception = addCause(exception, cause);
   throw(exception)
end
