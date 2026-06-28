function result = packFixtures(version, kwargs)
   %PACKFIXTURES Bundle committed demo fixture data into a versioned archive.
   %
   %  result = icemodel.verification.setup.packFixtures(version)
   %  result = icemodel.verification.setup.packFixtures("v0.1.0", ...
   %     staging_dir="/tmp/stage", root="/some/demo/data")
   %
   %  Producer side of the release-with-assets fixture strategy. Collects the
   %  heavy committed binary fixture DATA (demo/data/eval/**/*.mat plus the demo
   %  met/userdata forcing the unit suite loads, enumerated by
   %  icemodel.verification.setup.fixtureFileList) into a single compressed
   %  tar.gz archive plus a sidecar MANIFEST.json holding the file list, each
   %  file's SHA-256, and the version. The archive is verifiable: a consumer
   %  re-hashes the extracted files and compares against the manifest.
   %
   %  Both artifacts are written to a gitignored staging directory
   %  (release-staging/ by default) so the bundle is NEVER committed. The user
   %  then uploads them as a GitHub release asset (see README_FIXTURES.md); the
   %  consumer side is icemodel.verification.setup.fetchFixtures.
   %
   %  This pass does NOT delete or move the committed fixtures: packFixtures only
   %  reads them. The one-time flip (git rm the committed .mat, gitignore them,
   %  upload the asset) is the user's step, tracked as icemodel-1ps.17.
   %
   %  Input
   %    version : string
   %        Bundle version tag (e.g. "v0.1.0"). Stamped into the manifest and
   %        the artifact filenames so an asset is uniquely identifiable.
   %
   %  Name-value
   %    root : string (default demo/data root)
   %        Demo data root to bundle. Defaults to the canonical committed tree;
   %        tests pass a temporary root.
   %    staging_dir : string (default <repo>/release-staging)
   %        Gitignored output directory for the archive + manifest.
   %    overwrite : logical (default false)
   %        Overwrite existing artifacts of the same version instead of erroring.
   %    silent : logical (default false)
   %        Suppress the size-saving report printout.
   %
   %  Returns
   %    result : struct with fields
   %        version          : the bundle version tag
   %        archive_file     : absolute path to the .tar.gz
   %        manifest_file    : absolute path to the MANIFEST.json
   %        manifest         : the manifest struct (version, files, created)
   %        num_files        : number of fixture files bundled
   %        source_bytes     : total committed footprint of the fixtures
   %        archive_bytes    : compressed archive size
   %        saving_bytes     : source_bytes - archive_bytes
   %
   % See also: icemodel.verification.setup.fetchFixtures,
   %  icemodel.verification.setup.fixtureFileList,
   %  icemodel.verification.setup.fileSha256

   arguments
      version (1, 1) string
      kwargs.root (1, 1) string = string(icemodel.internal.fullpath('demo', 'data'))
      kwargs.staging_dir (1, 1) string = defaultStagingDir()
      kwargs.overwrite (1, 1) logical = false
      kwargs.silent (1, 1) logical = false
   end

   root = kwargs.root;

   % Enumerate the fixture set from the single source of truth so the bundle and
   % the later verification cannot drift.
   files = icemodel.verification.setup.fixtureFileList(root=root);
   if isempty(files)
      error('icemodel:verification:packFixtures:noFixtures', ...
         'No fixture data found under %s. Nothing to pack.', root);
   end

   % Ensure the gitignored staging directory exists.
   icemodel.helpers.ensureDirExists(kwargs.staging_dir);

   % Build the versioned artifact paths. The manifest filename mirrors the
   % archive name so an uploaded release carries an obvious matched pair.
   base = "icemodel-fixtures-" + version;
   archive_file = fullfile(kwargs.staging_dir, base + ".tar.gz");
   manifest_file = fullfile(kwargs.staging_dir, base + ".MANIFEST.json");

   % Refuse to clobber an existing bundle unless asked, so a prior release asset
   % staged on disk is not silently overwritten.
   if ~kwargs.overwrite && (isfile(archive_file) || isfile(manifest_file))
      error('icemodel:verification:packFixtures:exists', ...
         ['Bundle artifacts for %s already exist in %s. ' ...
         'Pass overwrite=true to replace them.'], version, kwargs.staging_dir);
   end

   % Hash every fixture file and total the committed footprint.
   n = numel(files);
   entries = repmat(struct('path', "", 'sha256', "", 'bytes', 0), n, 1);
   source_bytes = 0;
   for k = 1:n
      abspath = fullfile(root, files(k));
      info = dir(char(abspath));
      entries(k).path = files(k);
      entries(k).sha256 = icemodel.verification.setup.fileSha256(abspath);
      entries(k).bytes = info.bytes;
      source_bytes = source_bytes + info.bytes;
   end

   % Assemble the verifiable manifest (version + per-file path/sha256/bytes +
   % a creation timestamp for provenance).
   manifest = struct();
   manifest.bundle = base;
   manifest.version = version;
   manifest.created = string(datetime("now", TimeZone="UTC", ...
      Format="uuuu-MM-dd'T'HH:mm:ss'Z'"));
   manifest.num_files = n;
   manifest.source_bytes = source_bytes;
   manifest.files = entries;

   % Write the compressed archive of the fixture files, with paths relative to
   % root so extraction over demo/data restores the committed layout exactly.
   tar(char(archive_file), cellstr(files), char(root));

   % Write the sidecar manifest as pretty JSON for a reviewable asset.
   writeJson(manifest_file, manifest);

   % Report sizes so the saving (committed footprint vs compressed asset) is
   % visible. archive_bytes is read back after writing the real file.
   archive_info = dir(char(archive_file));
   archive_bytes = archive_info.bytes;

   result = struct( ...
      'version', version, ...
      'archive_file', string(archive_file), ...
      'manifest_file', string(manifest_file), ...
      'manifest', manifest, ...
      'num_files', n, ...
      'source_bytes', source_bytes, ...
      'archive_bytes', archive_bytes, ...
      'saving_bytes', source_bytes - archive_bytes);

   if ~kwargs.silent
      reportSaving(result);
   end
end

%% Local helpers
function pathname = defaultStagingDir()
   %DEFAULTSTAGINGDIR Canonical gitignored release-staging directory.
   %
   % Lives at <repo>/release-staging so the .gitignore rule that ignores it is
   % stable and the bundle never enters version control.
   pathname = string(icemodel.internal.fullpath('release-staging'));
end

function writeJson(pathname, data)
   %WRITEJSON Write a struct as pretty JSON.
   %
   % Mirrors icemodel.verification.setup.writeManifest so the staged manifest
   % reads the same as the committed family manifests.
   fid = fopen(pathname, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, jsonencode(data, PrettyPrint=true), 'char');
end

function reportSaving(result)
   %REPORTSAVING Print the committed-footprint vs compressed-asset comparison.
   mb = @(b) double(b) / 1024 / 1024;
   fprintf('\n');
   fprintf('=== icemodel fixture bundle packed ===\n');
   fprintf('Version:          %s\n', result.version);
   fprintf('Files bundled:    %d\n', result.num_files);
   fprintf('Archive:          %s\n', result.archive_file);
   fprintf('Manifest:         %s\n', result.manifest_file);
   fprintf('Committed source: %.1f MB\n', mb(result.source_bytes));
   fprintf('Compressed asset: %.1f MB\n', mb(result.archive_bytes));
   fprintf('Saving:           %.1f MB (%.0f%%)\n', ...
      mb(result.saving_bytes), ...
      100 * double(result.saving_bytes) / double(result.source_bytes));
   fprintf('\n');
end
