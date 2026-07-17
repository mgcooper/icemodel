function result = fetchFixtures(version, kwargs)
   %FETCHFIXTURES Verify or provision the demo fixture data from a release asset.
   %
   %  result = icemodel.verification.setup.fetchFixtures()
   %  result = icemodel.verification.setup.fetchFixtures("v0.1.0")
   %  result = icemodel.verification.setup.fetchFixtures("v0.1.0", ...
   %     archive="/path/icemodel-fixtures-v0.1.0.tar.gz", extract=true)
   %  result = icemodel.verification.setup.fetchFixtures(manifest="/path/MANIFEST.json")
   %
   %  Consumer side of the release-with-assets fixture strategy and the
   %  on-demand fixture provisioner. Mirrors the verify-cache contract of
   %  icemodel.verification.setup.fetchSumup / fetchEsmSnowmip: it verifies a
   %  local cache (here, the demo fixture data tree) against a bundle manifest
   %  and, when fixtures are missing, prints actionable instructions to download
   %  the GitHub release asset rather than auto-downloading.
   %
   %  Three modes, in priority order:
   %    1. archive given        : extract the local archive into the demo data
   %                              root, then verify every file against its
   %                              manifest SHA-256.
   %    2. manifest resolvable  : verify the on-disk fixtures (committed today,
   %                              extracted after a flip) against the manifest
   %                              SHA-256.
   %    3. no manifest          : presence-only check of the fixture file set.
   %                              This keeps the helper a non-breaking no-op when
   %                              the committed fixtures are already on disk and
   %                              no bundle manifest exists yet (the current
   %                              repo state), so wiring it into the test
   %                              bootstrap cannot break offline CI.
   %
   %  Input
   %    version : string (optional)
   %        Release version tag. Only used to compose the download-banner asset
   %        URL and to locate a staged manifest by name; verification itself is
   %        manifest-driven, not version-driven.
   %
   %  Name-value
   %    root : string (default demo/data root)
   %        Demo data root to verify / extract into.
   %    manifest : string (default "")
   %        Bundle MANIFEST.json to verify against. When empty, fetchFixtures
   %        looks for a staged manifest matching `version` under staging_dir;
   %        if none is found it falls back to a presence-only check.
   %    archive : string (default "")
   %        Local fixture archive (.tar.gz) to extract before verifying. The
   %        archive's sibling *.MANIFEST.json is used when `manifest` is unset.
   %    staging_dir : string (default <repo>/release-staging)
   %        Where to look for a staged manifest when `manifest` is unset.
   %    extract : logical (default true when archive given, else false)
   %        Extract `archive` into root before verifying.
   %    repo : string (default "mgcooper/icemodel")
   %        GitHub owner/repo for the download-instruction banner.
   %    strict : logical (default true)
   %        Error when fixtures are missing or fail SHA-256 verification.
   %    silent : logical (default false)
   %        Suppress the retrieval-instructions banner.
   %
   %  Returns
   %    result : struct with fields
   %        ok            : true when the fixture set verified (or is present)
   %        mode          : "verified" | "present" | "missing"
   %        root          : the demo data root checked
   %        manifest_file : manifest used ("" when presence-only)
   %        missing       : string vector of missing fixture paths
   %        mismatched    : string vector of files whose SHA-256 differs
   %
   % See also: icemodel.verification.setup.packFixtures,
   %  icemodel.verification.setup.fixtureFileList,
   %  icemodel.verification.setup.fetchSumup

   arguments
      version (1, 1) string = ""
      kwargs.root (1, 1) string = icemodel.internal.fullpath('demo', 'data')
      kwargs.manifest (1, 1) string = ""
      kwargs.archive (1, 1) string = ""
      kwargs.staging_dir (1, 1) string = icemodel.internal.fullpath('release-staging')
      kwargs.extract (1, 1) logical = true
      kwargs.repo (1, 1) string = "mgcooper/icemodel"
      kwargs.strict (1, 1) logical = true
      kwargs.silent (1, 1) logical = false
   end

   root = kwargs.root;

   % Ensure the demo data root exists so a download banner points at a path that
   % is already there.
   icemodel.helpers.ensureDirExists(root);

   % Mode 1: extract a local archive only when requested. With extract=false the
   % helper verifies the archive manifest against files already on disk.
   if strlength(kwargs.archive) > 0 && kwargs.extract
      if ~isfile(kwargs.archive)
         error('icemodel:verification:fetchFixtures:noArchive', ...
            'Archive not found: %s', kwargs.archive);
      end
      % untar restores the relative paths recorded at pack time over root.
      untar(char(kwargs.archive), char(root));
   end

   % Resolve the manifest to verify against. Priority: explicit manifest kwarg,
   % then the archive's sibling manifest, then a version-named staged manifest.
   manifest_file = resolveManifest(kwargs.manifest, kwargs.archive, ...
      kwargs.staging_dir, version);

   % Without a manifest, fall back to a presence-only check. This is the
   % non-breaking path used by the test bootstrap today: committed fixtures are
   % present, no bundle exists yet, so this is a clean no-op.
   if strlength(manifest_file) == 0
      result = presenceOnly(root, kwargs.strict, kwargs.silent, ...
         version, kwargs.repo);
      return
   end

   % Manifest-driven verification: re-hash each listed file on disk and compare
   % to the recorded SHA-256.
   manifest = jsondecode(fileread(manifest_file));
   [missing, mismatched] = verifyAgainstManifest(root, manifest);

   ok = isempty(missing) && isempty(mismatched);
   if ok
      result = struct('ok', true, 'mode', "verified", 'root', root, ...
         'manifest_file', string(manifest_file), ...
         'missing', strings(0, 1), 'mismatched', strings(0, 1));
      return
   end

   % Print the actionable download-and-extract banner when verification fails.
   if ~kwargs.silent
      printBanner(root, string(manifest_file), missing, mismatched, ...
         version, kwargs.repo);
   end

   if kwargs.strict
      error('icemodel:verification:fetchFixtures:incompleteFixtures', ...
         ['Demo fixture data under %s does not match the bundle manifest. ' ...
         'Missing: %d, mismatched: %d. See instructions above.'], ...
         root, numel(missing), numel(mismatched));
   end

   result = struct('ok', false, 'mode', "missing", 'root', root, ...
      'manifest_file', string(manifest_file), ...
      'missing', missing, 'mismatched', mismatched);
end

%% Local helpers
function manifest_file = resolveManifest(manifest, archive, staging_dir, version)
   %RESOLVEMANIFEST Pick the manifest path from explicit / sibling / staged.

   % Explicit manifest wins.
   if strlength(manifest) > 0
      if ~isfile(manifest)
         error('icemodel:verification:fetchFixtures:noManifest', ...
            'Manifest not found: %s', manifest);
      end
      manifest_file = manifest;
      return
   end

   % An archive carries a sibling *.MANIFEST.json (packFixtures writes both).
   if strlength(archive) > 0
      sibling = replace(string(archive), ".tar.gz", ".MANIFEST.json");
      if isfile(sibling)
         manifest_file = sibling;
         return
      end
   end

   % A version-named staged manifest under staging_dir (post-pack, pre-flip).
   if strlength(version) > 0
      staged = fullfile(staging_dir, ...
         "icemodel-fixtures-" + version + ".MANIFEST.json");
      if isfile(staged)
         manifest_file = staged;
         return
      end
   end

   % None found: signal presence-only fallback.
   manifest_file = "";
end

function [missing, mismatched] = verifyAgainstManifest(root, manifest)
   %VERIFYAGAINSTMANIFEST Re-hash each manifest file on disk and compare.

   entries = manifest.files;
   n = numel(entries);
   missing = strings(n, 1);
   mismatched = strings(n, 1);
   n_missing = 0;
   n_mismatch = 0;

   for k = 1:n
      relpath = string(entries(k).path);
      abspath = fullfile(root, relpath);
      if ~isfile(abspath)
         n_missing = n_missing + 1;
         missing(n_missing) = relpath;
         continue
      end
      % SHA-256 mismatch means a stale or corrupt fixture, not just absence.
      digest = icemodel.verification.setup.fileSha256(abspath);
      if ~strcmp(digest, string(entries(k).sha256))
         n_mismatch = n_mismatch + 1;
         mismatched(n_mismatch) = relpath;
      end
   end

   missing = missing(1:n_missing);
   mismatched = mismatched(1:n_mismatch);
end

function result = presenceOnly(root, strict, silent, version, repo)
   %PRESENCEONLY Presence-only check when no bundle manifest is available.
   %
   % Used as the non-breaking fallback: with the committed fixtures present and
   % no manifest yet, this returns ok with mode "present" (a no-op). It only
   % becomes actionable once the fixtures have been flipped to a release asset
   % and the tree is empty.

   files = icemodel.verification.setup.fixtureFileList(root=root);

   % An empty fixture set means neither committed fixtures nor an extracted
   % asset are present: surface the download instructions.
   if isempty(files)
      if ~silent
         printBanner(root, "", strings(0, 1), strings(0, 1), version, repo);
      end
      if strict
         error('icemodel:verification:fetchFixtures:noFixtures', ...
            ['No demo fixture data found under %s and no bundle manifest ' ...
            'available. See instructions above.'], root);
      end
      result = struct('ok', false, 'mode', "missing", 'root', root, ...
         'manifest_file', "", 'missing', strings(0, 1), ...
         'mismatched', strings(0, 1));
      return
   end

   % Fixtures present (committed, as today): clean no-op.
   result = struct('ok', true, 'mode', "present", 'root', root, ...
      'manifest_file', "", 'missing', strings(0, 1), ...
      'mismatched', strings(0, 1));
end

function printBanner(root, manifest_file, missing, mismatched, version, repo)
   %PRINTBANNER Stable, grep-able download-and-extract instructions.
   %
   % Mirrors the fetchSumup / fetchEsmSnowmip retrieval-banner contract: name
   % the asset URL pattern and the exact extraction step; never auto-download.

   tag = version;
   if strlength(tag) == 0
      tag = "<version>";
   end

   fprintf('\n');
   fprintf('=== icemodel demo fixture data incomplete ===\n');
   fprintf('Demo data root: %s\n', root);
   if strlength(manifest_file) > 0
      fprintf('Bundle manifest: %s\n', manifest_file);
   end
   if ~isempty(missing)
      fprintf('Missing fixture files (%d):\n', numel(missing));
      for j = 1:min(numel(missing), 10)
         fprintf('  - %s\n', missing(j));
      end
      if numel(missing) > 10
         fprintf('  ... and %d more\n', numel(missing) - 10);
      end
   end
   if ~isempty(mismatched)
      fprintf('SHA-256 mismatched fixture files (%d):\n', numel(mismatched));
      for j = 1:min(numel(mismatched), 10)
         fprintf('  - %s\n', mismatched(j));
      end
   end
   fprintf('\nRetrieval (GitHub release asset):\n');
   fprintf('  Release:  https://github.com/%s/releases/tag/fixtures-%s\n', ...
      repo, tag);
   fprintf('  Asset:    icemodel-fixtures-%s.tar.gz\n', tag);
   fprintf('  Download (gh CLI):\n');
   fprintf('    gh release download fixtures-%s --repo %s \\\n', tag, repo);
   fprintf('      --pattern "icemodel-fixtures-%s.tar.gz"\n', tag);
   fprintf('  Then extract into the demo data root and verify:\n');
   fprintf('    icemodel.verification.setup.fetchFixtures("%s", ...\n', tag);
   fprintf('      archive="icemodel-fixtures-%s.tar.gz", extract=true)\n', tag);
   fprintf('\n');
end
