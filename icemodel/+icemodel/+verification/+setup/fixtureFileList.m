function files = fixtureFileList(kwargs)
   %FIXTUREFILELIST Enumerate the committed demo fixture data files.
   %
   %  files = icemodel.verification.setup.fixtureFileList()
   %  files = icemodel.verification.setup.fixtureFileList(root="/some/demo/data")
   %
   %  Single source of truth for the set of binary fixture DATA files that the
   %  unit/verification suites load from the committed demo data tree and that
   %  the release-with-assets workflow bundles into a versioned archive. Both
   %  icemodel.verification.setup.packFixtures (bundle producer) and
   %  icemodel.verification.setup.fetchFixtures (verify/provision) call this so
   %  the packed set and the verified set can never drift apart.
   %
   %  The enumerated set is exactly the heavy committed .mat data:
   %    demo/data/eval/**/*.mat          observation / reference bundles
   %    demo/data/input/met/**/*.mat        forcing timetables
   %    demo/data/input/userdata/**/*.mat   userdata swap files
   %
   %  The lean, reviewable companions (per-family manifest.json, .gitkeep) are
   %  deliberately EXCLUDED: they stay committed even after the fixtures flip to
   %  a release asset, so they must not live inside the bundle.
   %
   %  Name-value
   %    root : string (default demo/data root)
   %        Demo data root to enumerate under. Defaults to the canonical
   %        repo-local demo/data so callers normally pass nothing; tests pass a
   %        temporary root to exercise pack/fetch without touching the repo.
   %
   %  Returns
   %    files : string column vector
   %        Repo-relative POSIX paths (relative to root), sorted, so the
   %        manifest and the on-disk verification are order-stable across
   %        platforms. Empty when the root holds no fixture data.
   %
   % See also: icemodel.verification.setup.packFixtures,
   %  icemodel.verification.setup.fetchFixtures

   arguments
      kwargs.root (1, 1) string = defaultDemoDataRoot()
   end

   root = kwargs.root;

   % The three fixture-data globs, each rooted at a known subtree so the bundle
   % layout mirrors the committed layout 1:1 (extracting the archive over
   % demo/data restores the exact paths).
   patterns = [ ...
      fullfile("eval", "**", "*.mat"); ...
      fullfile("input", "met", "**", "*.mat"); ...
      fullfile("input", "userdata", "**", "*.mat")];

   % Collect absolute matches per pattern, then convert to root-relative POSIX
   % paths so the manifest is portable.
   parts = cell(numel(patterns), 1);
   for k = 1:numel(patterns)
      hits = dir(char(fullfile(root, patterns(k))));
      hits = hits(~[hits.isdir]);
      rel = strings(numel(hits), 1);
      for j = 1:numel(hits)
         rel(j) = toRelativePosix(root, ...
            fullfile(string(hits(j).folder), string(hits(j).name)));
      end
      parts{k} = rel;
   end

   % Concatenate, drop any empties, sort for a deterministic manifest order.
   files = vertcat(parts{:});
   files = files(strlength(files) > 0);
   files = sort(unique(files));
end

%% Local helpers
function pathname = defaultDemoDataRoot()
   %DEFAULTDEMODATAROOT Canonical committed demo data root (<repo>/demo/data).
   %
   % The committed fixtures live under demo/data (the casename="demo"/"test"
   % config points ICEMODEL_INPUT_PATH / ICEMODEL_EVAL_PATH here), so the
   % release-asset workflow bundles and verifies that tree by default.
   pathname = string(icemodel.internal.fullpath('demo', 'data'));
end

function rel = toRelativePosix(root, abspath)
   %TORELATIVEPOSIX Strip the root prefix and normalize separators to '/'.
   %
   % Manifests store forward-slash relative paths so a bundle packed on macOS
   % verifies unchanged on Linux CI.
   root = char(root);
   abspath = char(abspath);
   if ~endsWith(root, filesep)
      root = [root filesep];
   end
   if startsWith(abspath, root)
      rel = string(abspath(numel(root) + 1:end));
   else
      rel = string(abspath);
   end
   rel = replace(rel, "\", "/");
end
