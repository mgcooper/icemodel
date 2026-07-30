function [source_dir, status] = fetchLaughTests(kwargs)
   %FETCHLAUGHTESTS Locate or verify the Laugh-Tests source checkout.
   %
   %  source_dir = icemodel.verification.setup.fetchLaughTests()
   %  source_dir = icemodel.verification.setup.fetchLaughTests( ...
   %     cache_dir="/some/other/path", strict=true)
   %
   %  Resolves the local Laugh-Tests source-bundle directory used by
   %  icemodel.verification.setup.importLaughTests to stage the
   %  Colbeck 1976 smoke verification artifacts.
   %
   %  Default cache directory:
   %    data/verification/laugh_tests/
   %
   %  This is a gitignored cache. Populate it by either:
   %    1. cloning the upstream repository into the cache dir, e.g.
   %         git clone https://github.com/KyleKlenk/Laugh-Tests \
   %              data/verification/laugh_tests
   %    2. or pointing kwargs.cache_dir at an existing local checkout.
   %       If a sibling Laugh-Tests checkout exists at the project
   %       root parent directory, it is auto-detected (see below).
   %
   %  Files required by importLaughTests (Colbeck case):
   %    test_cases/input_data/colbeck1976/colbeck1976_forcing.nc
   %    validation_data/m2_mac_Sept23/colbeck1976/
   %        colbeck1976-exp1_G1-1_timestep.nc
   %    validation_data/m2_mac_Sept23/colbeck1976/
   %        colbeck1976-exp2_G1-1_timestep.nc
   %    validation_data/m2_mac_Sept23/colbeck1976/
   %        colbeck1976-exp3_G1-1_timestep.nc
   %
   %  Behaviour
   %    - If cache_dir contains a Laugh-Tests checkout with all
   %      required files, returns the directory path so the caller
   %      can pass it to importLaughTests(source_dir).
   %    - If cache_dir is empty or missing files, prints actionable
   %      retrieval instructions (clone URL, expected paths) and
   %      either errors (strict=true, default) or returns the
   %      partial path with a warning (strict=false).
   %    - As a convenience, if the default cache_dir is empty AND a
   %      sibling ../Laugh-Tests checkout exists at the icemodel
   %      project root parent directory (e.g. ~/MATLAB/projects/Laugh-Tests
   %      when icemodel lives at ~/MATLAB/projects/icemodel), it is
   %      auto-detected and returned. The user may set kwargs.cache_dir
   %      explicitly to disable this fallback.
   %
   %  Role
   %    Validator. With create_cache_dir=true, the fetch helper creates the
   %    cache directory before checking the required Laugh-Tests files. With
   %    false, it reports the same status without mutation. A successful strict
   %    result lets importLaughTests assume the checkout is complete.
   %
   %  Name-value
   %    cache_dir : string (default data/verification/laugh_tests)
   %        Local Laugh-Tests checkout directory.
   %    strict : logical (default true)
   %        Error when the checkout is missing or incomplete.
   %    silent : logical (default false)
   %        Suppress the retrieval-instructions printout.
   %    create_cache_dir : logical (default true)
   %        Create the resolved cache directory before validation.
   %
   %  Returns
   %    source_dir : string
   %        Absolute path to the resolved Laugh-Tests checkout.
   %    status : struct array
   %        Shared fetch-status rows describing missing checkout items.
   %
   % See also: icemodel.verification.setup.importLaughTests,
   %  icemodel.verification.setup.fetchEsmSnowmip

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.strict   (1, 1) logical = true
      kwargs.silent   (1, 1) logical = false
      kwargs.create_cache_dir (1, 1) logical = true
   end

   cache_dir = icemodel.verification.setup.resolveFetchCacheDir( ...
      kwargs.cache_dir, defaultCacheDir());
   user_supplied = ~strcmp(cache_dir, defaultCacheDir());

   % Cache creation is explicit so dry-run callers can remain non-mutating.
   if kwargs.create_cache_dir
      icemodel.helpers.ensureDirExists(cache_dir);
   end

   % Required files (Colbeck case).
   required = colbeckRequiredFiles();

   % Probe the requested cache for completeness first.
   [ok, missing] = checkLaughTestsCheckout(cache_dir, required);

   % Sibling-checkout fallback: if the default cache is incomplete
   % AND the user did not pass cache_dir explicitly AND a sibling
   % ../Laugh-Tests directory has a complete checkout, prefer that.
   % This preserves the sibling-checkout workflow without requiring a
   % data/verification migration today.
   if ~ok && ~user_supplied
      data_root = icemodel.internal.fullpath('data');
      repo_root = fileparts(data_root);
      sibling = string(fullfile(repo_root, '..', 'Laugh-Tests'));
      if isfolder(sibling)
         [ok_sibling, missing_sibling] = checkLaughTestsCheckout( ...
            sibling, required);
         if ok_sibling
            status = icemodel.verification.setup.fetchMissingStatus( ...
               strings(0, 1));
            source_dir = sibling;
            return
         else
            % Sibling exists but is incomplete; fall through with
            % the sibling path so the retrieval banner names it.
            cache_dir = sibling;
            missing = missing_sibling;
         end
      end
   end

   % Build cache status before strict handling.
   if ok
      missing = strings(0, 1);
   end
   status = icemodel.verification.setup.fetchMissingStatus(missing);
   source_dir = icemodel.verification.setup.finishFetchStatus( ...
      cache_dir, status, strict=kwargs.strict, silent=kwargs.silent, ...
      error_id="icemodel:verification:fetchLaughTests:missingSources", ...
      error_label="Laugh-Tests", ...
      banner_callback=@(~, ~, items) printRetrievalBanner(cache_dir, items), ...
      error_callback=@(~, ~, items) throwMissingSources(cache_dir, items));
end

%% Local helpers
function printRetrievalBanner(cache_dir, missing)
   %PRINTRETRIEVALBANNER Print Laugh-Tests retrieval instructions.
   fprintf('\n');
   fprintf('=== Laugh-Tests source checkout incomplete ===\n');
   fprintf('Cache directory: %s\n', cache_dir);
   fprintf('Missing files (relative to cache directory):\n');
   for j = 1:numel(missing)
      fprintf('  - %s\n', missing(j));
   end
   fprintf('\nRetrieval:\n');
   fprintf('  Upstream: https://github.com/KyleKlenk/Laugh-Tests\n');
   fprintf('  Bundle:   m2_mac_Sept23 validation set\n');
   fprintf('  Manual workflow: clone the repository into the cache dir.\n');
   fprintf('    git clone https://github.com/KyleKlenk/Laugh-Tests %s\n', ...
      cache_dir);
   fprintf('\nAfter retrieval, re-run:\n');
   fprintf('  icemodel.verification.setup.fetchLaughTests()\n');
   fprintf('  icemodel.verification.setup.importLaughTests(source_dir, overwrite=true)\n');
   fprintf('\n');
end

function throwMissingSources(cache_dir, ~)
   %THROWMISSINGSOURCES Raise the Laugh-Tests missing-cache error.
   error('icemodel:verification:fetchLaughTests:missingSources', ...
      ['Laugh-Tests source checkout incomplete in %s. ' ...
      'See retrieval instructions above.'], cache_dir);
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical Laugh-Tests source-cache directory.
   % Repo-root developer resource (gitignored). Resolve with
   % icemodel.internal.fullpath('data') (always the repo-local data root), not
   % icemodel.getpath('data'): getpath returns the ICEMODEL_DATA_PATH env var,
   % which a scoped demo or test case sets to its owned data tree, so it points
   % away from the repo root under the active suite config.
   pathname = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'laugh_tests'));
end

function required = colbeckRequiredFiles()
   %COLBECKREQUIREDFILES Files needed by importLaughTests for colbeck1976.
   required = string([ ...
      "test_cases/input_data/colbeck1976/colbeck1976_forcing.nc"; ...
      "validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp1_G1-1_timestep.nc"; ...
      "validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp2_G1-1_timestep.nc"; ...
      "validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp3_G1-1_timestep.nc"]);
end

function [ok, missing] = checkLaughTestsCheckout(cache_dir, required)
   %CHECKLAUGHTESTSCHECKOUT Per-file presence check for a Laugh-Tests dir.

   if ~isfolder(cache_dir)
      missing = required;
      ok = false;
      return
   end

   exists = arrayfun(@(p) isfile(fullfile(cache_dir, p)), required);
   missing = required(~exists);
   ok = isempty(missing);
end
