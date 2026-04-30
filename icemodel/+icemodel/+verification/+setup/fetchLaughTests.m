function source_dir = fetchLaughTests(kwargs)
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
   %    data/verification/snow/laugh_tests/
   %
   %  This is a gitignored cache. Populate it by either:
   %    1. cloning the upstream repository into the cache dir, e.g.
   %         git clone https://github.com/KyleKlenk/Laugh-Tests \
   %              data/verification/snow/laugh_tests
   %    2. or pointing kwargs.cache_dir at an existing local checkout
   %       (Matt's machine has a co-located ../Laugh-Tests checkout
   %       used by the original development pass).
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
   %    - As a convenience for Matt's existing workflow, if the
   %      default cache_dir is empty AND a sibling ../Laugh-Tests
   %      checkout exists at the icemodel project root parent
   %      directory (e.g. ~/MATLAB/projects/Laugh-Tests when
   %      icemodel lives at ~/MATLAB/projects/icemodel), it is
   %      auto-detected and returned. The user may set kwargs.cache_dir
   %      explicitly to disable this fallback.
   %
   %  Name-value
   %    cache_dir : string (default data/verification/snow/laugh_tests)
   %        Local Laugh-Tests checkout directory.
   %    strict : logical (default true)
   %        Error when the checkout is missing or incomplete.
   %    silent : logical (default false)
   %        Suppress the retrieval-instructions printout.
   %
   %  Returns
   %    source_dir : string
   %        Absolute path to the resolved Laugh-Tests checkout.
   %
   % See also: icemodel.verification.setup.importLaughTests,
   %  icemodel.verification.setup.fetchEsmSnowmip

   arguments
      kwargs.cache_dir (1, 1) string = ""
      kwargs.strict   (1, 1) logical = true
      kwargs.silent   (1, 1) logical = false
   end

   % Resolve the default cache directory under the canonical data
   % root (<repo>/data/) returned by icemodel.getpath('data').
   data_root = icemodel.getpath('data');
   default_cache = string(fullfile(data_root, 'verification', ...
      'snow', 'laugh_tests'));

   user_supplied = kwargs.cache_dir ~= "";
   if user_supplied
      cache_dir = kwargs.cache_dir;
   else
      cache_dir = default_cache;
   end

   % Required files (Colbeck case).
   required = string([ ...
      "test_cases/input_data/colbeck1976/colbeck1976_forcing.nc"; ...
      "validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp1_G1-1_timestep.nc"; ...
      "validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp2_G1-1_timestep.nc"; ...
      "validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp3_G1-1_timestep.nc"]);

   % Probe the requested cache for completeness first.
   [ok, missing] = checkLaughTestsCheckout(cache_dir, required);

   % Sibling-checkout fallback: if the default cache is incomplete
   % AND the user did not pass cache_dir explicitly AND a sibling
   % ../Laugh-Tests directory has a complete checkout, prefer that.
   % This preserves Matt's existing workflow without requiring a
   % data/verification migration today. The repo root is data_root's
   % parent (data_root = <repo>/data).
   if ~ok && ~user_supplied
      repo_root = fileparts(data_root);
      sibling = string(fullfile(repo_root, '..', 'Laugh-Tests'));
      if exist(sibling, 'dir') == 7
         [ok_sibling, missing_sibling] = checkLaughTestsCheckout( ...
            sibling, required);
         if ok_sibling
            source_dir = sibling;
            return
         else
            % Fall through to the default cache message below; the
            % sibling checkout exists but is incomplete, which is
            % surfaced under the same retrieval banner.
            cache_dir = sibling;
            missing = missing_sibling;
         end
      end
   end

   if ok
      source_dir = string(cache_dir);
      return
   end

   % Print actionable retrieval instructions when files are missing.
   if ~kwargs.silent
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

   if kwargs.strict
      error('icemodel:verification:fetchLaughTests:missingSources', ...
         ['Laugh-Tests source checkout incomplete in %s. ' ...
         'See retrieval instructions above.'], ...
         cache_dir);
   end

   source_dir = string(cache_dir);
end

function [ok, missing] = checkLaughTestsCheckout(cache_dir, required)
   %CHECKLAUGHTESTSCHECKOUT Per-file presence check for a Laugh-Tests dir.

   missing = strings(0, 1);
   if exist(cache_dir, 'dir') ~= 7
      missing = required;
      ok = false;
      return
   end
   for i = 1:numel(required)
      pathname = fullfile(cache_dir, required(i));
      if exist(pathname, 'file') ~= 2
         missing(end + 1, 1) = required(i); %#ok<AGROW>
      end
   end
   ok = isempty(missing);
end
