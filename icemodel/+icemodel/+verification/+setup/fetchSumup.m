function source_dir = fetchSumup(kwargs)
   %FETCHSUMUP Locate or verify the SUMup firn source files.
   %
   %  source_dir = icemodel.verification.setup.fetchSumup()
   %  source_dir = icemodel.verification.setup.fetchSumup( ...
   %     cache_dir="/some/other/path", strict=true)
   %
   %  Resolves the local source-cache directory holding the SUMup firn
   %  observation files (density, accumulation/SMB, subsurface temperature)
   %  used by icemodel.verification.setup.importSumup to stage the per-point
   %  firn verification artifacts. By default the cache lives at
   %  data/verification/sumup/ (gitignored) and is populated by the
   %  user / developer following the retrieval instructions printed when
   %  files are missing.
   %
   %  SUMup is the firn observation source. It subsumes Humphrey 2012
   %  subsurface temperature and GreenTRACS accumulation, so no separate
   %  sourcing is required. FirnCover compaction strain is NOT in SUMup and
   %  is out of scope.
   %
   %  Expected files (SUMup 2024 release; CSV or NetCDF accepted):
   %    SUMup_*_density_*.{nc,csv}        firn/snow density profiles
   %    SUMup_*_accumulation_*.{nc,csv}   SMB / accumulation
   %    SUMup_*_temperature_*.{nc,csv}    subsurface temperature
   %  Globbing is intentionally loose because the release file names carry a
   %  year/version stamp and the user may keep CSV or NetCDF.
   %
   %  Source: SUMup 2024 release, NSIDC G02288.
   %    DOI:        https://doi.org/10.18739/A2M61BR5M
   %    Sub-DOIs:   density       10.18739/A2JH3D23R
   %                accumulation  10.18739/A2DR2P790
   %    NSIDC:      https://nsidc.org/data/g02288
   %
   %  Behaviour
   %    - For each required SUMup variable (density, accumulation,
   %      temperature), glob-match at least one CSV or NetCDF file.
   %    - On success, return the cache directory so the caller can pass it
   %      to icemodel.verification.setup.importSumup.
   %    - On any missing variable, print actionable NASA Earthdata retrieval
   %      instructions (DOI, NSIDC URL, Earthdata login) and either error
   %      (kwargs.strict=true, default) or return the partial cache
   %      directory (kwargs.strict=false).
   %    - Does NOT attempt automatic download. SUMup is access-gated behind
   %      a NASA Earthdata Login (registration), so the retrieval step is
   %      made explicit rather than a silent failure mode in CI.
   %
   %  Role
   %    Validator. The fetch helper guarantees the cache directory exists and
   %    that every required SUMup variable file is present, so downstream
   %    importers / builders can assume the layout is correct without
   %    repeating per-file checks. This mirrors fetchEsmSnowmip.
   %
   %  Name-value
   %    cache_dir : string (default data/verification/sumup)
   %        Local source-cache directory.
   %    variables : string vector (default density/accumulation/temperature)
   %        Required SUMup variable groups to verify.
   %    strict : logical (default true)
   %        Error when any required variable file is missing.
   %    silent : logical (default false)
   %        Suppress the retrieval-instructions printout when files are
   %        missing.
   %
   %  Returns
   %    source_dir : string
   %        Absolute path to the cache directory.
   %
   % See also: icemodel.verification.setup.importSumup,
   %  icemodel.verification.setup.fetchEsmSnowmip

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.variables (1, :) string = ...
         ["density", "accumulation", "temperature"]
      kwargs.strict   (1, 1) logical = true
      kwargs.silent   (1, 1) logical = false
   end

   cache_dir = kwargs.cache_dir;

   % Ensure the cache directory exists so users following the retrieval banner
   % can drop files into a path that is already there.
   icemodel.helpers.ensureDirExists(cache_dir);

   % Per-variable presence check. Each SUMup variable group needs at least one
   % CSV or NetCDF file matching the release naming pattern.
   missing = missingVariableFiles(cache_dir, kwargs.variables);

   ok = isempty(missing);
   if ok
      source_dir = string(cache_dir);
      return
   end

   % Print actionable retrieval instructions when files are missing. The banner
   % format is stable so developers can grep for it and so it cannot be
   % confused with a regular error.
   if ~kwargs.silent
      fprintf('\n');
      fprintf('=== SUMup firn source cache incomplete ===\n');
      fprintf('Cache directory: %s\n', cache_dir);
      fprintf('Missing variable file patterns:\n');
      for j = 1:numel(missing)
         fprintf('  - %s\n', missing(j));
      end
      fprintf('\nRetrieval (NASA Earthdata, access-gated):\n');
      fprintf('  Dataset:   SUMup 2024 release, NSIDC G02288\n');
      fprintf('  Data DOI:  https://doi.org/10.18739/A2M61BR5M\n');
      fprintf('  Sub-DOIs:  density      10.18739/A2JH3D23R\n');
      fprintf('             accumulation 10.18739/A2DR2P790\n');
      fprintf('  NSIDC:     https://nsidc.org/data/g02288\n');
      fprintf('  Login:     requires a (free) NASA Earthdata Login.\n');
      fprintf('             Register at https://urs.earthdata.nasa.gov/\n');
      fprintf('  Manual workflow: log in, download the density / accumulation\n');
      fprintf('  / temperature CSV or NetCDF files, and place them into the\n');
      fprintf('  cache directory above.\n');
      fprintf('\nAfter retrieval, re-run:\n');
      fprintf('  icemodel.verification.setup.fetchSumup()\n');
      fprintf('  icemodel.verification.setup.importSumup(source_dir, overwrite=true)\n');
      fprintf('\n');
   end

   if kwargs.strict
      error('icemodel:verification:fetchSumup:missingSources', ...
         ['SUMup firn source cache incomplete in %s. Missing: %s. ' ...
         'See retrieval instructions above.'], ...
         cache_dir, strjoin(missing, ', '));
   end

   source_dir = string(cache_dir);
end

%% Local helpers
function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical SUMup source-cache directory.
   %
   % icemodel.getpath('data') returns the canonical top-level data root
   % (<repo>/data/), which the verification source-cache layout extends under
   % data/verification/<dataset_family>/ (family-flat taxonomy; the same
   % layout used by ESM-SnowMIP and the other families).
   pathname = string(fullfile(icemodel.getpath('data'), ...
      'verification', 'sumup'));
end

function missing = missingVariableFiles(cache_dir, variables)
   %MISSINGVARIABLEFILES Collect the SUMup variable groups with no cache file.
   %
   % Each required variable contributes one expected glob pattern. A variable
   % is satisfied by any matching CSV or NetCDF file. The return array is
   % preallocated to the worst case (one per variable) then trimmed.

   variables = reshape(variables, 1, []);
   n_vars = numel(variables);
   missing = strings(n_vars, 1);
   n_missing = 0;

   for k = 1:n_vars
      v = variables(k);
      nc = dir(fullfile(cache_dir, sprintf('*%s*.nc', v)));
      csv = dir(fullfile(cache_dir, sprintf('*%s*.csv', v)));
      if isempty(nc) && isempty(csv)
         n_missing = n_missing + 1;
         missing(n_missing) = sprintf('*%s*.{nc,csv}', v);
      end
   end

   missing = missing(1:n_missing);
end
