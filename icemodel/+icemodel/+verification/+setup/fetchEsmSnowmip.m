function [source_dir, status] = fetchEsmSnowmip(kwargs)
   %FETCHESMSNOWMIP Locate or verify the ESM-SnowMIP source NetCDF files.
   %
   %  source_dir = icemodel.verification.setup.fetchEsmSnowmip()
   %  source_dir = icemodel.verification.setup.fetchEsmSnowmip( ...
   %     cache_dir="/some/other/path", strict=true)
   %
   %  Resolves the local source-cache directory holding the ESM-SnowMIP
   %  meteorological / observation NetCDF files used by
   %  icemodel.verification.setup.importEsmSnowmip to stage the per-site
   %  verification artifacts. By default the cache lives at
   %  data/verification/esm_snowmip/ (gitignored) and is populated
   %  by the user / developer following the retrieval instructions
   %  printed when files are missing.
   %
   %  Expected files (per site, all 10 ESM-SnowMIP reference sites):
   %    met_insitu_<sitename>_<years>.nc    (forcing)
   %    obs_insitu_<sitename>_<years>.nc    (observations)
   %
   %  Source: Menard et al. 2019, "Meteorological and evaluation
   %  datasets for snow modelling at 10 reference sites: description
   %  of in situ and bias-corrected reanalysis data", Earth Syst.
   %  Sci. Data, https://doi.org/10.5194/essd-11-865-2019.
   %
   %  Data DOI: https://doi.org/10.1594/PANGAEA.897575
   %
   %  Behaviour
   %    - For each site in the canonical namelist, glob-match one met
   %      and one obs NetCDF and confirm both are well-formed.
   %    - On success, return the cache directory so the caller can
   %      pass it to icemodel.verification.setup.importEsmSnowmip.
   %    - On any missing or unreadable file, print actionable
   %      retrieval instructions (DOI, URL, expected pattern) and
   %      either error (kwargs.strict=true, default) or return the
   %      partial cache directory (kwargs.strict=false).
   %    - Does NOT attempt automatic download. The PANGAEA dataset
   %      access surface is not stable for unattended fetch and may
   %      require user registration / acceptance of terms; making the
   %      retrieval step explicit is preferable to a silent failure
   %      mode in CI.
   %
   %  Role
   %    Validator. With create_cache_dir=true, the fetch helper creates the
   %    cache directory before checking every site's met/obs pair. With false,
   %    it reports the same status without mutation. A successful strict result
   %    lets downstream importers/builders assume the layout is complete.
   %    An empty stations selection returns an empty status without creating or
   %    scanning cache_dir.
   %
   %  Name-value
   %    cache_dir : string (default data/verification/esm_snowmip)
   %        Local source-cache directory.
   %    stations : string vector (default all 10 ESM-SnowMIP sites)
   %        Restrict the validation to a subset (used by builders that
   %        only need one site).
   %    strict : logical (default true)
   %        Error when any expected file is missing or unreadable.
   %    silent : logical (default false)
   %        Suppress the retrieval-instructions printout when files
   %        are missing.
   %    create_cache_dir : logical (default true)
   %        Create the resolved cache directory before validation.
   %
   %  Returns
   %    source_dir : string
   %        Absolute path to the cache directory.
   %    status : struct array
   %        Shared fetch-status rows describing missing cache items.
   %
   % See also: icemodel.verification.setup.importEsmSnowmip,
   %  icemodel.verification.setup.fetchLaughTests,
   %  icemodel.verification.namelists.snowmipsite

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.stations (1, :) string = ...
         icemodel.verification.namelists.snowmipsite()
      kwargs.strict   (1, 1) logical = true
      kwargs.silent   (1, 1) logical = false
      kwargs.create_cache_dir (1, 1) logical = true
   end

   cache_dir = icemodel.verification.setup.resolveFetchCacheDir( ...
      kwargs.cache_dir, defaultCacheDir());

   % Cache creation is explicit so dry-run callers can remain non-mutating.
   if kwargs.create_cache_dir && ~isempty(kwargs.stations)
      icemodel.helpers.ensureDirExists(cache_dir);
   end

   % Build cache status before strict handling.
   [missing, broken] = missingOrBrokenFiles(cache_dir, kwargs.stations);
   failures = [missing(:); broken(:)];
   status = icemodel.verification.setup.fetchMissingStatus(failures);
   source_dir = icemodel.verification.setup.finishFetchStatus( ...
      cache_dir, status, strict=kwargs.strict, silent=kwargs.silent, ...
      error_id="icemodel:verification:fetchEsmSnowmip:missingSources", ...
      error_label="ESM-SnowMIP", ...
      banner_callback=@(root, ~, ~) printRetrievalBanner(root, missing, broken), ...
      error_callback=@(root, ~, ~) throwMissingSources(root, missing, broken));
end

%% Local helpers
function printRetrievalBanner(cache_dir, missing, broken)
   %PRINTRETRIEVALBANNER Print ESM-SnowMIP retrieval instructions.
   fprintf('\n');
   fprintf('=== ESM-SnowMIP source cache incomplete ===\n');
   fprintf('Cache directory: %s\n', cache_dir);
   if ~isempty(missing)
      fprintf('Missing file patterns:\n');
      for j = 1:numel(missing)
         fprintf('  - %s\n', missing(j));
      end
   end
   if ~isempty(broken)
      fprintf('Unreadable NetCDF files (partial download):\n');
      for j = 1:numel(broken)
         fprintf('  - %s\n', broken(j));
      end
   end
   fprintf('\nRetrieval:\n');
   fprintf('  Reference: Menard et al. 2019, ESSD\n');
   fprintf('             https://doi.org/10.5194/essd-11-865-2019\n');
   fprintf('  Data DOI:  https://doi.org/10.1594/PANGAEA.897575\n');
   fprintf('  Bundle:    "All ESM-SnowMIP netCDF files in one zip" on PANGAEA.\n');
   fprintf('  Manual workflow: download the zip and extract met_insitu_*.nc\n');
   fprintf('  and obs_insitu_*.nc into the cache directory.\n');
   fprintf('\nAfter retrieval, re-run:\n');
   fprintf('  icemodel.verification.setup.fetchEsmSnowmip()\n');
   fprintf('  icemodel.verification.setup.importEsmSnowmip(source_dir, overwrite=true)\n');
   fprintf('\n');
end

function throwMissingSources(cache_dir, missing, broken)
   %THROWMISSINGSOURCES Raise the ESM-SnowMIP missing-cache error.
   error('icemodel:verification:fetchEsmSnowmip:missingSources', ...
      ['ESM-SnowMIP source cache incomplete in %s. Missing: %s. ' ...
      'Unreadable: %s. See retrieval instructions above.'], ...
      cache_dir, strjoin(missing, ', '), strjoin(broken, ', '));
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical ESM-SnowMIP source-cache directory.
   %
   % The gitignored raw source cache is a repo-root developer resource at
   % <repo>/data/verification/<dataset_family>/. Resolve it with
   % icemodel.internal.fullpath('data') (always the repo-local data root), NOT
   % icemodel.getpath('data'): getpath returns the ICEMODEL_DATA_PATH env var,
   % which a scoped demo or test case sets to its owned data tree rather than
   % the repo-root verification workspace. (When ICEMODEL_DATA_PATH is unset
   % the two coincide via internal.fullpath.)
   pathname = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'esm_snowmip'));
end

function [missing, broken] = missingOrBrokenFiles(cache_dir, stations)
   %MISSINGORBROKENFILES Collect missing-pattern / unreadable-file lists.
   %
   % Each site contributes one met and one obs expected pattern. The
   % return arrays are preallocated to the worst case (2 * n_sites) so
   % no growth in the loop, then trimmed.

   stations = reshape(stations, 1, []);
   n_sites = numel(stations);
   n_max = 2 * n_sites;
   missing = strings(n_max, 1);
   broken  = strings(n_max, 1);
   n_missing = 0;
   n_broken  = 0;

   patterns = strings(n_max, 1);
   for k = 1:n_sites
      patterns(2 * k - 1) = sprintf('met_insitu_%s_*.nc', stations(k));
      patterns(2 * k)     = sprintf('obs_insitu_%s_*.nc', stations(k));
   end

   for k = 1:n_max
      pattern = patterns(k);
      matches = dir(fullfile(cache_dir, pattern));
      if isempty(matches)
         n_missing = n_missing + 1;
         missing(n_missing) = pattern;
         continue
      end
      % Use the first match (the upstream bundle has a single file per
      % site; ambiguous matches are surfaced by the importer).
      try
         ncinfo(fullfile(matches(1).folder, matches(1).name));
      catch
         n_broken = n_broken + 1;
         broken(n_broken) = string(matches(1).name);
      end
   end

   missing = missing(1:n_missing);
   broken  = broken(1:n_broken);
end
