function source_dir = fetchEsmSnowmip(kwargs)
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
   %  data/verification/snow/esm_snowmip/ (gitignored) and is populated
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
   %    Validator. The fetch helper guarantees the cache directory
   %    exists and that every site's met/obs file pair is present and
   %    readable, so downstream importers / builders can assume the
   %    layout is correct without repeating per-file checks.
   %
   %  Name-value
   %    cache_dir : string (default data/verification/snow/esm_snowmip)
   %        Local source-cache directory.
   %    sitenames : string vector (default all 10 ESM-SnowMIP sites)
   %        Restrict the validation to a subset (used by builders that
   %        only need one site).
   %    strict : logical (default true)
   %        Error when any expected file is missing or unreadable.
   %    silent : logical (default false)
   %        Suppress the retrieval-instructions printout when files
   %        are missing.
   %
   %  Returns
   %    source_dir : string
   %        Absolute path to the cache directory.
   %
   % See also: icemodel.verification.setup.importEsmSnowmip,
   %  icemodel.verification.setup.fetchLaughTests,
   %  icemodel.verification.namelists.snowmipsite

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.sitenames (1, :) string = ...
         icemodel.verification.namelists.snowmipsite()
      kwargs.strict   (1, 1) logical = true
      kwargs.silent   (1, 1) logical = false
   end

   cache_dir = kwargs.cache_dir;

   % Ensure the cache directory exists so users following the retrieval
   % banner can drop files into a path that is already there.
   icemodel.helpers.ensureDirExists(cache_dir);

   % Per-site presence + readability check. Each ESM-SnowMIP site needs one met
   % (forcing) and one obs (evaluation) NetCDF; both are located by glob match
   % against the canonical PANGAEA naming pattern so year-stamped filenames do
   % not need to be hard-coded.
   [missing, broken] = missingOrBrokenFiles(cache_dir, kwargs.sitenames);

   ok = isempty(missing) && isempty(broken);
   if ok
      source_dir = string(cache_dir);
      return
   end

   % Print actionable retrieval instructions when files are missing
   % or broken. The banner format is stable so developers can grep
   % for it and so it cannot be confused with a regular error.
   if ~kwargs.silent
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

   if kwargs.strict
      error('icemodel:verification:fetchEsmSnowmip:missingSources', ...
         ['ESM-SnowMIP source cache incomplete in %s. Missing: %s. ' ...
         'Unreadable: %s. See retrieval instructions above.'], ...
         cache_dir, ...
         strjoin(missing, ', '), ...
         strjoin(broken, ', '));
   end

   source_dir = string(cache_dir);
end

%% Local helpers
function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical ESM-SnowMIP source-cache directory.
   %
   % icemodel.getpath('data') returns the canonical top-level data root
   % (<repo>/data/), which the verification source-cache layout extends
   % under data/verification/snow/<dataset_family>/.
   pathname = string(fullfile(icemodel.getpath('data'), ...
      'verification', 'snow', 'esm_snowmip'));
end

function [missing, broken] = missingOrBrokenFiles(cache_dir, sitenames)
   %MISSINGORBROKENFILES Collect missing-pattern / unreadable-file lists.
   %
   % Each site contributes one met and one obs expected pattern. The
   % return arrays are preallocated to the worst case (2 * n_sites) so
   % no growth in the loop, then trimmed.

   sitenames = reshape(sitenames, 1, []);
   n_sites = numel(sitenames);
   n_max = 2 * n_sites;
   missing = strings(n_max, 1);
   broken  = strings(n_max, 1);
   n_missing = 0;
   n_broken  = 0;

   patterns = strings(n_max, 1);
   for k = 1:n_sites
      patterns(2 * k - 1) = sprintf('met_insitu_%s_*.nc', sitenames(k));
      patterns(2 * k)     = sprintf('obs_insitu_%s_*.nc', sitenames(k));
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
