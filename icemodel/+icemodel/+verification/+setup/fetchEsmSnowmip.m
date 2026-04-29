function source_dir = fetchEsmSnowmip(kwargs)
   %FETCHESMSNOWMIP Locate or verify the ESM-SnowMIP source NetCDF files.
   %
   %  source_dir = icemodel.verification.setup.fetchEsmSnowmip()
   %  source_dir = icemodel.verification.setup.fetchEsmSnowmip( ...
   %     cache_dir="/some/other/path", strict=true)
   %
   %  Resolves the local source-cache directory holding the ESM-SnowMIP
   %  meteorological / observation NetCDF files used by
   %  icemodel.verification.setup.importEsmSnowmip to stage the cdp /
   %  wfj smoke verification artifacts. By default the cache lives at
   %  data/verification/snow/esm_snowmip/ (gitignored) and is populated
   %  by the user / developer following the retrieval instructions
   %  printed when files are missing.
   %
   %  Files expected:
   %    met_insitu_cdp_1994_2014.nc    (Col de Porte forcing)
   %    obs_insitu_cdp_1994_2014.nc    (Col de Porte observations)
   %    met_insitu_wfj_1996_2016.nc    (Weissfluhjoch forcing)
   %    obs_insitu_wfj_1996_2016.nc    (Weissfluhjoch observations)
   %
   %  Source: Menard et al. 2019, "Meteorological and evaluation
   %  datasets for snow modelling at 10 reference sites: description
   %  of in situ and bias-corrected reanalysis data", Earth Syst.
   %  Sci. Data, https://doi.org/10.5194/essd-11-865-2019.
   %
   %  Data DOI: https://doi.org/10.1594/PANGAEA.897575
   %
   %  Behaviour
   %    - If all four files exist in cache_dir and are valid NetCDFs,
   %      returns the cache directory path so the caller can pass it
   %      to icemodel.verification.setup.importEsmSnowmip(source_dir).
   %    - If any file is missing or unreadable, prints actionable
   %      retrieval instructions (DOI, URL, expected filename, target
   %      path, optional manual-download workflow) and either errors
   %      (kwargs.strict=true, default) or returns the partial cache
   %      directory (kwargs.strict=false) so callers can decide how
   %      to handle the missing-data state.
   %    - Does NOT attempt automatic download. The PANGAEA dataset
   %      access surface is not stable for unattended fetch and may
   %      require user registration / acceptance of terms; making the
   %      retrieval step explicit is preferable to a silent failure
   %      mode in CI.
   %
   %  Name-value
   %    cache_dir : string (default data/verification/snow/esm_snowmip)
   %        Local source-cache directory.
   %    strict : logical (default true)
   %        Error when any expected file is missing or unreadable.
   %        Set false to obtain the path with a warning instead.
   %    silent : logical (default false)
   %        Suppress the retrieval-instructions printout when files
   %        are missing. Use only in tests; user-facing calls should
   %        leave the actionable message visible.
   %
   %  Returns
   %    source_dir : string
   %        Absolute path to the cache directory.
   %
   % See also: icemodel.verification.setup.importEsmSnowmip,
   %  icemodel.verification.setup.importLaughTests

   arguments
      kwargs.cache_dir (1, 1) string = ""
      kwargs.strict   (1, 1) logical = true
      kwargs.silent   (1, 1) logical = false
   end

   % Resolve the default cache directory. Project root is found by
   % walking up from this file's location (icemodel/+icemodel/...
   % +setup/fetchEsmSnowmip.m), which is robust against the caller's
   % current working directory.
   if kwargs.cache_dir == ""
      cache_dir = string(fullfile(repoRoot(), 'data', 'verification', ...
         'snow', 'esm_snowmip'));
   else
      cache_dir = kwargs.cache_dir;
   end

   % Expected files. Sourced from defaultSpecs() in importEsmSnowmip:
   % each ESM-SnowMIP smoke-case site needs one met (forcing) NetCDF
   % and one obs (evaluation) NetCDF.
   expected = struct( ...
      'site',     {"cdp",                          "cdp",                          "wfj",                          "wfj"}, ...
      'role',     {"met",                          "obs",                          "met",                          "obs"}, ...
      'filename', {"met_insitu_cdp_1994_2014.nc",  "obs_insitu_cdp_1994_2014.nc",  "met_insitu_wfj_1996_2016.nc",  "obs_insitu_wfj_1996_2016.nc"});

   % Ensure the cache directory exists before checking files. The
   % check is non-destructive: if the directory does not exist yet,
   % it is created, but no files are written.
   if exist(cache_dir, 'dir') ~= 7
      icemodel.helpers.ensureDirExists(cache_dir);
   end

   % Per-file presence + readability check. ncinfo is the cheapest
   % way to confirm the NetCDF file is well-formed without loading
   % the bulk of its data.
   missing  = strings(0, 1);
   broken   = strings(0, 1);
   for i = 1:numel(expected)
      pathname = fullfile(cache_dir, expected(i).filename);
      if exist(pathname, 'file') ~= 2
         missing(end + 1, 1) = expected(i).filename; %#ok<AGROW>
         continue
      end
      try
         ncinfo(pathname);
      catch
         broken(end + 1, 1) = expected(i).filename; %#ok<AGROW>
      end
   end

   ok = isempty(missing) && isempty(broken);
   if ok
      source_dir = string(cache_dir);
      return
   end

   % Print actionable retrieval instructions when files are missing
   % or broken. Keep the message stable so future agents can grep
   % for it and so it cannot be confused with a regular error.
   if ~kwargs.silent
      fprintf('\n');
      fprintf('=== ESM-SnowMIP source cache incomplete ===\n');
      fprintf('Cache directory: %s\n', cache_dir);
      if ~isempty(missing)
         fprintf('Missing files:\n');
         for j = 1:numel(missing)
            fprintf('  - %s\n', missing(j));
         end
      end
      if ~isempty(broken)
         fprintf('Unreadable NetCDF files (likely partial download):\n');
         for j = 1:numel(broken)
            fprintf('  - %s\n', broken(j));
         end
      end
      fprintf('\nRetrieval:\n');
      fprintf('  Reference: Menard et al. 2019, ESSD\n');
      fprintf('             https://doi.org/10.5194/essd-11-865-2019\n');
      fprintf('  Data DOI:  https://doi.org/10.1594/PANGAEA.897575\n');
      fprintf('  Manual workflow: download the four files above from the\n');
      fprintf('  PANGAEA dataset and place them in the cache directory.\n');
      fprintf('  Filenames must match exactly.\n');
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

function root = repoRoot()
   %REPOROOT Resolve the icemodel project root.
   %
   %  Walks up from this file's path until a directory containing
   %  the canonical project marker (the inner 'icemodel' source
   %  directory) is found. This is robust against MATLAB's package-
   %  path resolution and against the caller's current working
   %  directory.
   here = fileparts(mfilename('fullpath'));
   root = here;
   while ~isempty(root)
      if exist(fullfile(root, 'icemodel', '+icemodel'), 'dir') == 7
         return
      end
      parent = fileparts(root);
      if strcmp(parent, root)
         break
      end
      root = parent;
   end
   error('icemodel:verification:fetchEsmSnowmip:noRepoRoot', ...
      'unable to resolve icemodel project root from %s', ...
      fileparts(mfilename('fullpath')));
end
