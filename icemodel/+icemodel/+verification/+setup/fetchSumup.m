function [source_dir, status] = fetchSumup(kwargs)
   %FETCHSUMUP Locate or verify the SUMup firn source files.
   %
   %  source_dir = icemodel.verification.setup.fetchSumup()
   %  source_dir = icemodel.verification.setup.fetchSumup( ...
   %     cache_dir="/some/other/path", strict=true)
   %
   %  Resolves the local source-cache directory holding the SUMup firn
   %  observation files (density, SMB, subsurface temperature) used by
   %  icemodel.verification.setup.importSumup to stage the per-point firn
   %  verification artifacts. By default the cache lives at
   %  data/verification/sumup/ (committed in this repo) and is populated by the
   %  user / developer following the retrieval instructions printed when
   %  files are missing.
   %
   %  SUMup is the firn observation source for density, SMB, and subsurface
   %  temperature. FirnCover compaction strain is not part of this cache.
   %
   %  Expected files (SUMup 2025 release; NetCDF, per ice sheet):
   %    SUMup_2025_density_greenland.nc       firn/snow density profiles
   %    SUMup_2025_SMB_greenland.nc           surface mass balance
   %    SUMup_2025_temperature_greenland.nc   subsurface temperature
   %  Antarctic counterparts (*_antarctica.nc) ship in the same release but are
   %  out of scope here. Globbing remains loose (per variable, *.nc/*.csv) so a
   %  future release stamp or a CSV export still satisfies the cache check.
   %
   %  Source: SUMup 2025 release, Arctic Data Center / NSIDC G02288.
   %    DOI:        https://doi.org/10.18739/A2M61BR5M
   %    NSIDC:      https://nsidc.org/data/g02288
   %    Contact:    Baptiste Vandecrux (bav@geus.dk)
   %
   %  Behaviour
   %    - For each required SUMup variable (density, SMB,
   %      temperature), glob-match at least one CSV or NetCDF file.
   %    - On success, return the cache directory so the caller can pass it
   %      to icemodel.verification.setup.importSumup.
   %    - On any missing variable, print actionable NASA Earthdata retrieval
   %      instructions (DOI, NSIDC URL, Earthdata login) and either error
   %      (kwargs.strict=true, default) or return the partial cache
   %      directory (kwargs.strict=false).
   %    - Does NOT attempt automatic download. SUMup is access-gated behind
   %      a NASA Earthdata Login (registration), so the retrieval step is
   %      made explicit rather than hidden behind an automatic download.
   %
   %  Role
   %    Validator. With create_cache_dir=true, the fetch helper creates the
   %    cache directory before checking every required SUMup variable. With
   %    false, it reports the same status without mutation. A successful strict
   %    result lets downstream importers/builders assume the layout is complete.
   %    An empty variables selection returns an empty status without creating or
   %    scanning cache_dir.
   %
   %  Name-value
   %    cache_dir : string (default data/verification/sumup)
   %        Local source-cache directory.
   %    variables : string vector (default density/SMB/temperature)
   %        Required SUMup variable groups to verify.
   %    region : string (default "greenland")
   %        Ice-sheet region tag used to scope the per-variable glob to the
   %        Greenland files (the SUMup release ships *_greenland.nc and
   %        *_antarctica.nc per variable).
   %    strict : logical (default true)
   %        Error when any required variable file is missing.
   %    silent : logical (default false)
   %        Suppress the retrieval-instructions printout when files are
   %        missing.
   %    create_cache_dir : logical (default true)
   %        Create the resolved cache directory before validation.
   %
   %  Returns
   %    source_dir : string
   %        Absolute path to the cache directory.
   %    status : struct array
   %        Shared fetch-status rows describing missing cache items.
   %
   % See also: icemodel.verification.setup.importSumup,
   %  icemodel.verification.setup.fetchEsmSnowmip

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.variables (1, :) string ...
         {mustBeKnownVariables(kwargs.variables)} = variableRegistry()
      kwargs.region   (1, 1) string = "greenland"
      kwargs.strict   (1, 1) logical = true
      kwargs.silent   (1, 1) logical = false
      kwargs.create_cache_dir (1, 1) logical = true
   end

   cache_dir = icemodel.verification.setup.resolveFetchCacheDir( ...
      kwargs.cache_dir, defaultCacheDir());

   % Cache creation is explicit so dry-run callers can remain non-mutating.
   if kwargs.create_cache_dir && ~isempty(kwargs.variables)
      icemodel.helpers.ensureDirExists(cache_dir);
   end

   % Build cache status before strict handling.
   missing = missingVariableFiles(cache_dir, kwargs.variables, kwargs.region);
   status = icemodel.verification.setup.fetchMissingStatus(missing);
   source_dir = icemodel.verification.setup.finishFetchStatus( ...
      cache_dir, status, strict=kwargs.strict, silent=kwargs.silent, ...
      error_id="icemodel:verification:fetchSumup:missingSources", ...
      error_label="SUMup", ...
      banner_callback=@(~, ~, items) printRetrievalBanner(cache_dir, items), ...
      error_callback=@(~, ~, items) throwMissingSources(cache_dir, items));
end

%% Local helpers
function printRetrievalBanner(cache_dir, missing)
   %PRINTRETRIEVALBANNER Print SUMup retrieval instructions.
   fprintf('\n');
   fprintf('=== SUMup firn source cache incomplete ===\n');
   fprintf('Cache directory: %s\n', cache_dir);
   fprintf('Missing variable file patterns:\n');
   for j = 1:numel(missing)
      fprintf('  - %s\n', missing(j));
   end
   fprintf('\nRetrieval (Arctic Data Center / NSIDC G02288):\n');
   fprintf('  Dataset:   SUMup 2025 release\n');
   fprintf('  Data DOI:  https://doi.org/10.18739/A2M61BR5M\n');
   fprintf('  NSIDC:     https://nsidc.org/data/g02288\n');
   fprintf('  Files:     SUMup_2025_density_greenland.nc\n');
   fprintf('             SUMup_2025_SMB_greenland.nc\n');
   fprintf('             SUMup_2025_temperature_greenland.nc\n');
   fprintf('  Manual workflow: download the density / SMB / temperature\n');
   fprintf('  NetCDF files for the region and place them into the cache\n');
   fprintf('  directory above.\n');
   fprintf('\nAfter retrieval, re-run:\n');
   fprintf('  icemodel.verification.setup.fetchSumup()\n');
   fprintf('  icemodel.verification.setup.importSumup(source_dir, overwrite=true)\n');
   fprintf('\n');
end

function throwMissingSources(cache_dir, missing)
   %THROWMISSINGSOURCES Raise the SUMup missing-cache error.
   error('icemodel:verification:fetchSumup:missingSources', ...
      ['SUMup firn source cache incomplete in %s. Missing: %s. ' ...
      'See retrieval instructions above.'], cache_dir, strjoin(missing, ', '));
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical SUMup source-cache directory.
   pathname = icemodel.forcing.helpers.verificationSourceDir("", "sumup");
end

function missing = missingVariableFiles(cache_dir, variables, region)
   %MISSINGVARIABLEFILES Collect the SUMup variable groups with no cache file.
   %
   % Each required variable contributes one expected glob pattern keyed by the
   % variable name and the region tag (e.g. *density*greenland*.nc). A variable
   % is satisfied by any matching CSV or NetCDF file. The return array is
   % preallocated to the worst case (one per variable) then trimmed.

   variables = reshape(variables, 1, []);
   region = char(region);
   n_vars = numel(variables);
   missing = strings(n_vars, 1);
   n_missing = 0;

   for k = 1:n_vars
      v = char(variables(k));
      nc = dir(fullfile(cache_dir, sprintf('*%s*%s*.nc', v, region)));
      csv = dir(fullfile(cache_dir, sprintf('*%s*%s*.csv', v, region)));
      if isempty(nc) && isempty(csv)
         n_missing = n_missing + 1;
         missing(n_missing) = sprintf('*%s*%s*.{nc,csv}', v, region);
      end
   end

   missing = missing(1:n_missing);
end

function variables = variableRegistry()
   %VARIABLEREGISTRY Return supported SUMup variable groups in fetch order.
   variables = ["density", "SMB", "temperature"];
end

function mustBeKnownVariables(variables)
   %MUSTBEKNOWNVARIABLES Reject typos before cache creation or discovery.
   bad = setdiff(reshape(string(variables), 1, []), ...
      variableRegistry(), 'stable');
   if ~isempty(bad)
      error('icemodel:verification:fetchSumup:unknownVariable', ...
         'unknown SUMup variable group(s): %s', strjoin(bad, ', '));
   end
end
