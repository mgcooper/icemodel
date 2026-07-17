function [observations, metadata] = buildEsmSnowmipObservations(sitename, kwargs)
   %BUILDESMSNOWMIPOBSERVATIONS Convert ESM-SnowMIP obs NetCDF to verification targets.
   %
   %  [observations, metadata] = buildEsmSnowmipObservations(sitename)
   %  [observations, metadata] = buildEsmSnowmipObservations(sitename, ...
   %     source_dir=..., startdate=..., enddate=...)
   %
   %  Reads one ESM-SnowMIP obs_insitu_<sitename>_*.nc file and returns
   %  the observed snow / surface variables as a verification-target
   %  timetable. Variable mapping uses the verification-target schema consumed
   %  by comparecase, plotcase, and the runIcemodelSnowCandidate adapter without
   %  a per-site branch.
   %
   %  Observation contract: this builder stages every snow / surface
   %  channel that the upstream obs file makes available, so downstream
   %  model-vs-observation comparison can use any of them. Canonical
   %  comparison columns (snow_depth_m, swe_kg_m2, surface_temp_C,
   %  soil_temp_<k>_C) are added when the corresponding upstream
   %  channel is present in the NetCDF, regardless of whether the
   %  smoke window happens to contain finite values for that channel
   %  (so plotcase still draws a sparse-marker axis when observations
   %  are gappy in the staged window).
   %
   %  Variable mapping (ESM-SnowMIP -> verification target):
   %     snd_auto / snd_man / snd_gap_auto / snd_gap1_auto
   %        -> canonical snow_depth_m  [m]   (best-channel selector)
   %        -> snd_*_m                 [m]   (all available variants)
   %     snw_auto / snw_man -> swe_kg_m2     [kg m-2]
   %     ts                 -> surface_temp_C [C]
   %     tsl(:, k)          -> soil_temp_<k>_C [C] for each soil depth
   %     albs               -> albedo        [1]   (when present)
   %
   %  Inputs
   %    sitename : string  ESM-SnowMIP site code (e.g. "cdp"). Validated
   %                       against icemodel.verification.namelists.snowmipsite.
   %
   %  Name-value
   %    source_dir : string (default data/verification/esm_snowmip)
   %    startdate  : datetime or "" (default "")
   %    enddate    : datetime or "" (default "")
   %
   %  Returns
   %    observations : timetable (verification target schema)
   %    metadata     : struct with provenance + soil depth metadata
   %
   %  Role
   %    Reusable per-site observation builder. Symmetric with
   %    buildEsmSnowmipForcing. Used by importEsmSnowmip during staging.
   %
   % See also: icemodel.verification.setup.buildEsmSnowmipForcing,
   %  icemodel.verification.setup.importEsmSnowmip,
   %  icemodel.verification.setup.readNetcdfVariable,
   %  icemodel.verification.setup.readBestSnowDepth

   arguments
      sitename (1, 1) string ...
         {icemodel.verification.validators.mustBeSnowmipSite}
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate  = ""
      kwargs.enddate    = ""
   end
   % Validate and normalize the optional interval before source discovery so
   % malformed public input cannot be obscured by an unrelated file error.
   [window_start, window_end, has_window] = ...
      icemodel.internal.pairedWindow(kwargs.startdate, kwargs.enddate);

   % Treat an explicit blank like omission; never reinterpret it as the CWD.
   source_dir = icemodel.forcing.helpers.verificationSourceDir( ...
      kwargs.source_dir, "esm_snowmip");

   % Locate the obs file (one per site).
   pattern = sprintf('obs_insitu_%s_*.nc', sitename);
   matches = dir(fullfile(source_dir, pattern));
   if isempty(matches)
      error('icemodel:verification:buildEsmSnowmipObservations:fileNotFound', ...
         'no file matching %s in %s', pattern, source_dir);
   end
   if numel(matches) > 1
      names = strjoin({matches.name}, ', ');
      error('icemodel:verification:buildEsmSnowmipObservations:ambiguousFile', ...
         'multiple files matching %s in %s: %s', pattern, source_dir, names);
   end
   obsfile = fullfile(matches.folder, matches.name);

   % --- Read the obs time axis ----------------------------------------
   obs_time = icemodel.verification.setup.readNetcdfTime(obsfile, 'time');
   ntime = numel(obs_time);

   % Discover which channels the upstream file provides. Canonical
   % comparison columns are added only when the matching upstream
   % channel(s) exist. ESM-SnowMIP sites are heterogeneous:
   %   - boreal forest sites (oas, obs, ojp) report snd_gap_auto in the
   %     canopy gap rather than snd_auto;
   %   - sap lacks tsl, so no soil_temp_<k>_C columns are added;
   %   - rme and sod lack albs (no observed albedo channel);
   %   - sap lacks the sdepth dimension entirely;
   %   - sod lacks ts, so surface_temp_C is omitted.
   info = ncinfo(obsfile);
   channels = string({info.Variables.Name});
   has_snow_depth = any(ismember(channels, ...
      ["snd_auto", "snd_gap_auto", "snd_gap1_auto", "snd_man"]));
   has_swe = any(ismember(channels, ["snw_auto", "snw_man"]));
   has_ts = any(channels == "ts");
   has_tsl = any(channels == "tsl");
   has_albs = any(channels == "albs");

   % --- Read canonical channels --------------------------------------

   % snow depth
   if has_snow_depth
      [snd_best, snow_depth_source] = ...
         icemodel.verification.setup.readBestSnowDepth(obsfile);
   else
      snd_best = nan(ntime, 1);
      snow_depth_source = "none";
   end

   % SWE
   if has_swe
      swe_auto = icemodel.verification.setup.readNetcdfVariable( ...
         obsfile, "snw_auto", optional=true, ntime=ntime);
      swe_man  = icemodel.verification.setup.readNetcdfVariable( ...
         obsfile, "snw_man", optional=true, ntime=ntime);
      swe_kg_m2 = icemodel.verification.setup.preferPrimary(swe_auto, swe_man);
   else
      swe_kg_m2 = nan(ntime, 1);
   end

   % surface temperature
   if has_ts
      ts = icemodel.verification.setup.readNetcdfVariable(obsfile, "ts");
   else
      ts = nan(ntime, 1);
   end

   % albedo
   if has_albs
      albs = icemodel.verification.setup.readNetcdfVariable(obsfile, "albs");
   else
      albs = nan(ntime, 1);
   end

   % MATLAB reads the native tsl dimensions as sdepth x time. Preserve that
   % source orientation when normalizing singleton dimensions; each row must
   % remain one physical soil depth, so this operation intentionally does not
   % transpose the source data. Sites without tsl produce zero soil-temperature
   % columns rather than a placeholder NaN soil layer.
   if has_tsl
      tsl = icemodel.verification.setup.readNetcdfVariable(obsfile, "tsl");
      tsl = reshape(tsl, [], ntime);
      % Source depths are stored as single precision in NetCDF. Round after
      % conversion so JSON/MAT metadata records useful physical depths rather
      % than binary floating-point noise such as 0.10000000149.
      sdepth = round(double(ncread(obsfile, 'sdepth')), 4);
   else
      tsl = zeros(0, ntime);
      sdepth = zeros(0, 1);
   end

   % Read every available snow-depth variant separately so model-obs
   % comparison can target any specific channel. Variants absent from
   % the file are skipped (we do not synthesise NaN columns for those).
   snd_variants = struct();
   variant_names = ["snd_auto", "snd_man", "snd_gap_auto", "snd_gap1_auto"];
   for c = variant_names
      if any(channels == c)
         snd_variants.(c) = icemodel.verification.setup.readNetcdfVariable( ...
            obsfile, c);
      end
   end

   % --- Optional time-window subset -----------------------------------
   idx = true(ntime, 1);
   if has_window
      idx = obs_time >= window_start & obs_time <= window_end;
   end
   if ~any(idx)
      error('icemodel:verification:buildEsmSnowmipObservations:emptyWindow', ...
         'no observation samples in the requested window for site %s', sitename);
   end

   Time = obs_time(idx);

   % --- Assemble the staged timetable --------------------------------
   % Build canonical columns first so downstream code sees the expected
   % ordering, then append site-specific channels. Count the present channels
   % first and preallocate the variable_names / values vectors.
   tsl_window = tsl(:, idx);
   n_soil = size(tsl_window, 1);
   n_variants = 0;
   for c = variant_names
      if isfield(snd_variants, c)
         n_variants = n_variants + 1;
      end
   end
   n_cols = double(has_snow_depth) + double(has_swe) + double(has_ts) ...
      + n_soil + n_variants + double(has_albs);

   variable_names  = strings(n_cols, 1);
   variable_values = cell(n_cols, 1);
   col = 0;

   if has_snow_depth
      col = col + 1;
      variable_names(col) = "snow_depth_m";
      variable_values{col} = snd_best(idx);
   end
   if has_swe
      col = col + 1;
      variable_names(col) = "swe_kg_m2";
      variable_values{col} = swe_kg_m2(idx);
   end
   if has_ts
      col = col + 1;
      variable_names(col) = "surface_temp_C";
      variable_values{col} = ts(idx);
   end

   for n = 1:n_soil
      col = col + 1;
      variable_names(col) = "soil_temp_" + string(n) + "_C";
      variable_values{col} = tsl_window(n, :)';
   end

   % All snow-depth variants the upstream file actually has.
   for c = variant_names
      if isfield(snd_variants, c)
         col = col + 1;
         variable_names(col) = c + "_m";
         variable_values{col} = snd_variants.(c)(idx);
      end
   end

   if has_albs
      col = col + 1;
      variable_names(col) = "albedo";
      variable_values{col} = albs(idx);
   end

   if n_cols == 0
      error('icemodel:verification:buildEsmSnowmipObservations:noChannels', ...
         'no recognised obs channels in %s', obsfile);
   end

   observations = timetable(Time, variable_values{:}, ...
      'VariableNames', cellstr(variable_names));
   observations.Time.TimeZone = 'UTC';

   metadata = struct( ...
      'sitename', sitename, ...
      'obs_file', string(obsfile), ...
      'n_rows', height(observations), ...
      'window_start', observations.Time(1), ...
      'window_end', observations.Time(end), ...
      'soil_depths_m', sdepth(:), ...
      'snow_depth_source', snow_depth_source, ...
      'swe_source', "snw_auto with snw_man fallback when available");
end
