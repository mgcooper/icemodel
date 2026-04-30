function [forcing, metadata] = buildEsmSnowmipForcing(sitename, kwargs)
   %BUILDESMSNOWMIPFORCING Convert ESM-SnowMIP NetCDF files to icemodel-native forcing.
   %
   %  [forcing, metadata] = buildEsmSnowmipForcing(sitename)
   %  [forcing, metadata] = buildEsmSnowmipForcing(sitename, ...
   %     source_dir=..., met_kind=..., startdate=..., enddate=...)
   %
   %  Reads one ESM-SnowMIP forcing NetCDF and one observation NetCDF
   %  for the requested site, converts the variables to icemodel-native
   %  forcing names and units, and returns the result as a timetable
   %  ready for icemodel.loadmet / processmet (and also as the staged
   %  forcing.mat content consumed by importEsmSnowmip).
   %
   %  Variable mapping (ESM-SnowMIP -> icemodel-native):
   %     LWdown [W m-2]             -> lwd        [W m-2]
   %     SWdown [W m-2]             -> swd        [W m-2]
   %     Tair   [K]                 -> tair       [K]
   %     Wind   [m s-1]             -> wspd       [m s-1]
   %     Psurf  [Pa]                -> psfc       [Pa]
   %     Qair   [kg kg-1] + Psurf   -> rh         [1]
   %     Rainf + Snowf [kg m-2 s-1] -> ppt        [m s-1]
   %     obs sdepth (snd_auto/man)  -> snow_depth [m]
   %     obs albs                   -> albedo     [1]
   %
   % Notes on variable conversions:
   % Qair converts to rh via icemodel.vapor.relative_humidity_from_specific_humidity
   % albedo is filled with 0.85 over snow / 0.45 over bare ground
   % rainf + snowf = ppt, divided by ro_liq to get m/s.
   %
   %  Inputs
   %    sitename : string
   %        ESM-SnowMIP site code (e.g. "cdp"). See
   %        icemodel.verification.namelists.snowmipsite for the full
   %        catalog.
   %
   %  Name-value
   %    source_dir : string (default data/verification/snow/esm_snowmip)
   %        Directory containing met_<kind>_<sitename>_*.nc and
   %        obs_insitu_<sitename>_*.nc files.
   %    met_kind : string (default "insitu")
   %        Forcing source flavour: "insitu" (in-situ measurements) or
   %        "gswp3c" (GSWP3 bias-corrected reanalysis, 1980-2010).
   %    startdate : datetime or "" (default "")
   %        Optional start of the forcing window. When omitted, the
   %        full available range is returned.
   %    enddate : datetime or "" (default "")
   %        Optional end of the forcing window.
   %
   %  Returns
   %    forcing : timetable
   %        icemodel-native forcing in hourly UTC samples. Columns
   %        match the schema produced by the existing met builders:
   %        tair, swd, lwd, albedo, wspd, rh, psfc, ppt, snow_depth.
   %    metadata : struct
   %        Provenance: source files, row count, time bounds, the
   %        physicalConstant ro_liq used to convert mass flux to
   %        volumetric flux, and notes on the albedo policy.
   %
   %  References
   %    Menard, C. B. et al. (2019). Meteorological and evaluation
   %      datasets for snow modelling at 10 reference sites. Earth
   %      Syst. Sci. Data, 11, 865-880.
   %      https://doi.org/10.5194/essd-11-865-2019
   %    Menard, C. B. et al. (2020). Scientific and human errors in
   %      a snow model intercomparison. BAMS, 102 (2), E61-E79.
   %
   % See also: icemodel.verification.setup.importEsmSnowmip,
   %  icemodel.verification.setup.fetchEsmSnowmip,
   %  icemodel.vapor.relative_humidity_from_specific_humidity,
   %  icemodel.physicalConstant

   arguments
      sitename (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.met_kind   (1, 1) string {mustBeMember(kwargs.met_kind, ...
         ["insitu", "gswp3c"])} = "insitu"
      kwargs.startdate  = ""
      kwargs.enddate    = ""
   end

   % Resolve the source directory. Defaults to the gitignored cache
   % under <repo>/data/verification/snow/esm_snowmip/ via getpath.
   if kwargs.source_dir == ""
      source_dir = string(fullfile(icemodel.getpath('data'), ...
         'verification', 'snow', 'esm_snowmip'));
   else
      source_dir = kwargs.source_dir;
   end

   % Locate the per-site met and obs files. Match the canonical naming
   % convention used by the upstream PANGAEA bundle:
   %   met_<met_kind>_<sitename>_<start>_<end>.nc
   %   obs_insitu_<sitename>_<start>_<end>.nc
   met_pattern = sprintf('met_%s_%s_*.nc', kwargs.met_kind, sitename);
   obs_pattern = sprintf('obs_insitu_%s_*.nc', sitename);
   metfile = locateUnique(source_dir, met_pattern);
   obsfile = locateUnique(source_dir, obs_pattern);

   % --- Read the source forcing channels --------------------------------
   met_time = readNetcdfTime(metfile, 'time');
   tair     = double(ncread(metfile, 'Tair'));      % [K]
   swdown   = double(ncread(metfile, 'SWdown'));    % [W m-2]
   lwdown   = double(ncread(metfile, 'LWdown'));    % [W m-2]
   wind     = double(ncread(metfile, 'Wind'));      % [m s-1]
   psfc     = double(ncread(metfile, 'Psurf'));     % [Pa]
   qair     = double(ncread(metfile, 'Qair'));      % [kg kg-1]
   rainf    = double(ncread(metfile, 'Rainf'));     % [kg m-2 s-1]
   snowf    = double(ncread(metfile, 'Snowf'));     % [kg m-2 s-1]

   % --- Read observation channels needed for forcing -------------------
   % Only fields that feed forcing are read here (snow_depth and
   % albedo). Comparison observations are read in the obs builder.
   %
   % ESM-SnowMIP sites are heterogeneous: boreal forest sites
   % (oas, obs, ojp) report snd_gap_auto / snd_gap1_auto in the
   % canopy gap rather than snd_auto; some sites lack albs entirely.
   % readBestSnowDepth and optionalObs hide that variability so the
   % rest of the builder can stay site-agnostic.
   obs_time  = readNetcdfTime(obsfile, 'time');
   snd_auto  = readBestSnowDepth(obsfile);
   snd_man   = optionalObs(obsfile, 'snd_man', numel(obs_time));
   albs      = optionalObs(obsfile, 'albs',    numel(obs_time));

   % If met and obs grids differ, align observations to met time. This
   % handles edge cases where ESM-SnowMIP obs records start one hour
   % later than the met record (the upstream files declare different
   % "hours since" reference timestamps for some sites).
   if numel(obs_time) ~= numel(met_time) || any(obs_time ~= met_time)
      [snd_auto, snd_man, albs] = alignToMetTime(met_time, obs_time, ...
         snd_auto, snd_man, albs);
   end

   % --- Apply optional time-window subsetting --------------------------
   idx = trueMask(numel(met_time));
   if ~strcmp(string(kwargs.startdate), "")
      idx = idx & met_time >= ensureUtc(kwargs.startdate);
   end
   if ~strcmp(string(kwargs.enddate), "")
      idx = idx & met_time <= ensureUtc(kwargs.enddate);
   end
   if ~any(idx)
      error('icemodel:verification:buildEsmSnowmipForcing:emptyWindow', ...
         'no forcing samples in the requested window for site %s', sitename);
   end

   % --- Convert to icemodel-native forcing variables -------------------
   % Specific humidity -> relative humidity via the canonical icemodel
   % vapor kernel (Romps / Ambaum). Centralising this avoids site-specific
   % humidity formulas drifting across importers.
   rh = icemodel.vapor.relative_humidity_from_specific_humidity( ...
      qair(idx), psfc(idx), tair(idx));

   % Total precipitation flux: rain + snow. Mass flux [kg m-2 s-1] is
   % converted to volumetric flux [m s-1] via division by liquid water
   % density, which icemodel exposes through physicalConstant.
   ro_liq = icemodel.physicalConstant('ro_liq');
   ppt = (rainf(idx) + snowf(idx)) / ro_liq;

   % Snow depth: prefer automatic gauge, fall back to manual.
   snow_depth = preferPrimary(snd_auto(idx), snd_man(idx));

   % Albedo: build a continuous series from observed daily albedo with
   % snow / bare-ground constants where observations are missing.
   albedo = buildAlbedo(albs(idx), snow_depth);

   % --- Assemble the forcing timetable ---------------------------------
   Time = met_time(idx);
   forcing = timetable(Time, tair(idx), swdown(idx), lwdown(idx), albedo, ...
      wind(idx), rh, psfc(idx), ppt, snow_depth, ...
      'VariableNames', {'tair', 'swd', 'lwd', 'albedo', 'wspd', ...
      'rh', 'psfc', 'ppt', 'snow_depth'});
   forcing.Time.TimeZone = 'UTC';

   % Provenance metadata describing how the forcing was produced.
   metadata = struct( ...
      'sitename',         sitename, ...
      'met_kind',         kwargs.met_kind, ...
      'met_file',         string(metfile), ...
      'obs_file',         string(obsfile), ...
      'n_rows',           height(forcing), ...
      'window_start',     forcing.Time(1), ...
      'window_end',       forcing.Time(end), ...
      'ro_liq_kg_per_m3', ro_liq, ...
      'humidity_kernel',  "icemodel.vapor.relative_humidity_from_specific_humidity", ...
      'albedo_policy',    "observed daily albedo + snow/ground fallback");
end

% =====================================================================
% Local helpers
% =====================================================================

function pathname = locateUnique(source_dir, pattern)
   matches = dir(fullfile(source_dir, pattern));
   if isempty(matches)
      error('icemodel:verification:buildEsmSnowmipForcing:fileNotFound', ...
         'no file matching %s in %s', pattern, source_dir);
   end
   if numel(matches) > 1
      names = strjoin({matches.name}, ', ');
      error('icemodel:verification:buildEsmSnowmipForcing:ambiguousFile', ...
         'multiple files matching %s in %s: %s', pattern, source_dir, names);
   end
   pathname = fullfile(matches.folder, matches.name);
end

function time = readNetcdfTime(pathname, varname)
   raw = double(ncread(pathname, varname));
   units = string(ncreadatt(pathname, varname, 'units'));
   tref = datetime(extractAfter(units, 'hours since '), ...
      'InputFormat', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
   time = tref + hours(raw);
   time = time(:);
end

function values = readObs(pathname, varname)
   values = double(ncread(pathname, varname));
   values(values <= -900) = NaN;
end

function values = optionalObs(pathname, varname, ntime)
   info = ncinfo(pathname);
   if any(string({info.Variables.Name}) == varname)
      values = readObs(pathname, varname);
   else
      values = nan(ntime, 1);
   end
end

function values = readBestSnowDepth(pathname)
   %READBESTSNOWDEPTH Site-aware snow-depth selector.
   %
   % ESM-SnowMIP boreal forest sites (oas, obs, ojp) report
   % snd_gap_auto / snd_gap1_auto in the canopy gap rather than
   % snd_auto. We prefer in this order:
   %   snd_auto > snd_gap_auto > snd_gap1_auto > snd_man
   % so the resulting series is the most representative open-area
   % snow-depth time series available for the site.
   info = ncinfo(pathname);
   names = string({info.Variables.Name});
   candidates = ["snd_auto", "snd_gap_auto", "snd_gap1_auto", "snd_man"];
   for c = candidates
      if any(names == c)
         values = readObs(pathname, char(c));
         return
      end
   end
   error('icemodel:verification:buildEsmSnowmipForcing:noSnowDepth', ...
      'no usable snow-depth channel in %s', pathname);
end

function merged = preferPrimary(primary, fallback)
   merged = primary;
   replace = ~isfinite(merged) & isfinite(fallback);
   merged(replace) = fallback(replace);
end

function albedo = buildAlbedo(raw_albedo, snow_depth)
   % Filter raw observed albedo to the [0, 1] range, then linearly
   % interpolate gaps. Persistent gaps fall back to snow / bare-ground
   % constants so the staged forcing is continuous (icemodel does not
   % currently model albedo evolution from state).
   albedo = raw_albedo;
   albedo(albedo < 0 | albedo > 1) = NaN;
   albedo = fillmissing(albedo, 'linear', 'EndValues', 'nearest');
   snow_mask = isfinite(snow_depth) & snow_depth > 0.05;
   albedo(~isfinite(albedo) & snow_mask) = 0.85;
   albedo(~isfinite(albedo)) = 0.45;
   albedo = min(max(albedo, 0.05), 0.95);
end

function [a, b, c] = alignToMetTime(met_time, obs_time, a_in, b_in, c_in)
   % Fast path: same length and uniformly offset by an integer hour.
   % Otherwise, fall back to nearest-time matching with NaN fill.
   N = numel(met_time);
   a = nan(N, 1);
   b = nan(N, 1);
   c = nan(N, 1);
   [hit, locs] = ismember(met_time, obs_time);
   a(hit) = a_in(locs(hit));
   b(hit) = b_in(locs(hit));
   c(hit) = c_in(locs(hit));
end

function mask = trueMask(n)
   mask = true(n, 1);
end

function t = ensureUtc(t)
   if isstring(t) || ischar(t)
      t = datetime(t, 'TimeZone', 'UTC');
   elseif isdatetime(t) && isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   end
end
