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
   %     Rainf [kg m-2 s-1]         -> rainf      [m s-1]   (mass / ro_liq)
   %     Snowf [kg m-2 s-1]         -> snowf      [m s-1]   (mass / ro_liq)
   %     Rainf + Snowf              -> ppt        [m s-1]   (rainf + snowf)
   %     obs sdepth (snd_auto/man)  -> snow_depth [m]
   %     obs albs                   -> albedo     [1]
   %
   %  Notes on variable conversions
   %    Qair converts to rh via
   %    icemodel.vapor.relative_humidity_from_specific_humidity. rainf / snowf
   %    are stored both as their own forcing channels (so a downstream model can
   %    consume rain and snow precipitation separately when available) and
   %    summed into ppt for the legacy single-channel consumers. Mass flux [kg
   %    m-2 s-1] is divided by the canonical liquid-water density to obtain
   %    volumetric flux [m s-1]. Albedo is filled with snow / bare-ground
   %    constants where observations are missing; this is a placeholder until
   %    the diagnostic / prognostic albedo kernel lands.
   %
   %  Inputs
   %    sitename : string
   %        ESM-SnowMIP site code (e.g. "cdp"). See
   %        icemodel.verification.helpers.snowmipinfo for the full
   %        catalog and snowmipsite for the bare site-name namelist.
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
   %        tair, swd, lwd, albedo, wspd, rh, psfc, ppt, rainf, snowf,
   %        snow_depth.
   %    metadata : struct
   %        Provenance: source files, row count, time bounds, the
   %        physicalConstant ro_liq used to convert mass flux to
   %        volumetric flux, and notes on the albedo policy.
   %
   %  Role
   %    Reusable per-site forcing builder. Called by importEsmSnowmip
   %    during staging and by any future on-the-fly icemodel run that
   %    needs to consume the upstream NetCDF directly.
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
   %  icemodel.verification.setup.readNetcdfTime,
   %  icemodel.verification.setup.readObsChannel,
   %  icemodel.verification.setup.readBestSnowDepth,
   %  icemodel.vapor.relative_humidity_from_specific_humidity,
   %  icemodel.physicalConstant

   arguments
      sitename (1, 1) string ...
         {icemodel.verification.validators.mustBeSnowmipSite}
      kwargs.source_dir (1, 1) string = defaultSourceDir()
      kwargs.met_kind   (1, 1) string {mustBeMember(kwargs.met_kind, ...
         ["insitu", "gswp3c"])} = "insitu"
      kwargs.startdate  = ""
      kwargs.enddate    = ""
   end

   source_dir = kwargs.source_dir;

   % Locate the per-site met and obs files. Match the canonical naming
   % convention used by the upstream PANGAEA bundle:
   %   met_<met_kind>_<sitename>_<start>_<end>.nc
   %   obs_insitu_<sitename>_<start>_<end>.nc
   met_pattern = sprintf('met_%s_%s_*.nc', kwargs.met_kind, sitename);
   obs_pattern = sprintf('obs_insitu_%s_*.nc', sitename);
   metfile = locateUnique(source_dir, met_pattern);
   obsfile = locateUnique(source_dir, obs_pattern);

   % --- Read the source forcing channels -------------------------------
   met_time = icemodel.verification.setup.readNetcdfTime(metfile, 'time');
   tair = double(ncread(metfile, 'Tair'));      % [K]
   swdown = double(ncread(metfile, 'SWdown'));  % [W m-2]
   lwdown = double(ncread(metfile, 'LWdown'));  % [W m-2]
   wind = double(ncread(metfile, 'Wind'));      % [m s-1]
   psfc = double(ncread(metfile, 'Psurf'));     % [Pa]
   qair = double(ncread(metfile, 'Qair'));      % [kg kg-1]
   rainf = double(ncread(metfile, 'Rainf'));    % [kg m-2 s-1]
   snowf = double(ncread(metfile, 'Snowf'));    % [kg m-2 s-1]

   % --- Read observation channels needed for forcing -------------------
   % ESM-SnowMIP sites are heterogeneous: boreal forest sites
   % (oas, obs, ojp) report snd_gap_auto / snd_gap1_auto in the
   % canopy gap rather than snd_auto; some sites lack albs entirely.
   % readBestSnowDepth and readObsChannel(optional=true) hide that
   % variability so the rest of the builder can stay site-agnostic.
   obs_time  = icemodel.verification.setup.readNetcdfTime(obsfile, 'time');
   [snd_auto, snow_depth_source] = ...
      icemodel.verification.setup.readBestSnowDepth(obsfile);
   snd_man = icemodel.verification.setup.readObsChannel(obsfile, "snd_man", ...
      optional=true, ntime=numel(obs_time));
   albs = icemodel.verification.setup.readObsChannel(obsfile, "albs", ...
      optional=true, ntime=numel(obs_time));

   % If met and obs grids differ, align observations to met time. This
   % handles edge cases where ESM-SnowMIP obs records start one hour
   % later than the met record (the upstream files declare different
   % "hours since" reference timestamps for some sites).
   if numel(obs_time) ~= numel(met_time) || any(obs_time ~= met_time)
      [snd_auto, snd_man, albs] = alignToMetTime(met_time, obs_time, ...
         snd_auto, snd_man, albs);
   end

   % --- Apply optional time-window subsetting --------------------------
   idx = true(numel(met_time), 1);
   if ~strcmp(string(kwargs.startdate), "")
      idx = idx & met_time >= ...
         icemodel.verification.setup.ensureUtc(kwargs.startdate);
   end
   if ~strcmp(string(kwargs.enddate), "")
      idx = idx & met_time <= ...
         icemodel.verification.setup.ensureUtc(kwargs.enddate);
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

   % Precipitation: rain and snow stored separately and as a sum (ppt).
   % Mass flux [kg m-2 s-1] -> volumetric flux [m s-1] via division by
   % liquid-water density, which icemodel exposes through physicalConstant.
   ro_liq = icemodel.physicalConstant('ro_liq');
   rainf = rainf(idx) / ro_liq;
   snowf = snowf(idx) / ro_liq;
   ppt = rainf + snowf;

   % Snow depth: prefer automatic gauge, fall back to manual.
   snow_depth = icemodel.verification.setup.preferPrimary( ...
      snd_auto(idx), snd_man(idx));

   % Albedo: build a continuous series from observed daily albedo with
   % snow / bare-ground constants where observations are missing.
   albedo = buildAlbedo(albs(idx), snow_depth);

   % --- Assemble the forcing timetable ---------------------------------
   Time = met_time(idx);
   forcing = timetable(Time, tair(idx), swdown(idx), lwdown(idx), albedo, ...
      wind(idx), rh, psfc(idx), ppt, rainf, snowf, snow_depth, ...
      'VariableNames', {'tair', 'swd', 'lwd', 'albedo', 'wspd', ...
      'rh', 'psfc', 'ppt', 'rainf', 'snowf', 'snow_depth'});
   forcing.Time.TimeZone = 'UTC';

   % Provenance metadata describing how the forcing was produced.
   metadata = struct( ...
      'sitename', sitename, ...
      'met_kind', kwargs.met_kind, ...
      'met_file', string(metfile), ...
      'obs_file', string(obsfile), ...
      'n_rows', height(forcing), ...
      'window_start', forcing.Time(1), ...
      'window_end', forcing.Time(end), ...
      'ro_liq_kg_per_m3', ro_liq, ...
      'humidity_kernel', "icemodel.vapor.relative_humidity_from_specific_humidity", ...
      'snow_depth_source', snow_depth_source, ...
      'albedo_policy', "observed daily albedo + snow/ground fallback", ...
      'precip_policy', "rainf and snowf stored separately; ppt = rainf + snowf");
end

%% Local helpers
function pathname = defaultSourceDir()
   %DEFAULTSOURCEDIR Default ESM-SnowMIP source-cache directory.
   %
   % Mirrors fetchEsmSnowmip's default cache layout:
   %   <repo>/data/verification/snow/esm_snowmip/
   pathname = string(fullfile(icemodel.getpath('data'), ...
      'verification', 'snow', 'esm_snowmip'));
end

function pathname = locateUnique(source_dir, pattern)
   %LOCATEUNIQUE Glob-match a file in source_dir; error on miss / ambiguity.
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

function albedo = buildAlbedo(raw_albedo, snow_depth)
   %BUILDALBEDO Continuous albedo with snow / bare-ground fallback.
   %
   % Filter raw observed albedo to the [0, 1] range, then linearly
   % interpolate gaps. Persistent gaps fall back to snow / bare-ground
   % constants so the staged forcing is continuous (icemodel does not
   % currently model albedo evolution from state). The diagnostic /
   % prognostic albedo kernel (icemodel-0gt.2) will replace this
   % placeholder once it lands.
   albedo = raw_albedo;
   albedo(albedo < 0 | albedo > 1) = NaN;
   albedo = fillmissing(albedo, 'linear', 'EndValues', 'nearest');
   snow_mask = isfinite(snow_depth) & snow_depth > 0.05;
   albedo(~isfinite(albedo) & snow_mask) = 0.85;
   albedo(~isfinite(albedo)) = 0.45;
   albedo = min(max(albedo, 0.05), 0.95);
end

function [a, b, c] = alignToMetTime(met_time, obs_time, a_in, b_in, c_in)
   %ALIGNTOMETTIME Align obs samples to met times by exact-match lookup.
   N = numel(met_time);
   a = nan(N, 1);
   b = nan(N, 1);
   c = nan(N, 1);
   [hit, locs] = ismember(met_time, obs_time);
   a(hit) = a_in(locs(hit));
   b(hit) = b_in(locs(hit));
   c(hit) = c_in(locs(hit));
end
