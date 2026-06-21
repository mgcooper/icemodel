function manifest = importPromiceSites(kwargs)
   %IMPORTPROMICESITES Stage co-located multi-model firn-evaluation bundles.
   %
   %  manifest = icemodel.verification.setup.importPromiceSites()
   %  manifest = icemodel.verification.setup.importPromiceSites( ...
   %     sites=["KAN_M","KAN_L"], overwrite=true)
   %  manifest = icemodel.verification.setup.importPromiceSites( ...
   %     sites="KAN_M", startdate="2013-01-01", enddate="2015-12-31")
   %
   %  Stages a co-located, multi-model forcing/Data bundle anchored on each
   %  requested PROMICE automatic-weather-station site, under
   %  demo/data/eval/promice/<site>/. PROMICE sites anchor the firn
   %  evaluation: downstream firn-model work swaps each model's albedo into
   %  the PROMICE met and runs the firn model at the PROMICE point, comparing
   %  against the RCMs and the station observations. This driver therefore
   %  produces, per site:
   %
   %    Forcing (icemodel input layout, data/input/met/, window naming):
   %      * PROMICE station met       (buildPromiceMet  -> writemet)
   %      * MAR met at the site point  (buildMarMet     -> writemet)
   %      * MERRA-2 met at the point   (buildMerraMet   -> writemet)
   %
   %    Data / met-swap userdata (data/input/userdata/, per-year naming):
   %      * PROMICE evaluation Data    (buildPromiceData -> writeuserdata)
   %      * RACMO Data at the point    (buildRacmoData   -> writeuserdata)
   %        (RACMO carries no met channels and is never a met source.)
   %
   %    Evaluation artifacts (promice family root, demo/data/eval/promice/):
   %      * evaluation.mat : PROMICE observed Data as comparison targets
   %      * reference.mat  : RACMO Data as the co-located RCM reference
   %      * per-site manifest.json fragment, rolled into the family manifest
   %
   %  Window resolution (per-leg, DECOUPLED)
   %    The requested study window defaults to 2009-01-01 .. 2022-12-31.
   %    Each leg is staged over the requested window intersected with that
   %    source's ON-DISK availability, probed live by
   %    icemodel.verification.setup.promiceSourceCoverage:
   %      * PROMICE met + eval : requested window (full station record spans it)
   %      * MAR met            : requested window cap MAR years on disk
   %      * MERRA-2 met        : requested window cap MERRA-2 years on disk
   %      * RACMO Data         : its OWN coverage (FGRN11 surface ~2012-2015,
   %                             subsurface ~2012-2018), INDEPENDENT of the met
   %                             window. RACMO is never clamped to the met span.
   %    A leg with zero overlap is skipped-with-reason (recorded in the
   %    manifest), never fabricated. Each leg's actual staged window is
   %    recorded at colocated_forcing.<model>.window. A per-source coverage
   %    table (requested vs actual, with the missing years) is printed at the
   %    start of every run. Explicit startdate/enddate override the default
   %    study window for the met/eval legs; RACMO keeps its own coverage.
   %
   %  Site coordinates
   %    Each site's lat/lon is read live from the v3 NetCDF metadata
   %    (readPromiceAws latitude/longitude) and converted WGS84 -> EPSG:3413
   %    with projfwd for record-keeping. The MAR / MERRA / RACMO point
   %    extractions sample at that [lat lon].
   %
   %  Name-value
   %    sites : string vector. Default is ALL PROMICE stations found under
   %        the source product (data/verification/promice/hour/*.nc). Pass an
   %        explicit list (e.g. ["KAN_L","KAN_M","KAN_U"]) to stage a subset.
   %        Any station id is accepted; uncurated sites use the generic
   %        ablation recipe and are staged when their NetCDF exists.
   %    models : string vector subset of ["promice","mar","merra","racmo"]
   %        (default all four). Drop a model to stage a partial bundle.
   %    startdate, enddate : datetime / string. Explicit study window for the
   %        met/eval legs; both or neither. When omitted, the default study
   %        window 2009-01-01 .. 2022-12-31 is used. RACMO ignores this and
   %        always uses its own on-disk coverage.
   %    output_root : base output root. When set, the committed-vs-research
   %        routing is explicit: evaluation artifacts go to <output_root>/eval
   %        and forcing/Data to <output_root>/input. Use
   %        icemodel.internal.fullpath("demo","data") for the committed CI
   %        fixtures (demo/data/eval + demo/data/input) and
   %        icemodel.internal.fullpath("data") for the gitignored research set
   %        (data/eval + data/input). When unset, evaluation_data_root /
   %        input_data_root (or the active config) resolve the roots.
   %    promice_dir, mar_dir, merra_dir, racmo_dir, modis_dir : raw-source
   %        directories for each model (default the gitignored caches /
   %        the S03 reference layouts encoded in each builder).
   %    evaluation_data_root, icemodel_config_casename : eval-root resolution
   %        (default the active config / "test").
   %    input_data_root : base input-data root holding met/ and userdata/
   %        (default config-derived). Forcing is staged under <root>/met/
   %        and met-swap Data under <root>/userdata/.
   %    dt_out : optional met output timestep ("15m") for the gridded met
   %        builders (PROMICE met stays at its native hourly axis).
   %    overwrite : logical (default false). Refresh staged artifacts.
   %    skip_missing : logical (default true). When a site cannot be staged
   %        (missing NetCDF or RCM coverage), record the reason and continue
   %        rather than erroring; set false to fail on the first gap.
   %
   %  Returns
   %    manifest : struct  Family manifest also written to manifest.json.
   %               manifest.skipped lists data-gated sites and the reason.
   %
   %  Role
   %    Setup/update tooling. Creates or refreshes staged data under
   %    demo/data/eval/promice and the standard icemodel input layout; not part
   %    of normal verification runs.
   %
   % See also: icemodel.verification.setup.importEsmSnowmip,
   %  icemodel.verification.helpers.promicesiteinfo,
   %  icemodel.verification.helpers.evaluationDataRoot,
   %  icemodel.forcing.buildPromiceMet, icemodel.forcing.buildPromiceData,
   %  icemodel.forcing.buildMarMet, icemodel.forcing.buildMerraMet,
   %  icemodel.forcing.buildRacmoData

   arguments
      kwargs.sites (1, :) string = strings(1, 0)
      kwargs.models (1, :) string {mustBeMember(kwargs.models, ...
         ["promice", "mar", "merra", "racmo"])} = ...
         ["promice", "mar", "merra", "racmo"]
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.promice_dir (1, 1) string = ""
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.dt_out (1, 1) string = ""
      kwargs.overwrite (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = true
   end

   models = reshape(kwargs.models, 1, []);

   % Site list. With no explicit sites, default to ALL PROMICE stations found
   % under the source product (the full research set). An explicit list stages
   % a subset (e.g. the committed KAN fixtures).
   if isempty(kwargs.sites)
      sites = discoverStations(kwargs.promice_dir);
   else
      sites = reshape(kwargs.sites, 1, []);
   end

   % Load the AWS site metadata once (location_type / project) so non-KAN
   % stations can resolve a surface_zone without re-reading per site.
   aws_sites = readAwsSitesMetadata(kwargs.promice_dir);

   % Resolve the study window for the met/eval legs. Explicit bounds must be
   % paired (mixing is an error, never a silent half-window). The default is
   % the 2009-2022 study window; RACMO is decoupled and uses its own coverage.
   has_startdate = ~strcmp(string(kwargs.startdate), "");
   has_enddate = ~strcmp(string(kwargs.enddate), "");
   if has_startdate ~= has_enddate
      error(['icemodel:verification:importPromiceSites:halfWindow ' ...
         'startdate and enddate must be provided together']);
   end
   if has_startdate && has_enddate
      window_start = icemodel.verification.setup.ensureUtc(kwargs.startdate);
      window_end = icemodel.verification.setup.ensureUtc(kwargs.enddate);
   else
      window_start = icemodel.verification.setup.ensureUtc("2009-01-01");
      window_end = icemodel.verification.setup.ensureUtc("2022-12-31 23:00:00");
   end

   % Probe each gridded source's on-disk coverage and DECOUPLE the per-leg
   % windows: MAR/MERRA met use the requested window capped to their on-disk
   % years; RACMO Data uses its own coverage independent of the met window.
   coverage = icemodel.verification.setup.promiceSourceCoverage(models, ...
      struct('mar', kwargs.mar_dir, 'merra', kwargs.merra_dir, ...
      'racmo', kwargs.racmo_dir));

   % Per-leg windows and the years each gridded builder is allowed to read.
   leg = resolveLegWindows(models, coverage, window_start, window_end);

   % Print the requested-vs-actual coverage table (honest gaps reporting).
   icemodel.verification.setup.reportPromiceCoverage(coverage, ...
      [year(window_start), year(window_end)], legReportWindows(leg, models));

   % Resolve the eval/input roots. output_root, when set, routes eval to
   % <output_root>/eval and forcing/Data to <output_root>/input (the explicit
   % committed-vs-research switch); otherwise the per-root kwargs / config win.
   [evaluation_data_root, input_root] = resolveRoots(kwargs);

   dataset_family = "promice";
   family_root = fullfile(evaluation_data_root, dataset_family);
   icemodel.helpers.ensureDirExists(family_root);
   manifest_file = fullfile(family_root, "manifest.json");

   met_outdir = fullfile(input_root, 'met');
   userdata_outdir = fullfile(input_root, 'userdata');
   icemodel.helpers.ensureDirExists(met_outdir);
   icemodel.helpers.ensureDirExists(userdata_outdir);

   proj = icemodel.forcing.helpers.psnProjection();

   case_entries = {};
   skipped = struct('site', {}, 'reason', {});

   for n = 1:numel(sites)
      site = sites(n);
      alias = lower(erase(site, "_"));

      try
         % Read the station metadata (lat/lon live from the L3 NetCDF) over
         % the requested window. This is also the first gate: a missing
         % station file or empty window throws here, before any staging.
         [~, aws_meta] = icemodel.forcing.readPromiceAws(site, ...
            source_dir=kwargs.promice_dir, timescale="hourly", ...
            startdate=window_start, enddate=window_end);

         lat = aws_meta.lat;
         lon = aws_meta.lon;
         [x3413, y3413] = projfwd(proj, lat, lon);
         point = [lat, lon];

         site_location = struct( ...
            'lat_wgs84', lat, ...
            'lon_wgs84', lon, ...
            'x_epsg3413', x3413, ...
            'y_epsg3413', y3413, ...
            'elev_m', aws_meta.elev);

         % Prepare the eval case folder (overwrite guard lives here).
         case_root = fullfile(family_root, alias);
         icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);

         % --- Build + stage the co-located bundle, model by model. ---
         colocated = struct();
         promice_data = [];

         if ismember("promice", models)
            % PROMICE met + eval over the full requested study window.
            promice_met = icemodel.forcing.buildPromiceMet(site, ...
               source_dir=kwargs.promice_dir, ...
               startdate=window_start, enddate=window_end);
            promice_met_files = icemodel.forcing.helpers.writemet( ...
               promice_met, alias, "promice", outdir=met_outdir, ...
               naming="window");

            [promice_data, promice_meta] = icemodel.forcing.buildPromiceData( ...
               site, source_dir=kwargs.promice_dir, ...
               startdate=window_start, enddate=window_end, frequency="daily");
            promice_data_files = icemodel.forcing.helpers.writeuserdata( ...
               promice_data, alias, "promice", outdir=userdata_outdir);

            colocated.promice = struct( ...
               'kind', 'station_met_and_eval', ...
               'met_files', relpaths(promice_met_files, met_outdir), ...
               'data_files', relpaths(promice_data_files, userdata_outdir), ...
               'ablation_source', promice_meta.ablation_source, ...
               'window', windowStruct(window_start, window_end));
         end

         if ismember("mar", models)
            if leg.mar.staged
               mar_met = icemodel.forcing.buildMarMet(point, leg.mar.years, ...
                  source_dir=kwargs.mar_dir, modis_dir=kwargs.modis_dir, ...
                  method="nearest", dt_out=kwargs.dt_out);
               mar_met = windowSubset(mar_met, leg.mar.start, leg.mar.end);
               mar_met_files = icemodel.forcing.helpers.writemet( ...
                  mar_met, alias, "mar", outdir=met_outdir, naming="window");
               colocated.mar = struct( ...
                  'kind', 'point_met', ...
                  'met_files', relpaths(mar_met_files, met_outdir), ...
                  'sample_method', 'nearest', ...
                  'window', windowStruct(leg.mar.start, leg.mar.end));
            else
               colocated.mar = skippedLeg('point_met', leg.mar.reason);
            end
         end

         if ismember("merra", models)
            if leg.merra.staged
               merra_met = icemodel.forcing.buildMerraMet(point, leg.merra.years, ...
                  source_dir=kwargs.merra_dir, modis_dir=kwargs.modis_dir, ...
                  method="nearest", dt_out=kwargs.dt_out);
               merra_met = windowSubset(merra_met, leg.merra.start, leg.merra.end);
               merra_met_files = icemodel.forcing.helpers.writemet( ...
                  merra_met, alias, "merra", outdir=met_outdir, naming="window");
               colocated.merra = struct( ...
                  'kind', 'point_met', ...
                  'met_files', relpaths(merra_met_files, met_outdir), ...
                  'sample_method', 'nearest', ...
                  'window', windowStruct(leg.merra.start, leg.merra.end));
            else
               colocated.merra = skippedLeg('point_met', leg.merra.reason);
            end
         end

         racmo_data = [];
         if ismember("racmo", models)
            if leg.racmo.staged
               % RACMO uses its OWN coverage, decoupled from the met window.
               [racmo_data, ~] = icemodel.forcing.buildRacmoData(point, ...
                  leg.racmo.years, source_dir=kwargs.racmo_dir, ...
                  modis_dir=kwargs.modis_dir, method="nearest", dt="1hr");
               racmo_data = windowSubset(racmo_data, leg.racmo.start, leg.racmo.end);
               racmo_data_files = icemodel.forcing.helpers.writeuserdata( ...
                  racmo_data, alias, "racmo", outdir=userdata_outdir);
               colocated.racmo = struct( ...
                  'kind', 'point_data_smb_eval', ...
                  'data_files', relpaths(racmo_data_files, userdata_outdir), ...
                  'sample_method', 'nearest', ...
                  'window', windowStruct(leg.racmo.start, leg.racmo.end), ...
                  'note', 'SMB/eval Data only; RACMO is not a met source.');
            else
               colocated.racmo = skippedLeg('point_data_smb_eval', leg.racmo.reason);
            end
         end

         % --- Evaluation + reference artifacts for the verification suite. ---
         % Targets are the PROMICE observed Data (firn evaluation anchor).
         % The reference is the co-located RACMO Data (the RCM the firn model
         % is compared against). Both fall back gracefully when their model
         % was excluded from the requested set.
         [targets, obs_vars, comparison_vars] = ...
            evaluationTargets(promice_data);
         reference = referenceArtifact(racmo_data);

         evaluation_file = fullfile(case_root, "evaluation.mat");
         reference_file = fullfile(case_root, "reference.mat");
         save(evaluation_file, 'targets');
         save(reference_file, 'reference');

         % --- Manifest entry. ---
         try
            anchor = icemodel.verification.helpers.promicesiteinfo(site);
            site_name = anchor.long_name;
            site_note = anchor.note;
            % KAN sites: surface_zone single-sourced from the curated catalog
            % so the manifest never drifts from promicesiteinfo (kanl=
            % lower_ablation, kanm=upper_ablation, kanu=lower_percolation).
            surface_zone = anchor.zone;
         catch
            % Non-curated stations: surface_zone from the AWS site metadata
            % location_type field (ice sheet / tundra / ...), normalized to the
            % canonical vocabulary; "unknown" when the CSV carries no usable
            % zone for this station.
            site_name = site;
            site_note = "Uncurated PROMICE station.";
            surface_zone = surfaceZoneFromAws(aws_sites, site);
         end

         case_values = { ...
            char(alias)
            'firn_observational'
            char(site)
            char(site_name)
            char(surface_zone)
            site_location
            char(fullfile(alias, "evaluation.mat"))
            char(fullfile(alias, "reference.mat"))
            'daily'
            struct('start', char(string(window_start)), ...
            'end', char(string(window_end)))
            cellstr(comparison_vars)
            obs_vars
            colocated
            char(site_note)};

         case_entries{end+1} = ...
            icemodel.verification.setup.makeFirnCaseManifestEntry(case_values); %#ok<AGROW>

      catch err
         if ~kwargs.skip_missing
            rethrow(err)
         end
         skipped(end+1) = struct('site', site, ...
            'reason', string(err.message)); %#ok<AGROW>
         warning('icemodel:verification:importPromiceSites:siteSkipped', ...
            'skipping %s: %s', site, err.message);
      end
   end

   % Family manifest. Provenance points at the multi-model raw sources;
   % the per-model DOIs/URLs live in each builder, so the family record
   % carries the PROMICE anchor reference and the model set.
   source_doi = "";   % multi-source family; per-model provenance in builders
   source_url = "https://promice.org";
   source_version = sprintf("colocated[%s]", strjoin(models, "+"));
   retrieval_date = string(datetime('today'));

   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      dataset_family, source_doi, source_url, source_version, ...
      retrieval_date, vertcat(case_entries{:}));

   % Attach the data-gated sites so a refresh records exactly what was not
   % staged and why (never fabricate a case for a missing site).
   manifest.skipped = skipped;

   icemodel.verification.setup.writeManifest(manifest_file, manifest);
end

%% Local helpers
function sites = discoverStations(promice_dir)
   %DISCOVERSTATIONS Full station list from the on-disk hourly NetCDF product.
   %
   % Resolves every <STATION>_hour.nc under the PROMICE product so a no-site
   % call stages the full research set. Mirrors readPromiceAws source-dir
   % resolution (the directory holding hour/ or the directory itself).
   source_dir = promice_dir;
   if source_dir == ""
      source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
         'verification', 'promice'));
   end
   if isfolder(fullfile(source_dir, 'hour'))
      source_dir = fullfile(source_dir, 'hour');
   end
   files = dir(fullfile(source_dir, '*_hour.nc'));
   if isempty(files)
      error('icemodel:verification:importPromiceSites:noStations', ...
         'no <STATION>_hour.nc files found under %s', source_dir)
   end
   sites = reshape(string(erase({files.name}, "_hour.nc")), 1, []);
end

function aws = readAwsSitesMetadata(promice_dir)
   %READAWSSITESMETADATA Read AWS_sites_metadata.csv (or empty on absence).
   source_dir = promice_dir;
   if source_dir == ""
      source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
         'verification', 'promice'));
   end
   if isfolder(fullfile(source_dir, 'hour')) ...
         && ~isfile(fullfile(source_dir, 'AWS_sites_metadata.csv'))
      % readPromiceAws accepts the hour/ dir; the CSV sits one level up.
   end
   csv = fullfile(source_dir, 'AWS_sites_metadata.csv');
   if ~isfile(csv)
      csv = fullfile(fileparts(char(source_dir)), 'AWS_sites_metadata.csv');
   end
   if ~isfile(csv)
      aws = table();
      return
   end
   aws = readtable(csv, 'TextType', 'string');
end

function zone = surfaceZoneFromAws(aws_sites, site)
   %SURFACEZONEFROMAWS Map the AWS location_type to a canonical surface_zone.
   %
   % The AWS_sites_metadata.csv location_type field is coarse (e.g. "ice
   % sheet", "tundra"); it does NOT resolve a firn zone. "tundra" maps to the
   % canonical "tundra"; everything else (incl. "ice sheet", which carries no
   % elevation-zone detail) is recorded as "unknown" rather than guessed.
   zone = "unknown";
   if isempty(aws_sites) || ~ismember('site_id', aws_sites.Properties.VariableNames)
      return
   end
   match = string(aws_sites.site_id) == string(site);
   if ~any(match) || ~ismember('location_type', aws_sites.Properties.VariableNames)
      return
   end
   loctype = lower(strtrim(string(aws_sites.location_type(find(match, 1)))));
   switch loctype
      case "tundra"
         zone = "tundra";
      otherwise
         zone = "unknown";
   end
end

function leg = resolveLegWindows(models, coverage, window_start, window_end)
   %RESOLVELEGWINDOWS Decouple each gridded leg's window from the met window.
   %
   % MAR/MERRA: requested window intersected with their on-disk years.
   % RACMO:     its OWN coverage, INDEPENDENT of the requested window.
   % A leg with zero overlap is marked not-staged with a reason.
   leg = struct();
   req_years = year(window_start):year(window_end);

   if ismember("mar", models)
      leg.mar = capLeg(coverage.mar, req_years, window_start, window_end, "MAR");
   end
   if ismember("merra", models)
      leg.merra = capLeg(coverage.merra, req_years, window_start, window_end, ...
         "MERRA-2");
   end
   if ismember("racmo", models)
      % RACMO is decoupled: stage its full on-disk coverage, ignoring the
      % requested met window entirely.
      leg.racmo = ownLeg(coverage.racmo, "RACMO");
   end
end

function L = capLeg(cov, req_years, window_start, window_end, label)
   %CAPLEG Met leg: requested window intersected with on-disk years.
   if isempty(cov.years)
      L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
         'reason', sprintf('%s absent (%s)', label, cov.reason));
      return
   end
   years = intersect(req_years, cov.years);
   if isempty(years)
      L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
         'reason', sprintf('%s on-disk %d-%d has no overlap with requested %d-%d', ...
         label, cov.year_min, cov.year_max, req_years(1), req_years(end)));
      return
   end
   % Clamp the window to the contiguous span the builders can read. The met
   % builders read whole calendar years, so the staged window starts no
   % earlier than max(requested start, first available year) and ends no
   % later than min(requested end, last available year).
   y0 = max(year(window_start), min(years));
   y1 = min(year(window_end), max(years));
   t1 = max(window_start, icemodel.verification.setup.ensureUtc( ...
      sprintf('%d-01-01', y0)));
   t2 = min(window_end, icemodel.verification.setup.ensureUtc( ...
      sprintf('%d-12-31 23:00:00', y1)));
   L = struct('staged', true, 'years', years(years >= y0 & years <= y1), ...
      'start', t1, 'end', t2, 'reason', "");
end

function L = ownLeg(cov, label)
   %OWNLEG RACMO leg: its full on-disk coverage, decoupled from the met window.
   if isempty(cov.years)
      L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
         'reason', sprintf('%s absent (%s)', label, cov.reason));
      return
   end
   t1 = icemodel.verification.setup.ensureUtc( ...
      sprintf('%d-01-01', cov.year_min));
   t2 = icemodel.verification.setup.ensureUtc( ...
      sprintf('%d-12-31 23:00:00', cov.year_max));
   L = struct('staged', true, 'years', cov.years, 'start', t1, 'end', t2, ...
      'reason', "");
end

function w = legReportWindows(leg, models)
   %LEGREPORTWINDOWS Flatten the per-leg windows for the coverage report.
   w = struct();
   for m = ["promice", "mar", "merra", "racmo"]
      if m == "promice" || ~ismember(m, models)
         continue
      end
      L = leg.(m);
      if L.staged
         w.(char(m)) = struct('start', L.start, 'end', L.end);
      else
         w.(char(m)) = sprintf('skipped: %s', L.reason);
      end
   end
end

function w = windowStruct(t1, t2)
   %WINDOWSTRUCT JSON-friendly window record for the manifest.
   w = struct('start', char(string(t1)), 'end', char(string(t2)));
end

function leg = skippedLeg(kind, reason)
   %SKIPPEDLEG Manifest entry for a leg with no on-disk coverage.
   leg = struct('kind', kind, 'staged', false, ...
      'reason', char(string(reason)));
end

function [eval_root, input_root] = resolveRoots(kwargs)
   %RESOLVEROOTS Resolve eval/input roots, honoring output_root when set.
   %
   % output_root is the explicit committed-vs-research switch: eval ->
   % <output_root>/eval, input -> <output_root>/input. When unset, the
   % per-root kwargs (or the active config) resolve each root independently.
   if kwargs.output_root ~= ""
      eval_root = fullfile(kwargs.output_root, 'eval');
      input_root = fullfile(kwargs.output_root, 'input');
      return
   end
   eval_root = icemodel.verification.helpers.evaluationDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   input_root = icemodel.verification.helpers.inputDataRoot( ...
      "input_data_root", kwargs.input_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
end

function tt = windowSubset(tt, t1, t2)
   %WINDOWSUBSET Clamp a timetable to [t1, t2] on a UTC-aware axis.
   %
   % The gridded builders extract whole calendar years; this trims the
   % staged bundle to the exact comparison window so every model in the
   % bundle shares the same span.
   t = tt.Time;
   if isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   end
   keep = t >= t1 & t <= t2;
   tt = tt(keep, :);
end

function rel = relpaths(filenames, base)
   %RELPATHS Reduce absolute staged paths to base-relative names for JSON.
   %
   % Returns a string row vector. A string array (not a cell) keeps the
   % enclosing provenance struct scalar - passing a multi-element cell to
   % struct() would expand it into a struct array, which jsonencode then
   % renders as a list of near-duplicate objects instead of one object
   % with a list-valued field.
   filenames = string(filenames);
   base = string(base);
   rel = erase(filenames, base + filesep);
   rel = reshape(rel, 1, []);
end

function [targets, obs_vars, comparison_vars] = evaluationTargets(promice_data)
   %EVALUATIONTARGETS Build the firn evaluation targets from PROMICE Data.
   %
   % The PROMICE observed Data is the firn evaluation anchor. The firn
   % comparison variables present in the staged record are recorded; the
   % obs metadata stays JSON-friendly. When PROMICE was excluded from the
   % requested model set the targets are an explicit empty placeholder.
   if isempty(promice_data)
      targets = struct('format', 'timeseries', 'data', [], ...
         'metadata', icemodel.verification.setup.metadataStruct({ ...
         'note', 'PROMICE excluded from requested model set'}));
      obs_vars = icemodel.verification.setup.metadataStruct({ ...
         'note', 'no PROMICE observations staged'});
      comparison_vars = strings(0, 1);
      return
   end

   % Canonical firn comparison axes present in the PROMICE Data record.
   present = string(promice_data.Properties.VariableNames);
   canonical = ["ablation"; "snow_depth"; "tsfc"; ...
      "tice1"; "tice2"; "tice3"; "tice4"; ...
      "tice5"; "tice6"; "tice7"; "tice8"];
   comparison_vars = canonical(ismember(canonical, present));

   targets = struct('format', 'timeseries', 'data', promice_data, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'observation_source', 'pypromice L3 AWS'
      'snow_depth_method', 'L3 snow_height channel (read, not derived)'
      'ablation_method', ...
      'L3 ice-surface height (z_ice_surf / z_surf_combined), read not derived'}));

   obs_vars = icemodel.verification.setup.metadataStruct({ ...
      'subsurface_temperature', 'tice1..tice8 (string thermistors)'
      'surface_ablation', 'ablation [m, positive down]'
      'snow_depth', 'snow_depth [m]'});
end

function reference = referenceArtifact(racmo_data)
   %REFERENCEARTIFACT Build the RCM reference from co-located RACMO Data.
   %
   % RACMO is the co-located RCM the firn model is compared against. When
   % RACMO was excluded the reference is an explicit placeholder rather
   % than fabricated data.
   if isempty(racmo_data)
      reference = struct('format', 'timeseries', 'data', [], ...
         'metadata', icemodel.verification.setup.metadataStruct({ ...
         'reference_kind', 'no_reference_staged'
         'note', 'RACMO excluded from requested model set'}));
      return
   end
   reference = struct('format', 'timeseries', 'data', racmo_data, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'reference_kind', 'colocated_rcm'
      'reference_source', 'RACMO2.3p3 FGRN11 (point extraction)'
      'note', 'SMB/eval Data only; RACMO carries no met channels.'}));
end
