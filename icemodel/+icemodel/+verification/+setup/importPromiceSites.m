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
   %  demo/data/eval/firn/promice/<site>/. PROMICE sites anchor the firn
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
   %    Evaluation artifacts (firn family root, demo/data/eval/firn/promice/):
   %      * evaluation.mat : PROMICE observed Data as comparison targets
   %      * reference.mat  : RACMO Data as the co-located RCM reference
   %      * per-site manifest.json fragment, rolled into the family manifest
   %
   %  Window resolution
   %    The default comparison window is the intersection of the model
   %    archives staged here. The binding constraint is the RACMO surface
   %    archive (FGRN11 2012-2015 on S03); MAR (1980-2019), MERRA-2
   %    (>=2012), and the PROMICE record all span it. With no startdate /
   %    enddate the default firn window 2012-01-01 .. 2015-12-31 is used;
   %    explicit bounds override it for every requested site.
   %
   %  Site coordinates
   %    Each site's lat/lon is read live from the v3 NetCDF metadata
   %    (readPromiceAws latitude/longitude) and converted WGS84 -> EPSG:3413
   %    with projfwd for record-keeping. The MAR / MERRA / RACMO point
   %    extractions sample at that [lat lon].
   %
   %  Name-value
   %    sites : string vector (default the curated anchor catalog,
   %        icemodel.verification.helpers.promicesiteinfo). Any v3 PROMICE
   %        station id is accepted; uncurated sites use the generic
   %        ablation recipe and are staged when their NetCDF + the RCM
   %        coverage exist.
   %    models : string vector subset of ["promice","mar","merra","racmo"]
   %        (default all four). Drop a model to stage a partial bundle.
   %    startdate, enddate : datetime / string. Explicit comparison window;
   %        both or neither. When omitted, the default firn window is used.
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
   %    demo/data/eval/firn and the standard icemodel input layout; not part
   %    of normal verification runs.
   %
   % See also: icemodel.verification.setup.importEsmSnowmip,
   %  icemodel.verification.helpers.promicesiteinfo,
   %  icemodel.verification.helpers.firnDataRoot,
   %  icemodel.forcing.buildPromiceMet, icemodel.forcing.buildPromiceData,
   %  icemodel.forcing.buildMarMet, icemodel.forcing.buildMerraMet,
   %  icemodel.forcing.buildRacmoData

   arguments
      kwargs.sites (1, :) string = ...
         string({icemodel.verification.helpers.promicesiteinfo().site})
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
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.dt_out (1, 1) string = ""
      kwargs.overwrite (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = true
   end

   sites = reshape(kwargs.sites, 1, []);
   models = reshape(kwargs.models, 1, []);

   % Resolve the comparison window. Explicit bounds must be paired (mixing
   % is an error, never a silent half-window). The default is the RACMO-
   % bound firn window (see header / Decision Log).
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
      window_start = icemodel.verification.setup.ensureUtc("2012-01-01");
      window_end = icemodel.verification.setup.ensureUtc("2015-12-31 23:00:00");
   end
   years = year(window_start):year(window_end);

   % Resolve the firn/promice family root and its manifest.
   dataset_family = "promice";
   firn_data_root = icemodel.verification.helpers.firnDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   family_root = fullfile(firn_data_root, dataset_family);
   icemodel.helpers.ensureDirExists(family_root);
   manifest_file = fullfile(family_root, "manifest.json");

   input_root = icemodel.verification.helpers.inputDataRoot( ...
      "input_data_root", kwargs.input_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   met_outdir = fullfile(input_root, 'met');
   userdata_outdir = fullfile(input_root, 'userdata');

   proj = icemodel.forcing.helpers.psnProjection();

   case_entries = {};
   skipped = struct('site', {}, 'reason', {});

   for n = 1:numel(sites)
      site = sites(n);
      alias = lower(erase(site, "_"));

      try
         % Read the station metadata (lat/lon live from the v3 NetCDF) over
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
               'ablation_recipe', promice_meta.ablation_recipe);
         end

         if ismember("mar", models)
            mar_met = icemodel.forcing.buildMarMet(point, years, ...
               source_dir=kwargs.mar_dir, modis_dir=kwargs.modis_dir, ...
               method="nearest", dt_out=kwargs.dt_out);
            mar_met = windowSubset(mar_met, window_start, window_end);
            mar_met_files = icemodel.forcing.helpers.writemet( ...
               mar_met, alias, "mar", outdir=met_outdir, naming="window");
            colocated.mar = struct( ...
               'kind', 'point_met', ...
               'met_files', relpaths(mar_met_files, met_outdir), ...
               'sample_method', 'nearest');
         end

         if ismember("merra", models)
            merra_met = icemodel.forcing.buildMerraMet(point, years, ...
               source_dir=kwargs.merra_dir, modis_dir=kwargs.modis_dir, ...
               method="nearest", dt_out=kwargs.dt_out);
            merra_met = windowSubset(merra_met, window_start, window_end);
            merra_met_files = icemodel.forcing.helpers.writemet( ...
               merra_met, alias, "merra", outdir=met_outdir, naming="window");
            colocated.merra = struct( ...
               'kind', 'point_met', ...
               'met_files', relpaths(merra_met_files, met_outdir), ...
               'sample_method', 'nearest');
         end

         racmo_data = [];
         if ismember("racmo", models)
            [racmo_data, ~] = icemodel.forcing.buildRacmoData(point, years, ...
               source_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
               method="nearest", dt="1hr");
            racmo_data = windowSubset(racmo_data, window_start, window_end);
            racmo_data_files = icemodel.forcing.helpers.writeuserdata( ...
               racmo_data, alias, "racmo", outdir=userdata_outdir);
            colocated.racmo = struct( ...
               'kind', 'point_data_smb_eval', ...
               'data_files', relpaths(racmo_data_files, userdata_outdir), ...
               'sample_method', 'nearest', ...
               'note', 'SMB/eval Data only; RACMO is not a met source.');
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
         catch
            site_name = site;
            site_note = "Uncurated PROMICE station (generic ablation recipe).";
         end

         case_values = { ...
            char(alias)
            'firn_observational'
            char(site)
            char(site_name)
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
      'observation_source', 'PROMICE v3 AWS'
      'snow_depth_method', ...
      'September-median boom-height reference minus boom height, clamped >= 0'
      'ablation_method', 'despiked transducer record, station-recipe corrected'}));

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
