function manifest = importPromiceSites(kwargs)
   %IMPORTPROMICESITES Stage PROMICE-anchored firn-evaluation cases.
   %
   %  manifest = icemodel.verification.setup.importPromiceSites()
   %  manifest = icemodel.verification.setup.importPromiceSites( ...
   %     sites=["KAN_M","KAN_L"], overwrite=true)
   %  manifest = icemodel.verification.setup.importPromiceSites( ...
   %     sites="KAN_M", startdate="2013-01-01", enddate="2015-12-31")
   %
   %  Stages a data-only observations.mat eval bundle per PROMICE automatic-
   %  weather-station site, the individual forcing/Data files, and a
   %  FORCING-AGNOSTIC case manifest. PROMICE sites anchor the firn evaluation:
   %  downstream firn-model work swaps each model's albedo into the PROMICE met
   %  and runs the firn model at the PROMICE point, comparing against the RCMs
   %  and the station observations.
   %
   %  EVAL IS FORCING-AGNOSTIC. The per-case eval target is a data-only
   %  observations.mat (same contract as ESM-SnowMIP/SUMup); the manifest records
   %  the eval contract (evaluation_file, comparison_variables) but NOT the
   %  forcing provider, so ANY forcing file may be used at the site at runtime
   %  without rewriting the eval metadata. No bundled evaluation.mat/reference.mat
   %  (forcing+obs together) is written; forcing lives in separate,
   %  runtime-discoverable met/userdata files in per-source subfolders.
   %
   %    Eval (data/eval/promice/<site>/):
   %      * observations.mat           (PROMICE obs target; buildPromiceData)
   %
   %    Forcing (data/input/met/<source>/, window naming via writemet):
   %      * PROMICE station met        (buildPromiceMet -> writemet)
   %      * MAR met at the site point  (buildMarMet     -> writemet)
   %      * MERRA-2 met at the point   (buildMerraMet   -> writemet)
   %
   %    Met-swap + reference Data (data/input/userdata/<source>/):
   %      * PROMICE Data  (met-swap source; observed channels swapped into a run -
   %                       the same buildPromiceData product as observations.mat)
   %      * RACMO Data    (eval/reference; buildRacmoData)
   %        The available RACMO 2.3p3 source files carry radiation (swd, lwd,
   %        derived albedo), turbulent fluxes (shf, lhf), precip and SMB
   %        components, but LACK the near-surface meteorological STATE variables
   %        (tair, wspd, rh, psfc) - so RACMO is staged as eval/reference Data,
   %        not a standalone met source. (RACMO in general carries these; only
   %        the available source files omit them - obtain the full set from the
   %        RACMO developers, or borrow tair/wspd/rh/psfc from MAR/MERRA/PROMICE
   %        at the point, if a RACMO-forced run is needed.)
   %
   %    Per-site forcing-agnostic manifest.json fragment, rolled into the family
   %    manifest.
   %
   %  Window resolution (per-leg, DECOUPLED; #15)
   %    PROMICE met + eval are NEVER gated by RCM coverage. A site with no
   %    MAR/MERRA/RACMO overlap still stages its FULL PROMICE met+eval record.
   %    When startdate/enddate are omitted, the PROMICE leg defaults to ALL
   %    AVAILABLE YEARS for that station (read live from the L3 record, INCLUDING
   %    partial years - e.g. KAN_U from Apr 2009). Explicit startdate/enddate
   %    override the PROMICE window. The RCM legs are independent and optional:
   %      * MAR met            : PROMICE window cap MAR years on disk
   %      * MERRA-2 met        : PROMICE window cap MERRA-2 years on disk
   %      * RACMO Data         : its OWN coverage, INDEPENDENT of the met window
   %    A leg with zero overlap is skipped-with-reason (recorded in the
   %    manifest), never fabricated. Each leg's actual staged window is recorded
   %    at colocation.<model>.window. A per-source coverage table (requested vs
   %    actual, with the missing years) is printed at the start of every run.
   %
   %  Site coordinates
   %    Each site's lat/lon is read live from the L3 NetCDF metadata
   %    (readPromiceAws latitude/longitude) and converted WGS84 -> EPSG:3413
   %    with projfwd for record-keeping. The MAR / MERRA / RACMO point
   %    extractions sample at that [lat lon].
   %
   %  Name-value
   %    sites : string vector. Default is ALL PROMICE stations found under the
   %        source product (data/verification/promice/hour/*.nc). Pass an
   %        explicit list (e.g. ["KAN_L","KAN_M","KAN_U"]) to stage a subset.
   %    models : string vector subset of ["promice","mar","merra","racmo"]
   %        (default all four). Drop a model to stage a partial bundle.
   %    startdate, enddate : datetime / string. OPTIONAL explicit PROMICE
   %        met/eval window; pass both or neither. The DEFAULT (omitted) is ALL
   %        AVAILABLE YEARS per station, read live from the L3 record (RR1) -
   %        there is no hidden 2009-2022 study window. RACMO ignores this and
   %        always uses its own on-disk coverage.
   %    output_root : base output root selecting WHICH eval tree is written.
   %        When set, eval manifest goes to <output_root>/eval and forcing/Data
   %        to <output_root>/input. The two real targets are:
   %          * COMMITTED demo fixtures : <repo>/demo/data/eval (+ .../input).
   %            Reviewed, version-controlled; only write here on purpose.
   %          * RESEARCH root (gitignored): <repo>/data/eval (+ .../input).
   %            The full research set; safe to churn.
   %        DEFAULT is SAFE: when output_root is unset the roots resolve via
   %        evaluation_data_root/input_data_root/icemodel_config_casename (the
   %        configured per-case research roots), NOT the committed demo tree -
   %        so "stage site X" never writes the committed fixtures unless an
   %        output_root/evaluation_data_root pointing at demo/data is passed.
   %    promice_dir, mar_dir, merra_dir, racmo_dir, modis_dir : raw-source
   %        directories for each model.
   %    evaluation_data_root, input_data_root, icemodel_config_casename :
   %        root resolution when output_root is unset.
   %    dt_out : optional met output timestep ("15m") for the gridded builders.
   %    overwrite : logical (default false). Refresh a REQUESTED site's own
   %        staged case folder. Other sites are never touched (see merge below).
   %    overwrite_family : logical (default false). Force a FULL rewrite of the
   %        family manifest from the requested sites alone, discarding other
   %        committed cases. The DEFAULT is MERGE.
   %    skip_missing : logical (default true). Record skip reasons and continue.
   %
   %  Incremental staging (MERGE by default)
   %    Staging one site ADDS or UPDATES only that site's case entry in the
   %    family manifest and PRESERVES every other site's committed case + files
   %    byte for byte (icemodel.verification.setup.writeFamilyManifestMerge).
   %    Re-staging the same site updates exactly its entry (idempotent). Set
   %    overwrite_family=true only to deliberately rebuild the family root.
   %
   %  Returns
   %    manifest : struct  Family manifest also written to manifest.json.
   %               manifest.skipped lists data-gated sites and the reason.
   %
   %  Role
   %    Setup/update tooling. Creates or refreshes the data-only observations.mat
   %    eval bundle, the separate forcing/Data files, and the forcing-agnostic
   %    manifest; not part of normal verification runs.
   %
   % See also: icemodel.verification.setup.importEsmSnowmip,
   %  icemodel.verification.helpers.promicesiteinfo,
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
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = true
   end

   models = reshape(kwargs.models, 1, []);

   % Site list. With no explicit sites, default to ALL PROMICE stations found
   % under the source product (the full research set).
   if isempty(kwargs.sites)
      sites = discoverStations(kwargs.promice_dir);
   else
      sites = reshape(kwargs.sites, 1, []);
   end

   % Load the AWS site metadata once (location_type) so non-curated stations
   % can resolve a first-pass surface_zone without re-reading per site.
   aws_sites = readAwsSitesMetadata(kwargs.promice_dir);

   % Resolve the explicit PROMICE met/eval window, if any. Explicit bounds must
   % be paired (mixing is an error). When OMITTED, each station defaults to its
   % FULL available record (resolved per site below from the L3 metadata), so
   % PROMICE is never gated by a hard-coded study window or by RCM coverage.
   has_startdate = ~strcmp(string(kwargs.startdate), "");
   has_enddate = ~strcmp(string(kwargs.enddate), "");
   if has_startdate ~= has_enddate
      error(['icemodel:verification:importPromiceSites:halfWindow ' ...
         'startdate and enddate must be provided together']);
   end
   explicit_window = has_startdate && has_enddate;
   if explicit_window
      req_start = icemodel.verification.setup.ensureUtc(kwargs.startdate);
      req_end = icemodel.verification.setup.ensureUtc(kwargs.enddate);
   else
      req_start = NaT;
      req_end = NaT;
   end

   % Probe each gridded source's on-disk coverage once.
   coverage = icemodel.verification.setup.promiceSourceCoverage(models, ...
      struct('mar', kwargs.mar_dir, 'merra', kwargs.merra_dir, ...
      'racmo', kwargs.racmo_dir));

   % Resolve eval/input roots.
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
         % Read station metadata over the requested window (or the full record
         % when no explicit window). This is the first gate: a missing station
         % file or empty window throws here, before any staging. PROMICE is
         % NEVER gated by RCM coverage.
         if explicit_window
            [~, aws_meta] = icemodel.forcing.readPromiceAws(site, ...
               source_dir=kwargs.promice_dir, timescale="hourly", ...
               startdate=req_start, enddate=req_end);
            promice_start = req_start;
            promice_end = req_end;
         else
            [~, aws_meta] = icemodel.forcing.readPromiceAws(site, ...
               source_dir=kwargs.promice_dir, timescale="hourly");
            % Full available record for this station, INCLUDING partial years.
            promice_start = aws_meta.window_start;
            promice_end = aws_meta.window_end;
         end

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

         % Per-leg windows decoupled from PROMICE. RCM legs use the PROMICE
         % window capped to on-disk years; RACMO uses its own coverage.
         leg = resolveLegWindows(models, coverage, promice_start, promice_end);

         % Print the requested-vs-actual coverage table for this site.
         icemodel.verification.setup.reportPromiceCoverage(coverage, ...
            [year(promice_start), year(promice_end)], legReportWindows(leg, models));

         % Prepare the eval case folder (overwrite guard lives here).
         case_root = fullfile(family_root, alias);
         icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);

         % --- Stage individual files, model by model, recording metadata. ---
         colocation = struct();
         forcing_sources = strings(0, 1);
         eval_sources = strings(0, 1);
         comparison_vars = strings(0, 1);
         obs_vars = struct();
         evaluation_file_rel = '';

         if ismember("promice", models)
            % EVAL leg FIRST: the PROMICE observations are the primary target and
            % do NOT need the met-file-required channels (e.g. lwd), so they stage
            % independently of the met leg (RR1: PROMICE eval is never gated).
            [promice_data, ~] = icemodel.forcing.buildPromiceData( ...
               site, source_dir=kwargs.promice_dir, ...
               startdate=promice_start, enddate=promice_end, frequency="daily");
            promice_data_files = icemodel.forcing.helpers.writeuserdata( ...
               promice_data, alias, "promice", outdir=userdata_outdir, ...
               naming="window");
            eval_sources(end + 1) = "promice_obs"; %#ok<AGROW>
            [comparison_vars, obs_vars] = firnComparisonContract(promice_data);

            % Forcing-AGNOSTIC eval bundle: the data-only observations.mat that
            % comparecase loads (same contract as ESM-SnowMIP/SUMup). The forcing
            % lives in separate met/userdata files, not here.
            targets = struct('format', 'timeseries', ...
               'data', promice_data, 'metadata', ...
               struct('source', 'promice_obs', 'site_id', char(site)));
            save(fullfile(case_root, 'observations.mat'), 'targets');
            evaluation_file_rel = char(fullfile(alias, 'observations.mat'));

            promice_co = struct('kind', 'station_met_and_eval');
            promice_co.data_files = ...
               icemodel.verification.setup.relpaths(promice_data_files, userdata_outdir);
            promice_co.window = windowStruct(promice_start, promice_end);

            % MET leg under its OWN guard: a station with no longwave sensor
            % (required met var lwd absent - common at firn/accumulation sites)
            % degrades to a SKIPPED met leg, never losing the eval obs above.
            try
               promice_met = icemodel.forcing.buildPromiceMet(site, ...
                  source_dir=kwargs.promice_dir, ...
                  startdate=promice_start, enddate=promice_end);
               promice_met_files = icemodel.forcing.helpers.writemet( ...
                  promice_met, alias, "promice", outdir=met_outdir, ...
                  naming="window");
               promice_co.met_files = ...
                  icemodel.verification.setup.relpaths(promice_met_files, met_outdir);
               forcing_sources(end + 1) = "promice"; %#ok<AGROW>
            catch met_err
               promice_co.met_files = strings(1, 0);
               promice_co.met_skipped_reason = string(met_err.message);
            end

            colocation.promice = promice_co;
         end

         % RCM legs are INDEPENDENT and OPTIONAL (#15 / RR1): a missing or
         % erroring MAR/MERRA/RACMO leg is recorded as a skipped leg and the
         % site's PROMICE met+eval still stages as a case. Each leg's staging
         % therefore runs under its OWN guard - the coverage probe gates the
         % expected "no on-disk coverage" case, and the try/catch additionally
         % catches a builder that throws (e.g. a source dir that vanished
         % between the probe and the build, or a read error) so the throw
         % degrades that one leg, never the whole site.
         if ismember("mar", models)
            try
               if leg.mar.staged
                  mar_met = icemodel.forcing.buildMarMet(point, leg.mar.years, ...
                     source_dir=kwargs.mar_dir, modis_dir=kwargs.modis_dir, ...
                     method="nearest", dt_out=kwargs.dt_out);
                  mar_met = windowSubset(mar_met, leg.mar.start, leg.mar.end);
                  mar_met_files = icemodel.forcing.helpers.writemet( ...
                     mar_met, alias, "mar", outdir=met_outdir, naming="window");
                  colocation.mar = struct( ...
                     'kind', 'point_met', ...
                     'met_files', icemodel.verification.setup.relpaths(mar_met_files, met_outdir), ...
                     'sample_method', 'nearest', ...
                     'window', windowStruct(leg.mar.start, leg.mar.end));
                  forcing_sources(end + 1) = "mar"; %#ok<AGROW>
               else
                  colocation.mar = skippedLeg('point_met', leg.mar.reason);
               end
            catch leg_err
               colocation.mar = skippedLeg('point_met', leg_err.message);
            end
         end

         if ismember("merra", models)
            try
               if leg.merra.staged
                  merra_met = icemodel.forcing.buildMerraMet(point, leg.merra.years, ...
                     source_dir=kwargs.merra_dir, modis_dir=kwargs.modis_dir, ...
                     method="nearest", dt_out=kwargs.dt_out);
                  merra_met = windowSubset(merra_met, leg.merra.start, leg.merra.end);
                  merra_met_files = icemodel.forcing.helpers.writemet( ...
                     merra_met, alias, "merra", outdir=met_outdir, naming="window");
                  colocation.merra = struct( ...
                     'kind', 'point_met', ...
                     'met_files', icemodel.verification.setup.relpaths(merra_met_files, met_outdir), ...
                     'sample_method', 'nearest', ...
                     'window', windowStruct(leg.merra.start, leg.merra.end));
                  forcing_sources(end + 1) = "merra"; %#ok<AGROW>
               else
                  colocation.merra = skippedLeg('point_met', leg.merra.reason);
               end
            catch leg_err
               colocation.merra = skippedLeg('point_met', leg_err.message);
            end
         end

         if ismember("racmo", models)
            try
               if leg.racmo.staged
                  % RACMO uses its OWN coverage, decoupled from the met window.
                  [racmo_data, ~] = icemodel.forcing.buildRacmoData(point, ...
                     leg.racmo.years, source_dir=kwargs.racmo_dir, ...
                     modis_dir=kwargs.modis_dir, method="nearest", dt="1hr");
                  racmo_data = windowSubset(racmo_data, leg.racmo.start, leg.racmo.end);
                  racmo_data_files = icemodel.forcing.helpers.writeuserdata( ...
                     racmo_data, alias, "racmo", outdir=userdata_outdir, ...
                     naming="window");
                  colocation.racmo = struct( ...
                     'kind', 'point_data_smb_eval', ...
                     'data_files', icemodel.verification.setup.relpaths(racmo_data_files, userdata_outdir), ...
                     'sample_method', 'nearest', ...
                     'window', windowStruct(leg.racmo.start, leg.racmo.end), ...
                     'note', 'SMB/eval Data only; RACMO is not a met source.');
                  eval_sources(end + 1) = "racmo"; %#ok<AGROW>
               else
                  colocation.racmo = skippedLeg('point_data_smb_eval', leg.racmo.reason);
               end
            catch leg_err
               colocation.racmo = skippedLeg('point_data_smb_eval', leg_err.message);
            end
         end

         % --- Forcing-agnostic manifest entry. ---
         % The eval target (observations.mat) is bundled above via
         % evaluation_file_rel; the forcing/Data sources are recorded by id only.
         anchor = siteCatalogEntry(site, aws_sites);

         case_values = { ...
            char(alias)
            'firn_observational'
            char(site)
            char(anchor.site_name)
            char(anchor.surface_zone)
            cellstr(anchor.eval_target)
            char(anchor.permafrost_zone)
            site_location
            struct('start', char(string(promice_start)), ...
            'end', char(string(promice_end)))
            evaluation_file_rel
            cellstr(forcing_sources)
            cellstr(eval_sources)
            cellstr(comparison_vars)
            obs_vars
            colocation
            'daily'
            char(anchor.note)};

         case_entries{end+1} = ...
            icemodel.verification.setup.makeFirnCaseManifestEntry(case_values); %#ok<AGROW>

      catch err
         % Only a PROMICE-leg or metadata failure reaches here and skips the
         % whole site (the RCM legs are guarded individually above and degrade
         % to skipped legs, never a skipped site - #15 / RR1). A missing
         % station file, an empty PROMICE window, or a manifest-entry error is
         % a legitimate whole-site skip.
         if ~kwargs.skip_missing
            rethrow(err)
         end
         skipped(end+1) = struct('site', site, ...
            'reason', string(err.message)); %#ok<AGROW>
         warning('icemodel:verification:importPromiceSites:siteSkipped', ...
            'skipping %s: %s', site, err.message);
      end
   end

   % Family manifest. Per-model DOIs/URLs live in each builder; the family
   % record carries the PROMICE anchor reference and the model set.
   source_doi = "";   % multi-source family; per-model provenance in builders
   source_url = "https://promice.org";
   source_version = sprintf("colocated[%s]", strjoin(models, "+"));
   retrieval_date = string(datetime('today'));

   if isempty(case_entries)
      cases = struct([]);
   else
      cases = vertcat(case_entries{:});
   end
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      dataset_family, source_doi, source_url, source_version, ...
      retrieval_date, cases);

   % Attach the data-gated sites so a refresh records exactly what was not
   % staged and why (never fabricate a case for a missing site).
   manifest.skipped = skipped;

   % MERGE by default: add/update only the requested sites' cases and preserve
   % every other committed case entry (and files) untouched. The requested-id
   % set is each requested site's compact alias, so a re-stage updates exactly
   % that case (idempotent) and a stale skip for a now-staged site clears.
   % overwrite_family forces a full rewrite of the family root.
   requested_ids = lower(erase(sites, "_"));
   manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
      manifest_file, manifest, requested_ids=requested_ids, ...
      overwrite_family=kwargs.overwrite_family);
end

%% Local helpers
function sites = discoverStations(promice_dir)
   %DISCOVERSTATIONS Full station list from the on-disk hourly NetCDF product.
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

function anchor = siteCatalogEntry(site, aws_sites)
   %SITECATALOGENTRY Resolve site_name/surface_zone/eval_target/note for a site.
   %
   % Curated and first-pass classifications both live in promicesiteinfo (the
   % single source of truth for surface_zone + eval_target). When a station is
   % not cataloged there (e.g. a brand-new L3 station absent from the AWS CSV),
   % fall back to "unknown" so the manifest still validates.
   try
      info = icemodel.verification.helpers.promicesiteinfo(site);
      anchor = struct('site_name', info.long_name, ...
         'surface_zone', info.surface_zone, ...
         'eval_target', string(info.eval_target), ...
         'permafrost_zone', info.permafrost_zone, ...
         'note', info.note);
   catch
      anchor = struct('site_name', site, ...
         'surface_zone', "unknown", ...
         'eval_target', strings(0, 1), ...
         'permafrost_zone', "unknown", ...
         'note', "Uncataloged PROMICE station; surface_zone unknown (review).");
      % Best-effort tundra recovery from the AWS CSV location_type.
      if ~isempty(aws_sites) ...
            && ismember('site_id', aws_sites.Properties.VariableNames) ...
            && ismember('location_type', aws_sites.Properties.VariableNames)
         match = string(aws_sites.site_id) == string(site);
         if any(match)
            loctype = lower(strtrim( ...
               string(aws_sites.location_type(find(match, 1)))));
            if loctype == "tundra"
               anchor.surface_zone = "tundra";
            end
         end
      end
   end
end

function leg = resolveLegWindows(models, coverage, window_start, window_end)
   %RESOLVELEGWINDOWS Decouple each gridded leg's window from the met window.
   leg = struct();

   if ismember("mar", models)
      leg.mar = capLeg(coverage.mar, window_start, window_end, "MAR");
   end
   if ismember("merra", models)
      leg.merra = capLeg(coverage.merra, window_start, window_end, "MERRA-2");
   end
   if ismember("racmo", models)
      % RACMO is decoupled: stage its full on-disk coverage, ignoring the met
      % window entirely.
      leg.racmo = ownLeg(coverage.racmo, "RACMO");
   end
end

function L = capLeg(cov, window_start, window_end, label)
   %CAPLEG Met leg: PROMICE window intersected with on-disk years.
   req_years = year(window_start):year(window_end);
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
   if kwargs.output_root ~= ""
      eval_root = fullfile(kwargs.output_root, 'eval');
      input_root = fullfile(kwargs.output_root, 'input');
      return
   end
   % If kwargs.output_root is empty, these resolve to demo/data/eval[input]
   eval_root = icemodel.verification.helpers.evaluationDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   input_root = icemodel.verification.helpers.inputDataRoot( ...
      "input_data_root", kwargs.input_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
end

function tt = windowSubset(tt, t1, t2)
   %WINDOWSUBSET Clamp a timetable to [t1, t2] on a UTC-aware axis.
   t = tt.Time;
   if isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   end
   keep = t >= t1 & t <= t2;
   tt = tt(keep, :);
end

function [comparison_vars, obs_vars] = firnComparisonContract(promice_data)
   %FIRNCOMPARISONCONTRACT Comparison axes + obs metadata from PROMICE Data.
   %
   % Records which firn comparison variables are present in the staged PROMICE
   % Data record and the observation metadata. This is the comparison CONTRACT
   % only - the observation data itself lives in the bundled observations.mat
   % eval target and the per-year userdata files, not in the manifest.
   if isempty(promice_data)
      comparison_vars = strings(0, 1);
      obs_vars = icemodel.verification.setup.metadataStruct({ ...
         'note', 'no PROMICE observations staged'});
      return
   end

   present = string(promice_data.Properties.VariableNames);
   canonical = ["ablation"; "snow_depth"; "tsfc"; ...
      "tice1"; "tice2"; "tice3"; "tice4"; ...
      "tice5"; "tice6"; "tice7"; "tice8"];
   comparison_vars = canonical(ismember(canonical, present));

   obs_vars = icemodel.verification.setup.metadataStruct({ ...
      'subsurface_temperature', 'tice1..tice8 (string thermistors)'
      'surface_ablation', 'ablation [m, positive down]'
      'snow_depth', 'snow_depth [m]'});
end
