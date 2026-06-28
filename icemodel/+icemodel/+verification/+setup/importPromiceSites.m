function manifest = importPromiceSites(source_dir, kwargs)
   %IMPORTPROMICESITES Stage PROMICE-anchored firn-evaluation cases.
   %
   %  manifest = icemodel.verification.setup.importPromiceSites()
   %  manifest = icemodel.verification.setup.importPromiceSites(source_dir, ...
   %     sites=["KAN_M","KAN_L"])
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
   %  Window resolution (per-leg, DECOUPLED)
   %    PROMICE met + eval are NEVER gated by RCM coverage. A site with no
   %    MAR/MERRA/RACMO overlap still stages its FULL PROMICE met+eval record.
   %    When startdate/enddate are omitted, the PROMICE leg defaults to ALL
   %    AVAILABLE YEARS for that station (read live from the L3 record, INCLUDING
   %    partial years - e.g. KAN_U from Apr 2009). Explicit startdate/enddate
   %    override the PROMICE window. The RCM legs are independent and optional:
   %      * MAR met            : PROMICE window cap MAR years on disk
   %      * MERRA-2 met        : PROMICE window cap MERRA-2 years on disk
   %      * RACMO Data         : PROMICE window cap RACMO years on disk;
   %                             skipped when the record has no RACMO overlap
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
   %        AVAILABLE YEARS per station, read live from the L3 record - there
   %        is no hidden study window. RACMO ignores this and
   %        always uses its own on-disk coverage.
   %    output_root : base output root selecting WHICH eval tree is written.
   %        When set, eval manifest goes to <output_root>/eval and forcing/Data
   %        to <output_root>/input. The two real targets are:
   %          * COMMITTED demo fixtures : <repo>/demo/data/eval (+ .../input).
   %            Reviewed, version-controlled; only write here on purpose.
   %          * RESEARCH root (gitignored): <repo>/data/eval (+ .../input).
   %            The full research set; safe to churn.
   %        DEFAULT (output_root unset): the roots resolve via
   %        evaluation_data_root/input_data_root/icemodel_config_casename, which
   %        for the default "test" casename point at the COMMITTED <repo>/demo/data
   %        tree (icemodel.config) - NOT a research root. To stage the gitignored
   %        research set pass output_root=<repo>/data explicitly. (Calling with no
   %        args therefore writes the committed demo tree - stage the research set
   %        only with output_root set.)
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
   %    build_forcing : logical (default false). When true (the convenience
   %        bundle), the co-located RCM forcing/Data is staged after the PROMICE
   %        observation import by delegating to
   %        icemodel.verification.setup.stageRcmForcing. When false, ONLY the
   %        PROMICE observations + station met are imported and the manifest is
   %        written; RCM forcing/Data can be built later, independently, by
   %        calling stageRcmForcing on the staged manifest (observation import is
   %        never gated on or coupled to RCM datasets).
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
      source_dir (1, 1) string = ""
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
      kwargs.build_forcing (1, 1) logical = false
   end

   models = reshape(kwargs.models, 1, []);
   kwargs.promice_dir = resolvePromiceDir(source_dir, kwargs.promice_dir);
   [kwargs.promice_dir, ~] = icemodel.verification.setup.fetchPromice( ...
      cache_dir=kwargs.promice_dir, stations=kwargs.sites, products="hour", ...
      strict=~kwargs.skip_missing, silent=true);

   % Site list. With no explicit sites, default to ALL PROMICE stations found
   % under the source product (the full research set), via the single-source-of-
   % truth auto-discovery namelist.
   if isempty(kwargs.sites)
      sites = icemodel.verification.namelists.promicesite(kwargs.promice_dir);
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

   % Preallocate the per-site state to numel(sites) and index by n; a site
   % yields at most one staged case OR one skip, so the buffers are compacted
   % at the end. The state captures everything the RCM batch passes and the
   % manifest assembly need.
   state = repmat(emptyState(), 1, numel(sites));
   alive = false(1, numel(sites));
   skipped = repmat(struct('site', "", 'reason', ""), 1, numel(sites));
   n_skipped = 0;

   % --- Pass 1: per-site PROMICE eval + met. ---
   % PROMICE met+eval is NEVER gated by RCM coverage; only a PROMICE metadata
   % or eval failure here skips the whole site. The RCM legs are staged in
   % Pass 2 over the union of staged points, then clipped per site.
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

         % Per-leg windows. RCM legs use the PROMICE window
         % capped to on-disk years. The shared resolver is the cheap fail-early
         % gate: a source with no overlap resolves staged=false here, before any
         % RCM build is attempted.
         leg = icemodel.verification.setup.resolveLegWindows( ...
            models, coverage, promice_start, promice_end);

         % Print the requested-vs-actual RCM coverage table, named by station
         % Only when forcing is actually staged - it is noise otherwise.
         if kwargs.build_forcing
            fprintf('[coverage] %s (PROMICE %d-%d):\n', site, ...
               year(promice_start), year(promice_end));
            icemodel.verification.setup.reportPromiceCoverage(coverage, ...
               [year(promice_start), year(promice_end)], ...
               legReportWindows(leg, models));
         end

         % Prepare the eval case folder (overwrite guard lives here).
         case_root = fullfile(family_root, alias);
         icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);

         % --- Stage individual files, model by model, recording metadata. ---
         colocation = struct();
         comparison_vars = strings(0, 1);
         obs_vars = struct();
         evaluation_file_rel = '';

         if ismember("promice", models)
            % EVAL leg FIRST: the PROMICE observations are the primary target and
            % do NOT need the met-file-required channels (e.g. lwd), so they stage
            % independently of the met leg.
            [promice_data, ~] = icemodel.forcing.buildPromiceData( ...
               site, source_dir=kwargs.promice_dir, ...
               startdate=promice_start, enddate=promice_end, frequency="daily");
            promice_data_files = icemodel.forcing.helpers.writeuserdata( ...
               promice_data, alias, "promice", outdir=userdata_outdir, ...
               naming="window");
            [comparison_vars, obs_vars] = firnComparisonContract(promice_data);

            % Forcing-AGNOSTIC eval bundle: the data-only observations.mat that
            % comparecase loads (same contract as ESM-SnowMIP/SUMup). The forcing
            % lives in separate met/userdata files, not here.
            targets = struct('format', 'timeseries', ...
               'data', promice_data, 'metadata', ...
               struct('source', 'promice_obs', 'site_id', char(site)));
            save(fullfile(case_root, 'observations.mat'), 'targets');
            evaluation_file_rel = char(fullfile(alias, 'observations.mat'));

            promice_co = struct('kind', 'station_met_and_eval', 'staged', true);
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
            catch met_err
               promice_co.met_files = strings(1, 0);
               promice_co.met_skipped_reason = string(met_err.message);
            end

            colocation.promice = promice_co;
         end

         % Record the alive site's state; the RCM legs are staged in Pass 2.
         state(n) = struct('site', site, 'alias', alias, 'point', point, ...
            'leg', leg, 'colocation', colocation, ...
            'site_location', site_location, ...
            'promice_start', promice_start, 'promice_end', promice_end, ...
            'comparison_vars', {comparison_vars}, 'obs_vars', obs_vars, ...
            'evaluation_file_rel', evaluation_file_rel, ...
            'anchor', siteCatalogEntry(site, aws_sites));
         alive(n) = true;

      catch err
         % Only a PROMICE-leg or metadata failure reaches here and skips the
         % whole site (the RCM legs are batched in Pass 2 and degrade to
         % skipped legs, never a skipped site. A missing station
         % file, an empty PROMICE window, or an eval-staging error is a
         % legitimate whole-site skip.
         if ~kwargs.skip_missing
            rethrow(err)
         end
         n_skipped = n_skipped + 1;
         skipped(n_skipped) = struct('site', site, ...
            'reason', string(err.message));
         warning('icemodel:verification:importPromiceSites:siteSkipped', ...
            'skipping %s: %s', site, err.message);
      end
   end

   % --- Pass 2: write the OBSERVATION manifest first, then stage
   % the co-located RCM forcing/Data as a SEPARATE, delegated step. ---
   % Observation import is complete now; persisting the manifest BEFORE any RCM
   % build means a killed/aborted forcing run never destroys the imported
   % observations. The RCM forcing/Data is delegated to
   % icemodel.verification.setup.stageRcmForcing, and observation import is
   % never gated on RCM presence.
   requested_ids = lower(erase(sites, "_"));
   manifest = assembleAndWrite(state, alive, sites, models, skipped, ...
      n_skipped, dataset_family, manifest_file, requested_ids, kwargs);

   % Convenience bundle (default): stage the co-located MAR/MERRA met+Data and
   % RACMO Data, merge the legs into each case, and rewrite the manifest. When
   % build_forcing=false the import is observations-only and stageRcmForcing can
   % be called later, independently, on the staged manifest.
   if kwargs.build_forcing
      % Persist the manifest after each RCM source so
      % a kill mid-forcing keeps every completed source's legs. The callback
      % re-assembles + MERGE-writes the manifest from the current state.
      persist = @(st) assembleAndWrite(st, alive, sites, models, skipped, ...
         n_skipped, dataset_family, manifest_file, requested_ids, kwargs);
      state = stageColocatedForcing(state, alive, models, ...
         met_outdir, userdata_outdir, kwargs, persist);
      manifest = persist(state);
   end
end

%% Local helpers
function s = emptyState()
   %EMPTYSTATE Prototype per-site staging state (preallocation seed).
   s = struct('site', "", 'alias', "", 'point', [NaN NaN], ...
      'leg', struct(), 'colocation', struct(), ...
      'site_location', struct(), 'promice_start', NaT, 'promice_end', NaT, ...
      'comparison_vars', {strings(0, 1)}, 'obs_vars', struct(), ...
      'evaluation_file_rel', '', 'anchor', struct());
end

function manifest = assembleAndWrite(state, alive, sites, models, skipped, ...
      n_skipped, dataset_family, manifest_file, requested_ids, kwargs)
   %ASSEMBLEANDWRITE Build the forcing-agnostic family manifest + MERGE-write it.
   %
   % Called twice: once with observations-only colocation (before
   % any RCM build) and once after the RCM legs are merged in. The forcing/eval
   % source lists derive from each case's colocation, so the second write simply
   % reflects the added mar/merra/racmo legs.
   case_entries = cell(1, nnz(alive));
   n_cases = 0;
   for n = 1:numel(sites)
      if ~alive(n)
         continue
      end
      s = state(n);
      [forcing_sources, eval_sources] = ...
         icemodel.verification.setup.colocationSourceLists(s.colocation);

      case_values = { ...
         char(s.alias)
         'firn_observational'
         char(s.site)
         char(s.anchor.site_name)
         char(s.anchor.surface_zone)
         cellstr(s.anchor.eval_target)
         char(s.anchor.permafrost_zone)
         s.site_location
         struct('start', char(string(s.promice_start)), ...
         'end', char(string(s.promice_end)))
         s.evaluation_file_rel
         cellstr(forcing_sources)
         cellstr(eval_sources)
         cellstr(s.comparison_vars)
         s.obs_vars
         s.colocation
         'daily'
         char(s.anchor.note)};

      n_cases = n_cases + 1;
      case_entries{n_cases} = ...
         icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);
   end
   case_entries = case_entries(1:n_cases);

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
   manifest.skipped = skipped(1:n_skipped);

   % MERGE by default: add/update only the requested sites' cases and preserve
   % every other committed case entry (and files) untouched. overwrite_family
   % forces a full rewrite of the family root.
   manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
      manifest_file, manifest, requested_ids=requested_ids, ...
      overwrite_family=kwargs.overwrite_family);
end

function state = stageColocatedForcing(state, alive, models, ...
      met_outdir, userdata_outdir, kwargs, persist)
   %STAGECOLOCATEDFORCING Delegate the RCM forcing/Data legs to stageRcmForcing.
   %
   % PROMICE is the station (already staged in Pass 1), not a gridded RCM leg, so
   % only mar/merra/racmo are staged here. Each source is staged in its OWN
   % stageRcmForcing call (still one file open per source-year), the legs are
   % merged into each site's state, and `persist` MERGE-writes the manifest after
   % each source so a partial forcing run keeps the completed sources. A
   % per-source progress line names the source and staged/skipped counts.
   % stageRcmForcing writes MAR/MERRA met+Data and RACMO Data, and degrades a
   % failing source's legs to skip-with-reason without losing the others.
   rcm_models = intersect(models, ["mar", "merra", "racmo"], "stable");
   alive_idx = find(alive);
   if isempty(alive_idx) || isempty(rcm_models)
      return
   end

   points = vertcat(state(alive_idx).point);
   for src = rcm_models
      fprintf('[staging] %s forcing for %d site(s)...\n', ...
         upper(char(src)), numel(alive_idx));

      legspec = repmat(legProto(src), 1, numel(alive_idx));
      for j = 1:numel(alive_idx)
         legspec(j).alias = state(alive_idx(j)).alias;
         legspec(j).(char(src)) = state(alive_idx(j)).leg.(char(src));
      end

      colocation = icemodel.verification.setup.stageRcmForcing(points, ...
         legspec=legspec, models=src, ...
         met_outdir=met_outdir, userdata_outdir=userdata_outdir, ...
         mar_dir=kwargs.mar_dir, merra_dir=kwargs.merra_dir, ...
         racmo_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
         method="nearest", dt_out=kwargs.dt_out);

      n_staged = 0;
      for j = 1:numel(alive_idx)
         state(alive_idx(j)).colocation = ...
            icemodel.verification.setup.mergeColocation( ...
            state(alive_idx(j)).colocation, colocation{j});
         if isfield(colocation{j}, char(src)) && colocation{j}.(char(src)).staged
            n_staged = n_staged + 1;
         end
      end
      fprintf('[staging] %s: %d staged, %d skipped\n', upper(char(src)), ...
         n_staged, numel(alive_idx) - n_staged);

      persist(state);   % Incremental manifest after this source.
   end
end

function L = legProto(models)
   %LEGPROTO Prototype legspec element (alias + one leg field per source).
   L = struct('alias', "");
   proto = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', "");
   for src = models
      L.(char(src)) = proto;
   end
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

function promice_dir = resolvePromiceDir(source_dir, promice_dir)
   %RESOLVEPROMICEDIR Resolve the positional source_dir / promice_dir alias.
   if source_dir ~= "" && promice_dir ~= "" && source_dir ~= promice_dir
      error('icemodel:verification:importPromiceSites:sourceAliasConflict', ...
         'source_dir and promice_dir refer to different PROMICE caches');
   end
   if source_dir ~= ""
      promice_dir = source_dir;
   end
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
   % tice10m (the L3 standardized 10 m-below-surface temperature) is the PRIMARY
   % subsurface comparison channel (D-QAQC-3); the raw tice1..ticeN string is
   % secondary/diagnostic. List tice10m first so it is recorded when staged.
   canonical = ["ablation"; "snow_depth"; "tsfc"; "tice10m"; ...
      "tice1"; "tice2"; "tice3"; "tice4"; ...
      "tice5"; "tice6"; "tice7"; "tice8"];
   comparison_vars = canonical(ismember(canonical, present));

   obs_vars = icemodel.verification.setup.metadataStruct({ ...
      'subsurface_temperature_primary', 'tice10m [K] (standardized 10 m below surface)'
      'subsurface_temperature_string', 'tice1..tice8 (raw string thermistors)'
      'surface_ablation', 'ablation [m, positive down]'
      'snow_depth', 'snow_depth [m]'});
end
