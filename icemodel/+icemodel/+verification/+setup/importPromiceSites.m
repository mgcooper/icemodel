function manifest = importPromiceSites(source_dir, kwargs)
   %IMPORTPROMICESITES Stage PROMICE-anchored firn-evaluation cases.
   %
   %  manifest = icemodel.verification.setup.importPromiceSites()
   %  manifest = icemodel.verification.setup.importPromiceSites(source_dir, ...
   %     case_ids=["KAN_M","KAN_L"])
   %  manifest = icemodel.verification.setup.importPromiceSites( ...
   %     case_ids=["KAN_M","KAN_L"], overwrite=true)
   %  manifest = icemodel.verification.setup.importPromiceSites( ...
   %     case_ids="KAN_M", startdate="2013-01-01", enddate="2015-12-31")
   %
   %  Role
   %    Setup/update tooling. Creates or refreshes the data-only observations.mat
   %    eval bundle, the separate forcing/Data files, and the forcing-agnostic
   %    manifest; not part of normal verification runs.
   %    forcing_sources selects runtime sources requested by the current call.
   %    Ordinary calls preserve omitted existing legs; overwrite_family=true
   %    deliberately replaces the whole family state.
   %    build_observations=false is a guarded non-dry fast path: requested cases
   %    must already exist in the target manifest, whose observation entry is
   %    reused while selected forcing is attached.
   %
   %  Default roots
   %    source_dir="" reads <repo>/data/verification/promice/hour. With no output_root,
   %    observations go to <repo>/data/eval/promice/<case_id>/observations.mat and
   %    native met/userdata go to <repo>/data/input/{met,userdata}/promice/.
   %    Explicit source_dir, output_root, evaluation_data_root, and
   %    input_data_root overrides are honored as-is.
   %
   %  Met and userdata
   %    Model met defaults to dt_out="15m"; pass dt_out="" for native cadence.
   %    Data/userdata defaults to hourly at the shared writer boundary.
   %    Native met schema completion is fixed at the importer boundary: absent
   %    required channels are retained as NaN placeholders.
   %    Call buildPromiceMet directly for strict source-schema validation.
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
   %  Window selection
   %    startdate/enddate optionally clamp each requested PROMICE hourly record.
   %    This is intended for short preview staging; omit both for full
   %    production artifacts.
   %
   %  Per-leg resolution (DECOUPLED)
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
   %    case_ids : string vector. Default is ALL PROMICE stations found under the
   %        source product (data/verification/promice/hour/*.nc). Pass an
   %        explicit list (e.g. ["KAN_L","KAN_M","KAN_U"]) to stage a subset.
   %    forcing_sources : string vector subset of
   %        ["promice","mar","merra","racmo"] (default "promice"). PROMICE
   %        observations are always the case definition when build_observations
   %        is true; forcing_sources selects only runtime met/userdata artifacts.
   %        It is a patch selector, not the complete desired source state: an
   %        existing site's omitted legs remain unchanged during ordinary merge
   %        updates and are removed only by explicit family replacement.
   %    startdate, enddate : datetime / string. OPTIONAL explicit PROMICE
   %        met/eval window; pass both or neither. The DEFAULT (omitted) is ALL
   %        AVAILABLE YEARS per station, read live from the L3 record - there
   %        is no hidden study window. RACMO ignores this and
   %        always uses its own on-disk coverage.
   %    output_root : base output root selecting WHICH eval tree is written.
   %        When set, eval manifest goes to <output_root>/eval and forcing/Data
   %        to <output_root>/input. The normal verification target is
   %        <repo>/data/eval (+ .../input). <repo>/demo/data is reserved for
   %        deliberate legacy/demo comparisons or fixture-maintenance work.
   %        DEFAULT (output_root unset): the roots resolve to the gitignored
   %        top-level <repo>/data tree. Pass output_root=<repo>/demo/data only
   %        when deliberately refreshing legacy/demo artifacts.
   %    promice_dir, mar_dir, merra_dir, racmo_dir, modis_dir : raw-source
   %        directories for each model.
   %    evaluation_data_root, input_data_root, icemodel_config_casename :
   %        root resolution when output_root is unset.
   %    dt_out : model-met output timestep (default "15m") for native and
   %        gridded forcing; pass "" to retain native model-met cadence.
   %        PROMICE userdata/Data remains at its native hourly cadence.
   %        Native met schema completion is fixed at the importer boundary:
   %        absent required channels are retained as NaN placeholders.
   %        Call buildPromiceMet directly for strict source-schema validation.
   %    overwrite : logical (default false). Existing requested artifacts are
   %        additive no-ops while missing artifacts/sources may be added. True
   %        explicitly refreshes a requested site's own artifacts. Other sites
   %        are never touched (see merge below).
   %    overwrite_family : logical (default false). Force a FULL rewrite of the
   %        family manifest from the requested sites alone, discarding other
   %        committed cases. The DEFAULT is MERGE.
   %    skip_missing : logical (default true). Record skip reasons and continue.
   %    dry_run : logical (default false). Return the manifest shape without
   %        writing eval/input artifacts or merging the family manifest.
   %    build_forcing : logical (default false). When true, runtime artifacts for
   %        forcing_sources are written. When false, no input/met or
   %        input/userdata artifacts are written; only observations are staged.
   %    build_observations : logical (default true). When false, reuse an
   %        existing manifest case and observations.mat contract, then update
   %        only requested forcing sources. This is a guarded fast path for
   %        already-staged cases, not a way to create new cases.
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
   %  See also: icemodel.verification.setup.importEsmSnowmip,
   %    icemodel.verification.setup.promiceSiteCatalog,
   %    icemodel.forcing.buildPromiceMet, icemodel.forcing.buildPromiceData,
   %    icemodel.forcing.buildMarMet, icemodel.forcing.buildMerraMet,
   %    icemodel.forcing.buildRacmoData

   arguments
      source_dir (1, 1) string = ""
      kwargs.case_ids (1, :) string = strings(1, 0)
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources, "promice")} = "promice"
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
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "15m"])} = "15m"
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = true
      kwargs.dry_run (1, 1) logical = false
      kwargs.build_forcing (1, 1) logical = false
      kwargs.build_observations (1, 1) logical = true
   end

   % Preserve the established PROMICE-specific catch boundary for empty builds.
   raw_forcing_sources = reshape(string(kwargs.forcing_sources), 1, []);
   if kwargs.build_forcing ...
         && all(strlength(strtrim(raw_forcing_sources)) == 0)
      error('icemodel:verification:importPromiceSites:emptyForcingSources', ...
         ['forcing_sources cannot be empty when build_forcing=true. ' ...
         'Use build_forcing=false for observation-only staging.'])
   end

   % Validate the optional clamp before any cache or staging side effect.
   [window_start, window_end, window_enabled] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);

   % Resolve the forcing sources.
   forcing_sources = ...
      icemodel.verification.setup.normalizeForcingSources( ...
      kwargs.forcing_sources, kwargs.build_forcing);
   kwargs.forcing_sources = forcing_sources;

   % Resolve the family identity and requested runtime source sets once.
   dataset_family = "promice";
   build_native_forcing = kwargs.build_forcing ...
      && ismember(dataset_family, forcing_sources);
   rcm_sources = intersect(forcing_sources, ...
      icemodel.verification.namelists.rcmsources(), "stable");
   build_rcm_forcing = kwargs.build_forcing && ~isempty(rcm_sources);

   % Resolve output roots and paths before raw sources. Forcing-only calls can
   % reuse the existing manifest without requiring observation/native caches.
   [evaluation_data_root, input_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   [family_root, manifest_file, met_outdir, userdata_outdir] = ...
      icemodel.verification.setup.datasetFamilyStagingPaths( ...
      evaluation_data_root, input_root, dataset_family);

   % Validate caches only when building observations or native runtime files.
   % Dry runs remain metadata-only; optional skips stay quiet while required
   % PROMICE products print their retrieval guidance before failing.
   needs_native_source = kwargs.build_observations || build_native_forcing;
   kwargs.promice_dir = resolvePromiceDir(source_dir, kwargs.promice_dir);
   if kwargs.dry_run
      kwargs.promice_dir = icemodel.forcing.helpers.verificationSourceDir( ...
         kwargs.promice_dir, dataset_family);
   elseif needs_native_source
      [kwargs.promice_dir, ~] = icemodel.verification.setup.fetchPromice( ...
         cache_dir=kwargs.promice_dir, stations=kwargs.case_ids, products="hour", ...
         strict=~kwargs.skip_missing, silent=kwargs.skip_missing);
   end

   % Resolve the case IDs.
   % Normal imports discover the full source catalog. A forcing-only call instead
   % derives omitted sites from the durable manifest without reading PROMICE data.
   if ~isempty(kwargs.case_ids)
      cases = reshape(kwargs.case_ids, 1, []);
   elseif ~kwargs.dry_run && ~kwargs.build_observations
      cases = priorPromiceSites(manifest_file);
   else
      cases = defaultPromiceSites(kwargs.promice_dir, kwargs.dry_run);
   end
   % PROMICE station IDs retain underscores; canonical manifest case IDs do not.
   requested_ids = lower(erase(cases, "_"));

   % AWS metadata is needed only while rebuilding observation/site contracts.
   % Fast paths preserve the existing manifest location and classification.
   aws_sites = table();
   if ~kwargs.dry_run && kwargs.build_observations
      aws_sites = readAwsSitesMetadata(kwargs.promice_dir);
   end

   coverage = struct();
   reuse_sources = strings(1, 0);

   if ~kwargs.dry_run && kwargs.build_observations
      icemodel.helpers.ensureDirExists(family_root);
   end
   if ~kwargs.dry_run && kwargs.build_forcing
      icemodel.helpers.ensureDirExists(met_outdir);
      icemodel.helpers.ensureDirExists(userdata_outdir);
   end

   % Resolve RCM coverage only for a real requested build.
   if ~kwargs.dry_run && build_rcm_forcing
      reuse_sources = rcm_sources;
      coverage = icemodel.verification.setup.promiceSourceCoverage( ...
         rcm_sources, struct('mar', kwargs.mar_dir, ...
         'merra', kwargs.merra_dir, 'racmo', kwargs.racmo_dir));
   end

   if ~kwargs.dry_run && ~kwargs.build_observations ...
         && ~build_native_forcing
      % Reuse the staged case entry so an RCM-only attachment does not require
      % observation or native-source caches.
      [state, alive, skipped] = ...
         icemodel.verification.setup.reuseDatasetFamilyCases( ...
         manifest_file, requested_ids, emptyState(), ...
         dataset_family=dataset_family, ...
         overwrite_family=kwargs.overwrite_family, ...
         forcing_sources=reuse_sources, coverage=coverage, ...
         startdate=kwargs.startdate, enddate=kwargs.enddate);
   else
      prior_cases = loadPriorPromiceCases(manifest_file, cases, kwargs);
      proj = icemodel.forcing.helpers.psnProjection();
      stage_callback = @(~, n) stagePromiceSite( ...
         cases(n), family_root, met_outdir, userdata_outdir, proj, aws_sites, ...
         window_enabled, window_start, window_end, coverage, rcm_sources, ...
         prior_cases, build_native_forcing, dataset_family, kwargs);

      % Stage each requested case into importer state.
      [state, alive, skipped] = ...
         icemodel.verification.setup.stageDatasetFamilyCases( ...
         1:numel(requested_ids), emptyState(), stage_callback, ...
         skip_missing=kwargs.skip_missing, ...
         warning_id="icemodel:verification:importPromiceSites:caseSkipped", ...
         label_callback=@(~, n) requested_ids(n));
   end

   % Record source provenance once so dry-run and persisted manifests match.
   source_doi = "";
   source_url = "https://promice.org";
   source_version = "pypromice-L3-hour";
   retrieval_date = string(datetime('today'));

   entry_callback = @icemodel.verification.setup.stateCaseEntry;

   % Return metadata without writing case or manifest artifacts.
   if kwargs.dry_run
      manifest = icemodel.verification.setup.runDatasetFamilyDryRun( ...
         state, alive, dataset_family=dataset_family, ...
         requested_ids=requested_ids, skipped=skipped, ...
         source_doi=source_doi, source_url=source_url, ...
         source_version=source_version, retrieval_date=retrieval_date, ...
         entry_callback=entry_callback);
      return
   end

   % Persist case entries first, then attach requested RCM source legs.
   [manifest, ~] = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family=dataset_family, ...
      manifest_file=manifest_file, requested_ids=requested_ids, ...
      skipped=skipped, ...
      source_doi=source_doi, source_url=source_url, ...
      source_version=source_version, retrieval_date=retrieval_date, ...
      overwrite_family=kwargs.overwrite_family, overwrite=kwargs.overwrite, ...
      entry_callback=entry_callback, ...
      build_forcing=build_rcm_forcing, forcing_sources=rcm_sources, ...
      leg_callback=@(s, src) s.leg.(char(src)), ...
      met_outdir=met_outdir, userdata_outdir=userdata_outdir, ...
      mar_dir=kwargs.mar_dir, merra_dir=kwargs.merra_dir, ...
      racmo_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
      method="nearest", dt_out=kwargs.dt_out);
end

%% Local helpers
function s = emptyState()
   %EMPTYSTATE Prototype dataset-family staging state.
   s = struct('site_id', "", 'case_id', "", ...
      'site_name', "", 'site_location', struct(), 'point', [NaN NaN], ...
      'period', struct('start', '', 'end', ''), 'evaluation_file_rel', '', ...
      'entry', struct(), 'colocation', struct(), 'leg', struct(), ...
      'comparison_variables', {strings(0, 1)}, ...
      'observation_variables', struct(), ...
      'surface_zone', "", 'eval_target', {strings(0, 1)}, ...
      'permafrost_zone', "", 'notes', "", 'reuse_entry', false, ...
      'dry_run', false);
end

function prior_cases = loadPriorPromiceCases(manifest_file, cases, kwargs)
   %LOADPRIORPROMICECASES Read existing contracts for merge and fast updates.
   prior_cases = struct([]);
   if kwargs.dry_run
      return
   end
   if ~isfile(manifest_file)
      % Observation-building calls can create new cases; the fast path cannot.
      if ~kwargs.build_observations
         error('icemodel:verification:importPromiceSites:missingPriorManifest', ...
            ['build_observations=false requires an existing PROMICE manifest: ' ...
            '%s'], manifest_file);
      end
      return
   end

   prior_entries = ...
      icemodel.verification.setup.loadPriorDatasetFamilyCases( ...
      manifest_file, overwrite_family=kwargs.overwrite_family, ...
      build_observations=kwargs.build_observations);
   requested = lower(erase(reshape(cases, 1, []), "_"));
   if isempty(prior_entries) || ~isfield(prior_entries, 'case_id')
      ids = strings(1, 0);
   else
      ids = string({prior_entries.case_id});
   end
   prior_cells = cell(1, numel(requested));
   n_prior = 0;
   for n = 1:numel(requested)
      hit = find(ids == requested(n), 1);
      if isempty(hit)
         % A normal observation build may add a new site during a merge update.
         if ~kwargs.build_observations
            error('icemodel:verification:importPromiceSites:missingPriorCase', ...
               ['build_observations=false requires existing case "%s" in ' ...
               '%s'], requested(n), manifest_file);
         end
         continue
      end
      if ~kwargs.build_observations
         assertPriorWindowCompatible(prior_entries(hit), kwargs);
      end
      n_prior = n_prior + 1;
      prior_cells{n_prior} = prior_entries(hit);
   end
   if n_prior > 0
      prior_cases = [prior_cells{1:n_prior}];
   end
end

function assertPriorWindowCompatible(prior_case, kwargs)
   %ASSERTPRIORWINDOWCOMPATIBLE Guard fast refresh against period drift.
   [window_start, window_end, window_enabled] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);
   if ~window_enabled
      return
   end

   [ok, prior_start, prior_end] = priorCasePeriod(prior_case);
   if ~ok
      error('icemodel:verification:importPromiceSites:unboundedPriorPeriod', ...
         ['build_observations=false requires a bounded existing period for ' ...
         'case "%s".'], string(prior_case.case_id));
   end

   if window_start < prior_start || window_end > prior_end
      error('icemodel:verification:importPromiceSites:priorWindowConflict', ...
         ['Requested window %s to %s is outside existing case "%s" period ' ...
         '%s to %s. Rebuild observations for a new period.'], ...
         string(kwargs.startdate), string(kwargs.enddate), ...
         string(prior_case.case_id), string(prior_case.period.start), ...
         string(prior_case.period.end));
   end
end

function s = stagePromiceSite(site, family_root, met_outdir, userdata_outdir, ...
      proj, aws_sites, window_enabled, window_start, window_end, coverage, ...
      rcm_sources, prior_cases, build_native_forcing, dataset_family, kwargs)
   %STAGEPROMICESITE Stage one PROMICE site and return importer state.
   case_id = lower(erase(site, "_"));
   prior_case = icemodel.verification.setup.priorCaseById( ...
      prior_cases, case_id);

   % A forcing-only refresh reuses the staged site/period contract and never
   % reads PROMICE metadata unless the native PROMICE leg itself was requested.
   if ~kwargs.build_observations && ~isempty(fieldnames(prior_case))
      [~, promice_start, promice_end] = priorCasePeriod(prior_case);
      site_location = prior_case.site_location;
      lat = site_location.lat_wgs84;
      lon = site_location.lon_wgs84;
      anchor = priorCaseAnchor(prior_case);
   else
      % Normal imports resolve source metadata before output writes so missing
      % files or empty windows remain whole-site skips.
      [aws_meta, promice_start, promice_end] = resolvePromiceSiteMetadata( ...
         site, window_enabled, window_start, window_end, kwargs);
      lat = aws_meta.lat;
      lon = aws_meta.lon;
      if all(isfinite([lat, lon]))
         [x3413, y3413] = projfwd(proj, lat, lon);
      else
         x3413 = NaN;
         y3413 = NaN;
      end
      site_location = struct( ...
         'lat_wgs84', lat, 'lon_wgs84', lon, ...
         'x_epsg3413', x3413, 'y_epsg3413', y3413, ...
         'elev_m', aws_meta.elev);
      anchor = siteCatalogEntry(site, aws_sites);
   end
   point = [lat, lon];

   leg = struct();
   if ~kwargs.dry_run && kwargs.build_forcing && ~isempty(rcm_sources)
      leg = icemodel.verification.setup.resolveLegWindows( ...
         rcm_sources, coverage, promice_start, promice_end);
      fprintf('[coverage] %s (PROMICE %d-%d):\n', site, ...
         year(promice_start), year(promice_end));
      icemodel.verification.setup.reportPromiceCoverage(coverage, ...
         [year(promice_start), year(promice_end)], ...
         legReportWindows(leg, rcm_sources));
   end

   case_root = fullfile(family_root, case_id);

   [colocation, comparison_variables, observation_variables, ...
      evaluation_file_rel] = ...
      stageNativePromice(site, case_id, case_root, userdata_outdir, ...
      met_outdir, promice_start, promice_end, prior_case, ...
      build_native_forcing, dataset_family, kwargs);

   s = struct('site_id', site, 'case_id', case_id, ...
      'site_name', anchor.site_name, 'site_location', site_location, ...
      'point', point, 'period', ...
      icemodel.verification.setup.manifestWindow( ...
      promice_start, promice_end), 'evaluation_file_rel', evaluation_file_rel, ...
      'entry', struct(), 'colocation', colocation, 'leg', leg, ...
      'comparison_variables', {comparison_variables}, ...
      'observation_variables', observation_variables, ...
      'surface_zone', anchor.surface_zone, ...
      'eval_target', {string(anchor.eval_target)}, ...
      'permafrost_zone', anchor.permafrost_zone, 'notes', anchor.note, ...
      'reuse_entry', false, ...
      'dry_run', kwargs.dry_run);
   s.entry = promiceCaseEntry(s);
end

function anchor = priorCaseAnchor(prior_case)
   %PRIORCASEANCHOR Preserve staged classification without reading raw sources.
   anchor = struct('site_name', string(prior_case.site_name), ...
      'surface_zone', string(prior_case.surface_zone), ...
      'eval_target', string(prior_case.eval_target), ...
      'permafrost_zone', string(prior_case.permafrost_zone), ...
      'note', string(prior_case.notes));
end

function [colocation, comparison_variables, observation_variables, ...
      evaluation_file_rel] = ...
      stageNativePromice(site, case_id, case_root, userdata_outdir, ...
      met_outdir, promice_start, promice_end, prior_case, ...
      build_native_forcing, dataset_family, kwargs)
   %STAGENATIVEPROMICE Stage native PROMICE observations and requested forcing.
   promice_data_files = strings(1, 0);
   % Dry runs and observation-reuse calls never enter the observation writer.
   write_observation = false;
   if ~kwargs.dry_run && kwargs.build_observations
      [promice_data, ~] = icemodel.forcing.buildPromiceData( ...
         site, source_dir=kwargs.promice_dir, ...
         startdate=promice_start, enddate=promice_end);
      observation_metadata = struct('source', 'promice_obs', ...
         'source_family', char(dataset_family), 'station', char(site), ...
         'site_id', char(site));
      requested_case = struct('period', ...
         icemodel.verification.setup.manifestWindow( ...
         promice_start, promice_end), ...
         'artifact_metadata', observation_metadata);
      write_observation = icemodel.verification.setup.prepareCaseRoot( ...
         case_root, kwargs.overwrite, "observations.mat", requested_case);
      if build_native_forcing
         promice_data_files = icemodel.forcing.helpers.writeuserdata( ...
            promice_data, case_id, dataset_family, outdir=userdata_outdir, ...
            naming="window", overwrite=kwargs.overwrite);
      end
      [comparison_variables, observation_variables] = ...
         firnComparisonContract(promice_data);
   elseif ~kwargs.dry_run && ~kwargs.build_observations
      [comparison_variables, observation_variables, evaluation_file_rel, ...
         promice_data_files] = ...
         priorObservationContract(prior_case);
      if build_native_forcing
         % Native forcing refreshes rebuild input Data without rewriting the
         % existing observations.mat evaluation contract.
         [promice_data, ~] = icemodel.forcing.buildPromiceData( ...
            site, source_dir=kwargs.promice_dir, ...
            startdate=promice_start, enddate=promice_end);
         promice_data_files = icemodel.forcing.helpers.writeuserdata( ...
            promice_data, case_id, dataset_family, outdir=userdata_outdir, ...
            naming="window", overwrite=kwargs.overwrite);
      end
   else
      [comparison_variables, observation_variables] = ...
         dryRunFirnComparisonContract();
   end

   if ~kwargs.dry_run && kwargs.build_observations && write_observation
      % Preserve an already-current observation bundle during ordinary
      % additive imports; explicit overwrite remains the refresh boundary.
      targets = struct('format', 'timeseries', ...
         'data', promice_data, 'metadata', observation_metadata);
      targets = icemodel.verification.setup.stampArtifactMetadata(targets);
      save(fullfile(case_root, 'observations.mat'), 'targets');
   end
   if kwargs.dry_run || kwargs.build_observations
      evaluation_file_rel = char(fullfile(case_id, 'observations.mat'));
   end

   preserve_prior_promice = ~kwargs.dry_run && ~build_native_forcing ...
      && ~kwargs.overwrite_family && hasPriorPromiceLeg(prior_case);
   if preserve_prior_promice
      % A merge refresh should add/update requested sources without erasing a
      % prior PROMICE leg merely because PROMICE was omitted from this call.
      promice_co = prior_case.colocation.promice;
   else
      % Evaluation availability is independent of native runtime staging. An
      % observation-only or RCM-only import owns observations.mat without
      % claiming nonexistent PROMICE met/Data artifacts.
      promice_co = struct('kind', 'station_met_and_eval', ...
         'staged', false, ...
         'eval_staged', ~kwargs.dry_run ...
         && isfile(fullfile(case_root, 'observations.mat')));
      promice_co.data_files = ...
         icemodel.verification.setup.relpaths(promice_data_files, userdata_outdir);
      promice_co.window = icemodel.verification.setup.manifestWindow( ...
         promice_start, promice_end);
      if kwargs.dry_run || ~build_native_forcing
         promice_co.met_files = strings(1, 0);
      else
         try
            promice_met = icemodel.forcing.buildPromiceMet(site, ...
               source_dir=kwargs.promice_dir, ...
               startdate=promice_start, enddate=promice_end, ...
               fillwithmissing=true);
            promice_met_files = icemodel.forcing.helpers.writemet( ...
               promice_met, case_id, dataset_family, outdir=met_outdir, ...
               naming="window", dt_out=kwargs.dt_out, ...
               overwrite=kwargs.overwrite);
            % Diagnose the exact returned path because no-overwrite staging may
            % select an existing exact or broader enclosing artifact.
            [forcing_ready, forcing_ready_reason, forcing_complete_windows] = ...
               icemodel.verification.setup.metArtifactReadiness(promice_met_files);
            promice_co.met_files = ...
               icemodel.verification.setup.relpaths(promice_met_files, met_outdir);
            promice_co.forcing_ready = logical(forcing_ready);
            promice_co.forcing_ready_reason = char(forcing_ready_reason);
            promice_co.forcing_complete_windows = forcing_complete_windows;
         catch met_err
            if ~isSkippablePromiceMetError(met_err)
               rethrow(met_err)
            end
            promice_co.met_files = strings(1, 0);
            promice_co.met_skipped_reason = string(met_err.message);
         end
      end

      % A native leg is staged only when at least one real runtime artifact is
      % selectable. This deliberately keeps a one-sample Data-only leg staged
      % when PROMICE met construction is inapplicable.
      promice_co.staged = ~isempty(promice_co.data_files) ...
         || ~isempty(promice_co.met_files);
   end

   colocation = struct('promice', promice_co);
end

function entry = promiceCaseEntry(s)
   %PROMICECASEENTRY Convert one PROMICE state record to a manifest case entry.
   [forcing_sources, eval_sources] = ...
      icemodel.verification.setup.colocationSourceLists(s.colocation);
   if s.dry_run && isempty(eval_sources)
      eval_sources = "promice_obs";
   end

   case_values = { ...
      char(s.case_id)
      'firn_observational'
      char(s.site_id)
      char(s.site_name)
      char(s.surface_zone)
      cellstr(s.eval_target)
      char(s.permafrost_zone)
      s.site_location
      s.period
      s.evaluation_file_rel
      cellstr(forcing_sources(:))
      cellstr(eval_sources(:))
      cellstr(s.comparison_variables(:))
      s.observation_variables
      s.colocation
      'hourly'
      char(s.notes)};

   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);
end

function [comparison_variables, observation_variables] = ...
      dryRunFirnComparisonContract()
   %DRYRUNFIRNCOMPARISONCONTRACT PROMICE manifest axes used without raw files.
   comparison_variables = ["ablation"; "snow_depth"; "tsfc"; "tice10m"; ...
      "tice1"; "tice2"; "tice3"; "tice4"; ...
      "tice5"; "tice6"; "tice7"; "tice8"];
   observation_variables = icemodel.verification.setup.metadataStruct({ ...
      'subsurface_temperature_primary', 'tice10m [K] (standardized 10 m below surface)'
      'subsurface_temperature_string', 'tice1..tice8 (raw string thermistors)'
      'surface_ablation', 'ablation [m, positive down]'
      'snow_depth', 'snow_depth [m]'});
end

function [comparison_variables, observation_variables, evaluation_file_rel, ...
      data_files] = ...
      priorObservationContract(prior_case)
   %PRIOROBSERVATIONCONTRACT Reuse an existing observation manifest contract.
   comparison_variables = string(prior_case.comparison_variables);
   comparison_variables = reshape(comparison_variables, [], 1);
   observation_variables = prior_case.observation_variables;
   evaluation_file_rel = char(prior_case.evaluation_file);
   data_files = strings(1, 0);
   if isfield(prior_case, 'colocation') ...
         && isfield(prior_case.colocation, 'promice') ...
         && isfield(prior_case.colocation.promice, 'data_files')
      data_files = reshape(string( ...
         prior_case.colocation.promice.data_files), 1, []);
   end
end

function tf = hasPriorPromiceLeg(prior_case)
   %HASPRIORPROMICELEG True when an existing manifest has a PROMICE leg.
   tf = isfield(prior_case, 'colocation') ...
      && isfield(prior_case.colocation, 'promice') ...
      && isstruct(prior_case.colocation.promice);
end

function [ok, t1, t2] = priorCasePeriod(prior_case)
   %PRIORCASEPERIOD Parse a bounded existing manifest period.
   ok = false;
   t1 = NaT('TimeZone', 'UTC');
   t2 = NaT('TimeZone', 'UTC');
   if ~isfield(prior_case, 'period') ...
         || ~isfield(prior_case.period, 'start') ...
         || ~isfield(prior_case.period, 'end') ...
         || strlength(string(prior_case.period.start)) == 0 ...
         || strlength(string(prior_case.period.end)) == 0
      return
   end
   [t1, t2] = ...
      icemodel.verification.setup.periodBounds(prior_case.period);
   ok = ~isnat(t1) && ~isnat(t2);
end

function tf = isSkippablePromiceMetError(err)
   %ISSKIPPABLEPROMICEMETERROR True for absent-source PROMICE met failures.
   ids = [
      "icemodel:forcing:readPromiceAws:sourceNotFound"
      "icemodel:forcing:readPromiceAws:stationNotFound"
      "icemodel:forcing:readPromiceAws:emptyWindow"
      "icemodel:forcing:validatemet:tooFewSamples"];
   tf = any(string(err.identifier) == ids);
end

function sites = defaultPromiceSites(promice_dir, dry_run)
   %DEFAULTPROMICESITES Resolve the default site list for import preflight.
   if dry_run
      catalog = icemodel.verification.setup.promiceSiteCatalog();
      sites = string({catalog.site});
   else
      sites = icemodel.verification.namelists.promicesite(promice_dir);
   end
end

function sites = priorPromiceSites(manifest_file)
   %PRIORPROMICESITES Return staged site ids without reading raw sources.
   if ~isfile(manifest_file)
      error('icemodel:verification:importPromiceSites:missingPriorManifest', ...
         ['build_observations=false requires an existing PROMICE manifest: ' ...
         '%s'], manifest_file)
   end
   manifest = jsondecode(fileread(manifest_file));
   if ~isfield(manifest, 'cases') || isempty(manifest.cases)
      error('icemodel:verification:importPromiceSites:missingPriorCase', ...
         'PROMICE manifest contains no cases: %s', manifest_file)
   end

   % site_id preserves the upstream station spelling; older fixtures may carry
   % only case_id, which is still sufficient for a previously staged case.
   cases = manifest.cases;
   if isfield(cases, 'site_id')
      sites = reshape(string({cases.site_id}), 1, []);
   else
      sites = reshape(string({cases.case_id}), 1, []);
   end
end

function [aws_meta, promice_start, promice_end] = ...
      resolvePromiceSiteMetadata(site, window_enabled, window_start, ...
      window_end, kwargs)
   %RESOLVEPROMICESITEMETADATA Resolve source metadata or dry-run placeholders.
   if kwargs.dry_run
      aws_meta = struct('lat', NaN, 'lon', NaN, 'elev', NaN);
      if window_enabled
         promice_start = window_start;
         promice_end = window_end;
      else
         promice_start = NaT;
         promice_end = NaT;
      end
      return
   end

   if window_enabled
      [~, aws_meta] = icemodel.forcing.readPromiceAws(site, ...
         source_dir=kwargs.promice_dir, timescale="hourly", ...
         startdate=window_start, enddate=window_end);
      promice_start = aws_meta.window_start;
      promice_end = aws_meta.window_end;
   else
      [~, aws_meta] = icemodel.forcing.readPromiceAws(site, ...
         source_dir=kwargs.promice_dir, timescale="hourly");
      promice_start = aws_meta.window_start;
      promice_end = aws_meta.window_end;
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
   % Curated and first-pass classifications both live in promiceSiteCatalog (the
   % single source of truth for surface_zone + eval_target). When a station is
   % not cataloged there (e.g. a brand-new L3 station absent from the AWS CSV),
   % fall back to "unknown" so the manifest still validates.
   try
      info = icemodel.verification.setup.promiceSiteCatalog(site);
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

function w = legReportWindows(leg, sources)
   %LEGREPORTWINDOWS Flatten the per-leg windows for the coverage report.
   w = struct();
   for m = icemodel.verification.namelists.rcmsources()
      if ~ismember(m, sources)
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

function [comparison_variables, observation_variables] = ...
      firnComparisonContract(promice_data)
   %FIRNCOMPARISONCONTRACT Comparison axes + obs metadata from PROMICE Data.
   %
   % Records which firn comparison variables are present in the staged PROMICE
   % Data record and the observation metadata. This is the comparison CONTRACT
   % only - the observation data itself lives in the bundled observations.mat
   % eval target and the per-year userdata files, not in the manifest.
   if isempty(promice_data)
      comparison_variables = strings(0, 1);
      observation_variables = icemodel.verification.setup.metadataStruct({ ...
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
   comparison_variables = canonical(ismember(canonical, present));

   observation_variables = icemodel.verification.setup.metadataStruct({ ...
      'subsurface_temperature_primary', 'tice10m [K] (standardized 10 m below surface)'
      'subsurface_temperature_string', 'tice1..tice8 (raw string thermistors)'
      'surface_ablation', 'ablation [m, positive down]'
      'snow_depth', 'snow_depth [m]'});
end
