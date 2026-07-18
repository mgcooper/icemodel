function manifest = importImau(source_dir, kwargs)
   %IMPORTIMAU Stage the IMAU hourly AWS verification family.
   %
   %  manifest = icemodel.verification.setup.importImau(source_dir)
   %  manifest = icemodel.verification.setup.importImau(source_dir, ...
   %     case_ids="S21", dry_run=true)
   %
   %  Role
   %    First-pass IMAU staging hook for the hourly PANGAEA S21/S22/S23 network.
   %    The daily 19-station product remains QA/provenance input, not a first-pass
   %    case inventory.
   %    forcing_sources selects runtime sources requested by the current call.
   %    Ordinary calls preserve omitted existing legs; overwrite_family=true
   %    deliberately replaces the whole family state.
   %    build_observations=false is a guarded non-dry fast path: requested cases
   %    must already exist in the target manifest, whose observation entry is
   %    reused while selected forcing is attached.
   %
   %  Default roots
   %    source_dir="" reads <repo>/data/verification/imau. With no output_root,
   %    observations go to <repo>/data/eval/imau/<case_id>/observations.mat and
   %    native met/userdata go to <repo>/data/input/{met,userdata}/imau/.
   %    Explicit source_dir, output_root, evaluation_data_root, and
   %    input_data_root overrides are honored as-is.
   %
   %  Met and userdata
   %    Model met defaults to dt_out="15m"; pass dt_out="" for native cadence.
   %    Data/userdata defaults to hourly at the shared writer boundary.
   %    Native met schema completion is fixed at the importer boundary: absent
   %    required channels are retained as NaN placeholders.
   %    Call buildImauHourlyMet directly for strict source-schema validation.
   %
   %  Window selection
   %    startdate/enddate optionally clamp each requested IMAU hourly record.
   %    This is intended for short preview staging; omit both for full
   %    production artifacts.
   %
   %  Name-value
   %    case_ids : string vector  IMAU site cases; blank selects the catalog.
   %    forcing_sources : string vector  Native/RCM legs requested by this call.
   %    startdate, enddate : paired optional observation/forcing clamp.
   %    output_root : string  Convenience root for eval and input artifacts.
   %    overwrite : logical  Refresh requested case artifacts (default false).
   %    overwrite_family : logical  Replace the whole family (default false).
   %    skip_missing : logical  Record source-local skips instead of failing.
   %    dry_run : logical  Return metadata without source or output writes.
   %    build_forcing : logical  Stage requested runtime legs (default false).
   %    build_observations : logical  Build or reuse the case observation entry.
   %
   %  Returns
   %    manifest : struct  Final or dry-run family manifest.
   %
   %  See also: icemodel.verification.setup.importPromiceSites,
   %    icemodel.verification.setup.importRetmip

   arguments
      source_dir (1, 1) string = ""
      kwargs.case_ids (1, :) string = strings(1, 0)
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources, "imau")} = "imau"
      kwargs.startdate = ""
      kwargs.enddate = ""
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

   % Resolve the case IDs.
   cases = icemodel.verification.setup.imauSiteCatalog(kwargs.case_ids);
   requested_ids = string({cases.case_id});

   % Validate the optional clamp before any cache or staging side effect.
   % IMAU builders consume the original name-values directly, so only paired
   % input validation is needed here.
   icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);

   % Resolve the forcing sources.
   forcing_sources = ...
      icemodel.verification.setup.normalizeForcingSources( ...
      kwargs.forcing_sources, kwargs.build_forcing);
   kwargs.forcing_sources = forcing_sources;

   % Resolve the family identity and requested runtime source sets once.
   dataset_family = "imau";
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

   prior_cases = struct([]);
   coverage = struct();
   reuse_sources = strings(1, 0);

   if ~kwargs.dry_run
      prior_cases = ...
         icemodel.verification.setup.loadPriorDatasetFamilyCases( ...
         manifest_file, overwrite_family=kwargs.overwrite_family, ...
         build_observations=kwargs.build_observations);
   end

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
      % Validate caches only when building observations or native runtime files.
      % Dry runs remain metadata-only; optional skips stay quiet while required
      % IMAU products print their retrieval guidance before failing.
      cache_status = struct();
      cache_products = "hourly";
      validate_daily_qa = ~kwargs.dry_run && kwargs.build_observations;
      if validate_daily_qa
         cache_products = ["hourly", "daily"];
      end
      if ~kwargs.dry_run
         strict_cache = ~kwargs.skip_missing;
         [source_dir, cache_status] = icemodel.verification.setup.fetchImau( ...
            cache_dir=icemodel.forcing.helpers.verificationSourceDir( ...
            source_dir, dataset_family), products=cache_products, ...
            strict=strict_cache, ...
            silent=kwargs.skip_missing, create_cache_dir=true);
      end
      if ~kwargs.dry_run && ~kwargs.skip_missing
         preflightCases(cases, source_dir);
      end
      daily_qa = dailyQaStatus(source_dir, [cases.site_id], ...
         validate=validate_daily_qa, ...
         strict=(validate_daily_qa && ~kwargs.dry_run && ~kwargs.skip_missing));
      stage_callback = @(~, n) stageImauCase( ...
         cases(n), source_dir, cache_status, family_root, met_outdir, ...
         userdata_outdir, kwargs, daily_qa, coverage, rcm_sources, ...
         prior_cases, build_native_forcing, dataset_family);

      % Stage each requested case into importer state.
      [state, alive, skipped] = ...
         icemodel.verification.setup.stageDatasetFamilyCases( ...
         1:numel(requested_ids), emptyState(), stage_callback, ...
         skip_missing=kwargs.skip_missing, ...
         warning_id="icemodel:verification:importImau:caseSkipped", ...
         label_callback=@(~, n) requested_ids(n));
   end

   % Record source provenance once so dry-run and persisted manifests match.
   source_doi = "10.1594/PANGAEA.971647;10.1594/PANGAEA.970127";
   source_url = "https://doi.org/10.1594/PANGAEA.971647";
   source_version = "hourly-s21-s22-s23+daily-qa";
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

function preflightCases(cases, source_dir)
   %PREFLIGHTCASES Validate every strict IMAU case before staging begins.
   for k = 1:numel(cases)
      site_id = cases(k).site_id;
      try
         icemodel.forcing.helpers.locateImauHourlyFile(source_dir, site_id);
      catch err
         if ~strcmp(err.identifier, ...
               'icemodel:forcing:locateImauHourlyFile:fileNotFound')
            rethrow(err)
         end
         error('icemodel:verification:importImau:missingHourlySource', ...
            '%s', err.message)
      end
   end
end

function state = stageImauCase(site, source_dir, cache_status, family_root, ...
      met_outdir, userdata_outdir, kwargs, daily_qa, coverage, rcm_sources, ...
      prior_cases, build_native_forcing, dataset_family)
   %STAGEIMAUCASE Stage one IMAU hourly case and return manifest state.
   if kwargs.dry_run
      state = dryRunImauCase(site, source_dir, cache_status, daily_qa);
      return
   end

   case_id = string(site.case_id);
   prior_case = icemodel.verification.setup.priorCaseById( ...
      prior_cases, case_id);
   if ~kwargs.build_observations && isempty(prior_case)
      error('icemodel:verification:importImau:missingPriorCase', ...
         ['build_observations=false requires existing case "%s" in the ' ...
         'target IMAU manifest.'], case_id);
   end
   if ~kwargs.build_observations && isfield(prior_case, 'colocation') ...
         && isfield(prior_case.colocation, 'daily_qa')
      % Native-only refreshes preserve prior daily QA without reopening that cache.
      daily_qa = prior_case.colocation.daily_qa;
   end

   % A native-only refresh defaults to the staged observation window instead of
   % silently widening the runtime artifact to the full source record.
   build_start = kwargs.startdate;
   build_end = kwargs.enddate;
   if ~kwargs.build_observations && strlength(string(build_start)) == 0
      build_start = string(prior_case.period.start);
      build_end = string(prior_case.period.end);
   end
   [met, metadata, Data] = icemodel.forcing.buildImauHourlyMet( ...
      site.site_id, source_dir=source_dir, startdate=build_start, ...
      enddate=build_end, fillwithmissing=true);

   % Resolve the saved/manifest contract before deciding whether an existing
   % fixed-name observation artifact still covers this request.
   site_location = icemodel.forcing.helpers.projectLocation( ...
      metadata.site_location);
   metadata.site_location = site_location;
   period = struct('start', ...
      icemodel.verification.setup.formatManifestTime(Data.Time(1)), ...
      'end', icemodel.verification.setup.formatManifestTime(Data.Time(end)));
   requested_case = struct('period', period, ...
      'site_location', site_location, 'artifact_metadata', metadata);

   write_observation = false;
   if kwargs.build_observations
      case_root = fullfile(family_root, case_id);
      write_observation = icemodel.verification.setup.prepareCaseRoot( ...
         case_root, kwargs.overwrite, "observations.mat", requested_case);
   end

   met_files = strings(1, 0);
   data_files = strings(1, 0);
   forcing_ready = false;
   if kwargs.build_forcing
      forcing_ready_reason = ...
         dataset_family + " was not requested in forcing_sources";
   else
      forcing_ready_reason = "native forcing disabled because build_forcing=false";
   end
   % No-forcing calls publish no saved-artifact window diagnostics; a compatible
   % prior forcing leg can selectively restore them below.
   forcing_complete_windows = repmat(struct( ...
      'start_time', "", 'end_time', "", 'sample_count', 0), 0, 1);
   if build_native_forcing
      % Native input artifacts obey the suite-wide forcing toggle; the same
      % source Data can still back observations.mat when forcing is disabled.
      met_files = icemodel.forcing.helpers.writemet( ...
         met, case_id, dataset_family, outdir=met_outdir, naming="window", ...
         dt_out=kwargs.dt_out, ...
         overwrite=kwargs.overwrite);
      % Diagnose the exact returned path because no-overwrite staging may
      % select an existing exact or broader enclosing artifact.
      [forcing_ready, forcing_ready_reason, forcing_complete_windows] = ...
         icemodel.verification.setup.metArtifactReadiness(met_files);
      data_files = icemodel.forcing.helpers.writeuserdata(Data, case_id, ...
         dataset_family, outdir=userdata_outdir, naming="window", ...
         overwrite=kwargs.overwrite);
   end

   if write_observation
      % A repeated non-overwrite import retains the current observation bytes.
      targets = struct('format', 'timeseries', 'data', Data, ...
         'metadata', metadata);
      targets = icemodel.verification.setup.stampArtifactMetadata(targets);
      save(fullfile(case_root, 'observations.mat'), 'targets');
   end
   evaluation_file_rel = char(fullfile(case_id, 'observations.mat'));
   if ~kwargs.build_observations
      evaluation_file_rel = char(prior_case.evaluation_file);
   end

   colocation = imauColocation(source_dir, metadata, met_files, data_files, ...
      met_outdir, userdata_outdir, evaluation_file_rel, period, ...
      site.source_association, daily_qa, forcing_ready, ...
      forcing_ready_reason, forcing_complete_windows);
   [colocation, identity_conflict] = preservePriorImauLeg( ...
      colocation, prior_case, metadata, build_native_forcing, kwargs);
   if identity_conflict
      % Keep the fresh evaluation provenance but make the missing native rebuild
      % explicit; the prior files stay on disk without remaining manifest links.
      colocation.imau.kind = 'hourly_aws_eval';
      colocation.imau.forcing_ready = false;
      colocation.imau.forcing_ready_reason = ...
         ['native source identity changed; rerun with build_forcing=true ' ...
         'and forcing_sources including imau'];
   end
   leg = struct();
   if kwargs.build_forcing && ~isempty(rcm_sources)
      leg = icemodel.verification.setup.resolveLegWindows( ...
         rcm_sources, coverage, Data.Time(1), Data.Time(end));
   end

   state = struct('site_id', site.site_id, 'case_id', site.case_id, ...
      'site_name', site.site_name, ...
      'site_location', site_location, ...
      'point', [site_location.lat_wgs84, site_location.lon_wgs84], ...
      'period', period, 'evaluation_file_rel', evaluation_file_rel, ...
      'entry', struct(), 'colocation', colocation, 'leg', leg, ...
      'comparison_variables', {imauComparisonVariables(Data)}, ...
      'observation_variables', ...
      imauObservationVariables(metadata, Data, daily_qa), ...
      'surface_zone', site.surface_zone, ...
      'eval_target', {string(site.eval_target)}, ...
      'permafrost_zone', site.permafrost_zone, ...
      'notes', imauCaseNote(site), ...
      'reuse_entry', false, 'dry_run', false);
   if kwargs.build_observations
      state.entry = imauCaseEntry(state);
   else
      state.entry = prior_case;
   end
end

function [colocation, identity_conflict] = preservePriorImauLeg( ...
      colocation, prior_case, metadata, build_native_forcing, kwargs)
   %PRESERVEPRIORIMAULEG Merge native artifacts into a fresh observation leg.
   identity_conflict = false;
   if build_native_forcing || kwargs.overwrite_family || isempty(prior_case)
      return
   end
   if ~isfield(prior_case, 'colocation') ...
         || ~isfield(prior_case.colocation, 'imau')
      return
   end
   prior = prior_case.colocation.imau;

   % Native artifacts remain reusable only across the same stable upstream
   % product. Cache paths and source filenames may change during an observation
   % refresh, so identity rests on the product kind and DOI fields.
   identity_fields = ["kind", "doi", "bundle_doi"];
   for fieldname = identity_fields
      name = char(fieldname);
      if isfield(prior, name) && isfield(colocation.imau, name) ...
            && string(prior.(name)) ~= string(colocation.imau.(name))
         identity_conflict = true;
         return
      end
   end

   % Compare canonical and native IMAU producer names through the shared scalar
   % identity rules. Downstream source_association describes RetMIP FA and is not
   % the producer identity for these met/Data artifacts.
   fresh_identity = struct('family', metadata.source_family, ...
      'source_id', metadata.station);
   prior_metadata = struct();
   if isfield(prior, 'artifact_metadata') && ~isempty(prior.artifact_metadata)
      if ~isstruct(prior.artifact_metadata) ...
            || ~isscalar(prior.artifact_metadata)
         identity_conflict = true;
         return
      end
      prior_metadata = prior.artifact_metadata;
      if ~icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
            prior_metadata, fresh_identity)
         identity_conflict = true;
         return
      end
      native_identity = struct();
      if isfield(prior_metadata, 'source_family')
         native_identity.family = prior_metadata.source_family;
      end
      if isfield(prior_metadata, 'station')
         native_identity.source_id = prior_metadata.station;
      end
      if ~icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
            native_identity, fresh_identity)
         identity_conflict = true;
         return
      end
   end
   if isfield(prior_case, 'observation_variables') ...
         && isstruct(prior_case.observation_variables)
      prior_obs = prior_case.observation_variables;
      obs_identity = struct();
      if isfield(prior_obs, 'source_family')
         obs_identity.family = prior_obs.source_family;
      end
      if isfield(prior_obs, 'hourly_site')
         obs_identity.source_id = prior_obs.hourly_site;
      end
      if ~icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
            obs_identity, fresh_identity)
         identity_conflict = true;
         return
      end
   end

   % Point semantics remain local: every known prior case/artifact point must
   % agree with the fresh producer, while missing or nonfinite legacy points are
   % compatible. Accept direct or nested artifact metadata location records.
   fresh_point = [metadata.site_location.lat_wgs84, ...
      metadata.site_location.lon_wgs84];
   prior_points = nan(3, 2);
   n_points = 0;
   if isfield(prior_case, 'site_location') ...
         && isstruct(prior_case.site_location) ...
         && all(isfield(prior_case.site_location, ...
         ["lat_wgs84", "lon_wgs84"]))
      n_points = n_points + 1;
      prior_points(n_points, :) = [ ...
         prior_case.site_location.lat_wgs84, ...
         prior_case.site_location.lon_wgs84];
   end
   if all(isfield(prior_metadata, ["lat_wgs84", "lon_wgs84"]))
      n_points = n_points + 1;
      prior_points(n_points, :) = [prior_metadata.lat_wgs84, ...
         prior_metadata.lon_wgs84];
   end
   if isfield(prior_metadata, 'site_location') ...
         && isstruct(prior_metadata.site_location) ...
         && all(isfield(prior_metadata.site_location, ...
         ["lat_wgs84", "lon_wgs84"]))
      n_points = n_points + 1;
      prior_points(n_points, :) = [ ...
         prior_metadata.site_location.lat_wgs84, ...
         prior_metadata.site_location.lon_wgs84];
   end
   for n = 1:n_points
      if all(isfinite(prior_points(n, :))) && all(isfinite(fresh_point)) ...
            && any(abs(prior_points(n, :) - fresh_point) > 1e-8)
         identity_conflict = true;
         return
      end
   end

   % Retain only forcing-owned runtime references, status, and artifact
   % provenance. Fresh source/window/evaluation fields intentionally win.
   native_fields = ["met_files", "data_files", "forcing_ready", ...
      "forcing_ready_reason", "forcing_complete_windows", ...
      "artifact_metadata"];
   for fieldname = native_fields
      name = char(fieldname);
      if isfield(prior, name)
         colocation.imau.(name) = prior.(name);
      end
   end
end

function state = dryRunImauCase(site, source_dir, cache_status, daily_qa)
   %DRYRUNIMAUCASE Build one non-writing IMAU state record.
   imau = struct('kind', 'hourly_aws', 'staged', false, ...
      'source_dir', char(source_dir), 'cache_status', cache_status);
   colocation = struct('imau', imau, ...
      'source_association', site.source_association, ...
      'daily_qa', daily_qa);
   state = struct('site_id', site.site_id, 'case_id', site.case_id, ...
      'site_name', site.site_name, ...
      'site_location', struct('lat_wgs84', site.lat_wgs84, ...
      'lon_wgs84', site.lon_wgs84, 'x_epsg3413', site.x_epsg3413, ...
      'y_epsg3413', site.y_epsg3413, 'elev_m', site.elev_m), ...
      'point', [site.lat_wgs84, site.lon_wgs84], ...
      'period', struct('start', '', 'end', ''), ...
      'evaluation_file_rel', '', 'entry', struct(), ...
      'colocation', colocation, 'leg', struct(), ...
      'comparison_variables', {["tair", "rh", "wspd", "wdir", "psfc", ...
      "swd", "swu", "lwd", "lwu", "albedo", "tsfc", "surface_height"]}, ...
      'observation_variables', struct('hourly_site', site.site_id, ...
      'daily_qa_bundle', daily_qa.bundle_doi), ...
      'surface_zone', site.surface_zone, ...
      'eval_target', {string(site.eval_target)}, ...
      'permafrost_zone', site.permafrost_zone, ...
      'notes', imauCaseNote(site), ...
      'reuse_entry', false, 'dry_run', true);
   state.entry = imauCaseEntry(state);
end

function entry = imauCaseEntry(s)
   %IMAUCASEENTRY Convert one IMAU state record to a manifest case entry.
   [forcing_sources, eval_sources] = ...
      icemodel.verification.setup.colocationSourceLists(s.colocation);
   if s.dry_run && isempty(eval_sources)
      eval_sources = "imau_obs";
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
      '1hr'
      char(s.notes)};

   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);
end

function note = imauCaseNote(site)
   %IMAUCASENOTE Return the staged IMAU case description.
   note = sprintf('IMAU hourly AWS site %s. %s', site.site_id, site.note);
end

function colocation = imauColocation(source_dir, metadata, met_files, ...
      data_files, met_outdir, userdata_outdir, evaluation_file_rel, period, ...
      source_association, daily_qa, forcing_ready, forcing_ready_reason, ...
      forcing_complete_windows)
   %IMAUCLOCATION Build staged IMAU colocation metadata.
   imau = struct('kind', 'hourly_aws_met_and_eval', 'staged', true, ...
      'source_dir', char(source_dir), ...
      'source_file', char(metadata.source_file), ...
      'doi', char(metadata.doi), ...
      'bundle_doi', char(metadata.bundle_doi), ...
      'evaluation_file', evaluation_file_rel, ...
      'met_files', ...
      icemodel.verification.setup.relpaths(met_files, met_outdir), ...
      'data_files', ...
      icemodel.verification.setup.relpaths(data_files, userdata_outdir), ...
      'forcing_ready', logical(forcing_ready), ...
      'forcing_ready_reason', char(forcing_ready_reason), ...
      'window', period, ...
      'precip_policy', char(metadata.precip_policy));
   imau.forcing_complete_windows = forcing_complete_windows;
   colocation = struct('imau', imau, ...
      'source_association', source_association, ...
      'daily_qa', daily_qa);
end

function vars = imauComparisonVariables(Data)
   %IMAUCOMPARISONVARIABLES Return staged IMAU comparison/eval axes.
   candidates = ["tair", "rh", "wspd", "wdir", "psfc", "swd", "swu", ...
      "lwd", "lwu", "swn", "lwn", "netr", "albedo", "tsfc", ...
      "boom_height", "surface_height"];
   present = string(Data.Properties.VariableNames);
   vars = candidates(ismember(candidates, present));
end

function obs = imauObservationVariables(metadata, Data, daily_qa)
   %IMAUOBSERVATIONVARIABLES Keep compact hourly/daily provenance metadata.
   obs = struct( ...
      'hourly_site', metadata.station, ...
      'source_family', metadata.source_family, ...
      'doi', metadata.doi, ...
      'bundle_doi', metadata.bundle_doi, ...
      'native_variables', {cellstr(string(Data.Properties.VariableNames(:)))}, ...
      'daily_qa_validated', daily_qa.validated, ...
      'daily_qa_bundle', daily_qa.bundle_doi, ...
      'precip_policy', metadata.precip_policy);
end

function qa = dailyQaStatus(source_dir, site_ids, kwargs)
   %DAILYQASTATUS Validate daily SEB QA files without promoting them to cases.
   arguments
      source_dir (1, 1) string
      site_ids (1, :) string
      kwargs.validate (1, 1) logical = true
      kwargs.strict (1, 1) logical = false
   end

   records = repmat(dailyQaRecord(), 1, numel(site_ids));
   for k = 1:numel(site_ids)
      if kwargs.validate
         records(k) = inspectDailyQaFile(source_dir, site_ids(k), kwargs.strict);
      else
         records(k).site_id = site_ids(k);
         records(k).reason = "daily QA cache not read for this call";
      end
   end

   present = [records.present];
   valid = [records.valid];
   qa = struct( ...
      'kind', "daily_19_site_seb_qa", ...
      'staged', false, ...
      'validated', all(valid), ...
      'bundle_doi', "10.1594/PANGAEA.970127", ...
      'landing_url', "https://doi.org/10.1594/PANGAEA.970127", ...
      'required_sites', {cellstr(site_ids)}, ...
      'present_sites', {cellstr(site_ids(present))}, ...
      'missing_sites', {cellstr(site_ids(~present))}, ...
      'invalid_sites', {cellstr(site_ids(present & ~valid))}, ...
      'records', records);
end

function record = dailyQaRecord()
   %DAILYQARECORD Prototype one daily QA validation record.
   record = struct('site_id', "", 'present', false, 'valid', false, ...
      'file', "", 'doi', "", 'landing_url', "", 'rows_checked', 0, ...
      'reason', "");
end

function record = inspectDailyQaFile(source_dir, site_id, strict)
   %INSPECTDAILYQAFILE Check header/sample content for one daily QA file.
   record = dailyQaRecord();
   record.site_id = site_id;
   filename = locateDailyFile(source_dir, site_id);
   if filename == ""
      record.reason = sprintf('missing daily QA file for %s', site_id);
      if strict
         error('icemodel:verification:importImau:missingDailyQa', ...
            '%s', record.reason)
      end
      return
   end

   record.present = true;
   record.file = filename;
   lines = readlines(filename);
   metadata_end = find(startsWith(lines, "*/"), 1, 'first');
   if isempty(metadata_end) || metadata_end == numel(lines)
      record.reason = sprintf('daily QA file has no tabular header: %s', ...
         filename);
      if strict
         error('icemodel:verification:importImau:malformedDailyQa', ...
            '%s', record.reason)
      end
      return
   end

   header = lines(metadata_end + 1);
   data_lines = lines(metadata_end + 2:end);
   first_data = find(strlength(strip(data_lines)) > 0, 1, 'first');
   has_site = any(contains(lines(1:metadata_end), "(" + site_id + ")"));
   record.valid = contains(header, "Date/Time") && ~isempty(first_data) ...
      && has_site;
   record.rows_checked = double(~isempty(first_data));
   text = strjoin(lines(1:metadata_end), newline);
   record.doi = icemodel.verification.setup.regexpOnce(text, ...
      'Citation:.*?https://doi.org/([0-9.]+/PANGAEA\.[0-9]+)');
   record.landing_url = "https://doi.org/" + string(record.doi);
   if ~record.valid
      record.reason = sprintf('daily QA file failed header/site/sample check: %s', ...
         filename);
      if strict
         error('icemodel:verification:importImau:malformedDailyQa', ...
            '%s', record.reason)
      end
   end
end

function filename = locateDailyFile(source_dir, site_id)
   %LOCATEDAILYFILE Find the daily QA file for one IMAU hourly site.
   patterns = [ ...
      fullfile(source_dir, 'daily', "GRL_" + site_id + "_AWS.tab")
      fullfile(source_dir, "GRL_" + site_id + "_AWS.tab")
      fullfile(source_dir, '**', "GRL_" + site_id + "_AWS.tab")];
   filename = "";
   for pattern = reshape(patterns, 1, [])
      hits = dir(pattern);
      hits = hits(~[hits.isdir]);
      if ~isempty(hits)
         filename = string(fullfile(hits(1).folder, hits(1).name));
         return
      end
   end
end
