function manifest = importImau(source_dir, kwargs)
   %IMPORTIMAU Stage the IMAU hourly AWS verification family.
   %
   %  manifest = icemodel.verification.setup.importImau(source_dir)
   %  manifest = icemodel.verification.setup.importImau(source_dir, dry_run=true)
   %
   % Role
   %  First-pass IMAU staging hook for the hourly PANGAEA S21/S22/S23 network.
   %  The daily 19-station product remains QA/provenance input, not a first-pass
   %  case inventory.
   %  forcing_sources selects the optional RCM sources requested by the current
   %  call. Ordinary calls preserve omitted existing RCM legs;
   %  overwrite_family=true is the explicit whole-family replacement boundary.
   %  build_observations=false is a guarded non-dry fast path: requested cases
   %  must already exist in the target manifest, whose observation/native entry
   %  is reused while selected RCM forcing is attached.
   %
   % Default roots
   %  source_dir="" reads <repo>/data/verification/imau. With no output_root,
   %  observations go to <repo>/data/eval/imau/<case_id>/observations.mat and
   %  native met/userdata go to <repo>/data/input/{met,userdata}/imau/.
   %  Model met defaults to dt_out="15m"; pass dt_out="" for native cadence.
   %  Data/userdata defaults to hourly at the shared writer boundary.
   %  Native met schema completion is fixed at the importer boundary: absent
   %  required channels are retained as NaN placeholders.
   %  Call buildImauHourlyMet directly for strict source-schema validation.

   arguments
      source_dir (1, 1) string = ""
      kwargs.site_ids (1, :) string = strings(1, 0)
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources)} = ...
         icemodel.verification.namelists.rcmsources()
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
      kwargs.dry_run (1, 1) logical = false
      kwargs.strict (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = true
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.build_forcing (1, 1) logical = false
      kwargs.build_observations (1, 1) logical = true
   end

   % Validate the optional clamp before any cache or staging side effect.
   icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);

   sites = icemodel.verification.setup.imauSiteCatalog(kwargs.site_ids);
   dataset_family = "imau";
   requested_ids = lower([sites.site_id]);
   rcm_sources = intersect(kwargs.forcing_sources, ...
      icemodel.verification.namelists.rcmsources(), "stable");
   coverage = struct();
   reuse_sources = strings(1, 0);

   if ~kwargs.dry_run
      % Non-dry runs write observations.mat under eval/imau and native runtime
      % files under input/met/imau and input/userdata/imau before optional RCMs.
      [evaluation_data_root, input_root] = ...
         icemodel.verification.setup.resolveStagingRoots( ...
         output_root=kwargs.output_root, ...
         evaluation_data_root=kwargs.evaluation_data_root, ...
         input_data_root=kwargs.input_data_root, ...
         icemodel_config_casename=kwargs.icemodel_config_casename);
      family_root = fullfile(evaluation_data_root, dataset_family);
      manifest_file = fullfile(family_root, "manifest.json");
      prior_cases = ...
         icemodel.verification.setup.loadPriorDatasetFamilyCases( ...
         manifest_file, build_forcing=kwargs.build_forcing, ...
         overwrite_family=kwargs.overwrite_family);

      met_outdir = fullfile(input_root, 'met');
      userdata_outdir = fullfile(input_root, 'userdata');
      if kwargs.build_observations
         icemodel.helpers.ensureDirExists(family_root);
      end
      if kwargs.build_forcing
         icemodel.helpers.ensureDirExists(met_outdir);
         icemodel.helpers.ensureDirExists(userdata_outdir);
      end

      if kwargs.build_forcing && ~isempty(rcm_sources)
         reuse_sources = rcm_sources;
         coverage = icemodel.verification.setup.promiceSourceCoverage( ...
            rcm_sources, struct('mar', kwargs.mar_dir, ...
            'merra', kwargs.merra_dir, 'racmo', kwargs.racmo_dir));
      end
   end

   if ~kwargs.dry_run && ~kwargs.build_observations
      % Reuse the staged hourly-site entry so an RCM-only attachment does not
      % require the PANGAEA hourly or daily-QA caches.
      [state, alive, skipped] = ...
         icemodel.verification.setup.reuseDatasetFamilyCases( ...
         manifest_file, requested_ids, emptyState(), ...
         forcing_sources=reuse_sources, coverage=coverage, ...
         startdate=kwargs.startdate, enddate=kwargs.enddate);
   else
      % Validate caches only when building observations/native runtime files.
      % Dry runs remain metadata-only; optional skips stay quiet while required
      % PANGAEA products print their DOI guidance before failing.
      strict_cache = kwargs.strict || (~kwargs.dry_run && ~kwargs.skip_missing);
      [source_dir, cache_status] = icemodel.verification.setup.fetchImau( ...
         cache_dir=icemodel.forcing.helpers.verificationSourceDir( ...
         source_dir, "imau"), strict=strict_cache, ...
         silent=kwargs.skip_missing, create_cache_dir=~kwargs.dry_run);
      if ~kwargs.dry_run && ~kwargs.skip_missing
         preflightHourlyFiles(source_dir, [sites.site_id]);
      end
      daily_qa = dailyQaStatus(source_dir, [sites.site_id], ...
         strict=(~kwargs.dry_run && ~kwargs.skip_missing));
      stage_callback = @(~, n) dryRunSite( ...
         sites(n), source_dir, cache_status, daily_qa);
      if ~kwargs.dry_run
         stage_callback = @(~, n) stageHourlySite( ...
            sites(n), source_dir, family_root, met_outdir, userdata_outdir, ...
            kwargs, daily_qa, coverage, rcm_sources, prior_cases);
      end

      % Stage each requested case into importer state.
      [state, alive, skipped] = ...
         icemodel.verification.setup.stageDatasetFamilyCases( ...
         1:numel(sites), emptyState(), stage_callback, ...
         skip_missing=kwargs.skip_missing, ...
         warning_id="icemodel:verification:importImau:siteSkipped", ...
         label_callback=@(~, n) sites(n).case_id);
   end

   % Record source provenance once so dry-run and persisted manifests match.
   source_doi = "10.1594/PANGAEA.971647;10.1594/PANGAEA.970127";
   source_url = "https://doi.org/10.1594/PANGAEA.971647";
   source_version = "hourly-s21-s22-s23+daily-qa";
   retrieval_date = string(datetime('today'));

   if kwargs.dry_run
      manifest = icemodel.verification.setup.runDatasetFamilyDryRun( ...
         state, alive, dataset_family=dataset_family, ...
         requested_ids=requested_ids, skipped=skipped, ...
         source_doi=source_doi, source_url=source_url, ...
         source_version=source_version, retrieval_date=retrieval_date, ...
         entry_callback=@caseEntry);
      return
   end

   leg_callback = [];
   if kwargs.build_forcing && ~isempty(rcm_sources)
      leg_callback = @(s, src) s.leg.(char(src));
   end

   [manifest, ~] = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family=dataset_family, ...
      manifest_file=manifest_file, requested_ids=requested_ids, ...
      skipped=skipped, ...
      source_doi=source_doi, source_url=source_url, ...
      source_version=source_version, retrieval_date=retrieval_date, ...
      overwrite_family=kwargs.overwrite_family, overwrite=kwargs.overwrite, ...
      entry_callback=@caseEntry, ...
      build_forcing=kwargs.build_forcing, forcing_sources=rcm_sources, ...
      leg_callback=leg_callback, met_outdir=met_outdir, ...
      userdata_outdir=userdata_outdir, mar_dir=kwargs.mar_dir, ...
      merra_dir=kwargs.merra_dir, racmo_dir=kwargs.racmo_dir, ...
      modis_dir=kwargs.modis_dir, method="nearest", dt_out=kwargs.dt_out);
end

function preflightHourlyFiles(source_dir, site_ids)
   %PREFLIGHTHOURLYFILES Fail before staging when requested IMAU files are absent.
   for site_id = reshape(site_ids, 1, [])
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

function s = emptyState()
   %EMPTYSTATE Prototype per-site staging state.
   s = struct('site_id', "", 'case_id', "", 'alias', "", ...
      'site_name', "", 'site_location', struct(), 'point', [NaN NaN], ...
      'period', struct('start', '', 'end', ''), 'evaluation_file_rel', '', ...
      'entry', struct(), 'colocation', struct(), 'leg', struct(), ...
      'comparison_vars', {strings(0, 1)}, 'obs_vars', struct(), ...
      'surface_zone', "", 'target', {strings(0, 1)}, ...
      'permafrost_zone', "", 'note', "", 'reuse_entry', false, ...
      'dry_run', false);
end

function state = stageHourlySite(site, source_dir, family_root, met_outdir, ...
      userdata_outdir, kwargs, daily_qa, coverage, rcm_sources, prior_cases)
   %STAGEHOURLYSITE Stage one IMAU hourly site and return manifest state.
   alias = lower(site.site_id);
   [met, metadata, Data] = icemodel.forcing.buildImauHourlyMet( ...
      site.site_id, source_dir=source_dir, startdate=kwargs.startdate, ...
      enddate=kwargs.enddate, fillwithmissing=true);

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

   case_root = fullfile(family_root, alias);
   write_observation = icemodel.verification.setup.prepareCaseRoot( ...
      case_root, kwargs.overwrite, "observations.mat", requested_case);

   met_files = strings(1, 0);
   data_files = strings(1, 0);
   forcing_ready = false;
   forcing_ready_reason = "build_forcing=false";
   % No-forcing calls publish no saved-artifact window diagnostics; a compatible
   % prior forcing leg can selectively restore them below.
   forcing_complete_windows = repmat(struct( ...
      'start_time', "", 'end_time', "", 'sample_count', 0), 0, 1);
   if kwargs.build_forcing
      % Native input artifacts obey the suite-wide forcing toggle; the same
      % source Data can still back observations.mat when forcing is disabled.
      met_files = icemodel.forcing.helpers.writemet(met, alias, "imau", ...
         outdir=met_outdir, naming="window", dt_out=kwargs.dt_out, ...
         overwrite=kwargs.overwrite);
      % Diagnose the exact returned path because no-overwrite staging may
      % select an existing exact or broader enclosing artifact.
      [forcing_ready, forcing_ready_reason, forcing_complete_windows] = ...
         icemodel.verification.setup.metArtifactReadiness(met_files);
      data_files = icemodel.forcing.helpers.writeuserdata(Data, alias, ...
         "imau", outdir=userdata_outdir, naming="window", ...
         overwrite=kwargs.overwrite);
   end

   if write_observation
      % A repeated non-overwrite import retains the current observation bytes.
      targets = struct('format', 'timeseries', 'data', Data, ...
         'metadata', metadata);
      targets = icemodel.verification.setup.stampArtifactMetadata(targets);
      save(fullfile(case_root, 'observations.mat'), 'targets');
   end
   evaluation_file_rel = char(fullfile(alias, 'observations.mat'));

   colocation = imauColocation(source_dir, metadata, met_files, data_files, ...
      met_outdir, userdata_outdir, evaluation_file_rel, period, ...
      site.source_association, daily_qa, forcing_ready, ...
      forcing_ready_reason, forcing_complete_windows);
   [colocation, identity_conflict] = preservePriorImauLeg( ...
      colocation, prior_cases, site.case_id, metadata, kwargs);
   if identity_conflict
      % Keep the fresh evaluation provenance but make the missing native rebuild
      % explicit; the prior files stay on disk without remaining manifest links.
      colocation.imau.kind = 'hourly_aws_eval';
      colocation.imau.forcing_ready = false;
      colocation.imau.forcing_ready_reason = ...
         'native source identity changed; rerun with build_forcing=true';
   end
   leg = struct();
   if kwargs.build_forcing && ~isempty(rcm_sources)
      leg = icemodel.verification.setup.resolveLegWindows( ...
         rcm_sources, coverage, Data.Time(1), Data.Time(end));
   end

   state = struct('site_id', site.site_id, 'case_id', site.case_id, ...
      'alias', alias, 'site_name', site.site_name, ...
      'site_location', site_location, ...
      'point', [site_location.lat_wgs84, site_location.lon_wgs84], ...
      'period', period, 'evaluation_file_rel', evaluation_file_rel, ...
      'entry', struct(), 'colocation', colocation, 'leg', leg, ...
      'comparison_vars', {imauComparisonVariables(Data)}, ...
      'obs_vars', imauObservationVariables(metadata, Data, daily_qa), ...
      'surface_zone', site.surface_zone, ...
      'target', {string(site.eval_target)}, ...
      'permafrost_zone', site.permafrost_zone, 'note', site.note, ...
      'reuse_entry', false, 'dry_run', false);
end

function [colocation, identity_conflict] = preservePriorImauLeg( ...
      colocation, prior_cases, case_id, metadata, kwargs)
   %PRESERVEPRIORIMAULEG Merge native artifacts into a fresh observation leg.
   identity_conflict = false;
   if kwargs.build_forcing || kwargs.overwrite_family || isempty(prior_cases)
      return
   end
   hit = find(string({prior_cases.case_id}) == string(case_id), 1);
   if isempty(hit) || ~isfield(prior_cases(hit), 'colocation') ...
         || ~isfield(prior_cases(hit).colocation, 'imau')
      return
   end
   prior = prior_cases(hit).colocation.imau;

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
   if isfield(prior_cases(hit), 'observation_variables') ...
         && isstruct(prior_cases(hit).observation_variables)
      prior_obs = prior_cases(hit).observation_variables;
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
   if isfield(prior_cases(hit), 'site_location') ...
         && isstruct(prior_cases(hit).site_location) ...
         && all(isfield(prior_cases(hit).site_location, ...
         ["lat_wgs84", "lon_wgs84"]))
      n_points = n_points + 1;
      prior_points(n_points, :) = [ ...
         prior_cases(hit).site_location.lat_wgs84, ...
         prior_cases(hit).site_location.lon_wgs84];
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

function state = dryRunSite(site, source_dir, cache_status, daily_qa)
   %DRYRUNSITE Build one non-writing IMAU state record.
   imau = struct('kind', 'hourly_aws', 'staged', false, ...
      'source_dir', char(source_dir), 'cache_status', cache_status);
   colocation = struct('imau', imau, ...
      'source_association', site.source_association, ...
      'daily_qa', daily_qa);
   state = struct('site_id', site.site_id, 'case_id', site.case_id, ...
      'alias', lower(site.site_id), 'site_name', site.site_name, ...
      'site_location', struct('lat_wgs84', site.lat_wgs84, ...
      'lon_wgs84', site.lon_wgs84, 'x_epsg3413', site.x_epsg3413, ...
      'y_epsg3413', site.y_epsg3413, 'elev_m', site.elev_m), ...
      'point', [site.lat_wgs84, site.lon_wgs84], ...
      'period', struct('start', '', 'end', ''), ...
      'evaluation_file_rel', '', 'entry', struct(), ...
      'colocation', colocation, 'leg', struct(), ...
      'comparison_vars', {["tair", "rh", "wspd", "wdir", "psfc", ...
      "swd", "swu", "lwd", "lwu", "albedo", "tsfc", "surface_height"]}, ...
      'obs_vars', struct('hourly_site', site.site_id, ...
      'daily_qa_bundle', daily_qa.bundle_doi), ...
      'surface_zone', site.surface_zone, ...
      'target', {string(site.eval_target)}, ...
      'permafrost_zone', site.permafrost_zone, 'note', site.note, ...
      'reuse_entry', false, 'dry_run', true);
end

function entry = caseEntry(s)
   %CASEENTRY Convert one staged IMAU state record to a manifest entry.
   if s.reuse_entry
      % Fast-path entries retain all family-specific metadata; only requested
      % colocation legs and their derived source lists may change.
      entry = s.entry;
      entry.colocation = s.colocation;
      [forcing_sources, eval_sources] = ...
         icemodel.verification.setup.colocationSourceLists(entry.colocation);
      entry.forcing_sources = cellstr(forcing_sources);
      entry.eval_sources = cellstr(eval_sources);
      return
   end
   [forcing_sources, eval_sources] = ...
      icemodel.verification.setup.colocationSourceLists(s.colocation);
   if isfield(s, 'dry_run') && s.dry_run && isempty(eval_sources)
      eval_sources = "imau_obs";
   end

   values = { ...
      char(s.case_id)
      'firn_observational'
      char(s.site_id)
      char(s.site_name)
      char(s.surface_zone)
      cellstr(s.target)
      char(s.permafrost_zone)
      s.site_location
      s.period
      s.evaluation_file_rel
      cellstr(forcing_sources(:))
      cellstr(eval_sources(:))
      cellstr(s.comparison_vars(:))
      s.obs_vars
      s.colocation
      '1hr'
      sprintf('IMAU hourly AWS site %s. %s', s.site_id, s.note)};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
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
      kwargs.strict (1, 1) logical = false
   end

   records = repmat(dailyQaRecord(), 1, numel(site_ids));
   for k = 1:numel(site_ids)
      records(k) = inspectDailyQaFile(source_dir, site_ids(k), kwargs.strict);
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
