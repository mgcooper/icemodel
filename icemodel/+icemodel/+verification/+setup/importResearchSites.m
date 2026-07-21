function manifest = importResearchSites(source_dir, kwargs)
   %IMPORTRESEARCHSITES Stage generic research-site firn targets.
   %
   %  manifest = icemodel.verification.setup.importResearchSites(source_dir)
   %  manifest = icemodel.verification.setup.importResearchSites(source_dir, ...
   %     case_ids="humphrey", dry_run=true)
   %
   %  Role
   %    Setup/update tooling. Creates or refreshes research_site manifest entries
   %    and, when not a dry run, data-only observations.mat targets. Dry runs are
   %    metadata-only and do not read source caches or write staged artifacts.
   %    forcing_sources selects RCM sources requested by the current call.
   %    Ordinary calls preserve omitted existing legs; overwrite_family=true
   %    deliberately replaces the whole family state.
   %    build_observations=false is a guarded non-dry fast path: requested cases
   %    must already exist in the target manifest, whose SUMup-derived observation
   %    entry is reused while selected RCM forcing is attached.
   %    Delegated RCM model met defaults to dt_out="15m"; pass dt_out="" for
   %    native model-met cadence. RCM Data/userdata defaults to hourly.
   %
   %  Default roots
   %    source_dir="" reads <repo>/data/verification/sumup. With no output_root,
   %    observations go to
   %    <repo>/data/eval/research_site/<case_id>/observations.mat and RCM
   %    met/userdata go to <repo>/data/input/{met,userdata}/<source>/.
   %
   %  Minimal first implementation: Humphrey is represented as a generic
   %  research_site anchor, with observations sourced from nearby SUMup records
   %  via buildSumupObservations. It records no native station met. Optional RCM
   %  forcing/Data legs are delegated to stageRcmForcing after observations stage.
   %
   %  Name-value
   %    case_ids : string vector. Research-site cases to stage (default
   %        "humphrey"). Pass an explicit list to stage a subset.
   %    forcing_sources : string vector subset of ["mar","merra","racmo"]
   %        (default all three). SUMup observations are always the case definition
   %        when build_observations is true; forcing_sources selects only runtime
   %        met/userdata artifacts. It is a patch selector, not the complete
   %        desired source state: an existing case's omitted legs remain
   %        unchanged during ordinary merge updates and are removed only by
   %        explicit family replacement.
   %    startdate, enddate : datetime / string. Optional explicit observation and
   %        forcing window; pass both or neither. With both omitted, a fresh
   %        observation import uses the source-authored case period and a forcing-
   %        only update reuses the staged case period. Every selected leg is
   %        limited to the overlap between this case window and its own on-disk
   %        coverage.
   %    family, observation_source : fixed "research_site" and "sumup" source
   %        contract. These options document rather than vary the importer type.
   %    radius_km : double (default 7.5). SUMup point-selection radius.
   %    mar_dir, merra_dir, racmo_dir, modis_dir : raw-source directories for
   %        the gridded source builders.
   %    output_root : string. Convenience root for eval and input artifacts.
   %        When set, eval artifacts go to <output_root>/eval and runtime files
   %        go to <output_root>/input.
   %    evaluation_data_root, input_data_root, icemodel_config_casename : lower-
   %        level root resolution when output_root is unset. With all roots
   %        unset, the importer targets the gitignored top-level <repo>/data tree.
   %    dt_out : model-met output timestep (default "15m") for gridded forcing;
   %        pass "" to retain native model-met cadence. Data/userdata uses the
   %        shared hourly writer default.
   %    overwrite : logical (default false). Refresh requested case artifacts;
   %        other cases are never touched.
   %    overwrite_family : logical (default false). Replace the family manifest
   %        from the requested cases alone. The default is merge.
   %    skip_missing : logical (default true). Record source-local skips in
   %        manifest.skipped instead of aborting.
   %    dry_run : logical (default false). Return the manifest shape without
   %        reading source caches, writing artifacts, or merging the manifest.
   %    build_forcing : logical (default false). Stage runtime artifacts for
   %        forcing_sources. When false, only observations are staged.
   %    build_observations : logical (default true). When false, requested cases
   %        must already exist in the target manifest. Their observation entries
   %        are reused while forcing_sources are attached without rebuilding the
   %        observation artifact.
   %
   %  Incremental staging (MERGE by default)
   %    Staging one case adds or updates only that case in the family manifest
   %    and preserves every other committed case and file. Re-staging the same
   %    case updates exactly its entry. Set overwrite_family=true only to
   %    deliberately rebuild the family root.
   %
   %  Returns
   %    manifest : struct  Final or dry-run family manifest.
   %
   %  See also: icemodel.verification.setup.researchSiteCatalog,
   %    icemodel.verification.setup.buildSumupObservations,
   %    icemodel.verification.setup.stageRcmForcing

   arguments
      source_dir (1, 1) string = ""
      kwargs.case_ids (1, :) string = "humphrey"
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources)} = ...
         icemodel.verification.namelists.rcmsources()
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.family (1, 1) string ...
         {mustBeMember(kwargs.family, "research_site")} = "research_site"
      kwargs.observation_source (1, 1) string ...
         {mustBeMember(kwargs.observation_source, "sumup")} = "sumup"
      kwargs.radius_km (1, 1) double {mustBePositive} = 7.5
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
   dataset_family = kwargs.family;
   rcm_sources = intersect(forcing_sources, ...
      icemodel.verification.namelists.rcmsources(), "stable");
   build_rcm_forcing = kwargs.build_forcing && ~isempty(rcm_sources);

   cases = icemodel.verification.setup.researchSiteCatalog(kwargs.case_ids);
   requested_ids = string({cases.case_id});

   % Use one normalized observation/forcing window across every staged site.
   if window_enabled
      period = struct('start', ...
         icemodel.verification.setup.formatManifestTime(window_start), ...
         'end', icemodel.verification.setup.formatManifestTime(window_end));
   else
      % SUMup readers use blank values to request all available observations.
      window_start = "";
      window_end = "";
      period = struct('start', '', 'end', '');
   end

   % Resolve output roots and paths before raw sources. Forcing-only calls can
   % reuse the existing manifest without requiring observation caches.
   [evaluation_data_root, input_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   [family_root, manifest_file, met_outdir, userdata_outdir] = ...
      icemodel.verification.setup.datasetFamilyStagingPaths( ...
      evaluation_data_root, input_root, dataset_family);
   source_dir = icemodel.verification.setup.sumupCacheDir(source_dir);

   coverage = struct();
   reuse_sources = strings(1, 0);
   % Resolve RCM coverage only for a real requested build.
   if ~kwargs.dry_run && build_rcm_forcing
      reuse_sources = rcm_sources;
      coverage = icemodel.verification.setup.rcmSourceCoverage( ...
         rcm_sources, struct('mar', kwargs.mar_dir, ...
         'merra', kwargs.merra_dir, 'racmo', kwargs.racmo_dir));
   end

   if ~kwargs.dry_run && ~kwargs.build_observations
      % Reuse the staged case entry so an RCM-only attachment does not require
      % source caches.
      [state, alive, skipped] = ...
         icemodel.verification.setup.reuseDatasetFamilyCases( ...
         manifest_file, requested_ids, emptyState(), ...
         dataset_family=dataset_family, ...
         overwrite_family=kwargs.overwrite_family, ...
         forcing_sources=reuse_sources, coverage=coverage, ...
         startdate=kwargs.startdate, enddate=kwargs.enddate);
   else
      % Validate caches only when building observations.
      % Dry runs remain metadata-only; optional skips stay quiet while required
      % SUMup products print their retrieval guidance before failing.
      if ~kwargs.dry_run
         source_dir = icemodel.verification.setup.fetchSumup( ...
            cache_dir=source_dir, strict=~kwargs.skip_missing, ...
            silent=kwargs.skip_missing);
      end

      stage_callback = @(~, n) stageCase( ...
         cases(n), source_dir, family_root, window_start, window_end, ...
         period, kwargs, coverage, rcm_sources);

      % Stage each requested case into importer state.
      [state, alive, skipped] = ...
         icemodel.verification.setup.stageDatasetFamilyCases( ...
         1:numel(requested_ids), emptyState(), stage_callback, ...
         skip_missing=kwargs.skip_missing, ...
         warning_id="icemodel:verification:importResearchSites:caseSkipped", ...
         label_callback=@(~, n) requested_ids(n));
   end

   % Record source provenance once so dry-run and persisted manifests match.
   source_doi = "10.18739/A2M61BR5M";
   source_url = "https://nsidc.org/data/g02288";
   source_version = "sumup-derived-research-site";
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

   if kwargs.build_observations
      icemodel.helpers.ensureDirExists(family_root);
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
      'permafrost_zone', "", 'eval_source', "", ...
      'native_timestep', "", 'notes', "", 'reuse_entry', false, ...
      'dry_run', false);
end

function state = stageCase(site, source_dir, family_root, obs_start, ...
      obs_end, period, kwargs, coverage, rcm_sources)
   %STAGECASE Stage one requested case and return importer state.
   if kwargs.dry_run
      state = dryRunCase(site, period, kwargs);
      return
   end

   [period, used_site_period] = casePeriod(site, period);
   if used_site_period
      % The source-authored default period also bounds SUMup selection.
      obs_start = period.start;
      obs_end = period.end;
   end
   site_location = caseLocation(site);
   case_id = string(site.case_id);

   % Build observations before touching the output tree, so skip-missing calls
   % cannot remove an existing case or leave an empty case directory.
   [observations, obs_meta] = ...
      icemodel.verification.setup.buildSumupObservations( ...
      [site.lat_wgs84, site.lon_wgs84], source_dir=source_dir, ...
      radius_km=kwargs.radius_km, startdate=obs_start, enddate=obs_end);
   if ~icemodel.verification.setup.hasObservationRecords(observations)
      error('icemodel:verification:importResearchSites:missingObservation', ...
         'no SUMup observations found for research site %s', case_id)
   end
   comparison_variables = ...
      icemodel.verification.setup.sumupComparisonVariables(observations);

   % A future site without source-authored bounds adopts its observation span.
   if isBlankPeriod(period)
      period = observationPeriod(observations);
   end
   requested_case = struct('period', period, ...
      'site_location', site_location, 'artifact_metadata', obs_meta);
   case_root = fullfile(family_root, case_id);
   write_observation = icemodel.verification.setup.prepareCaseRoot( ...
      case_root, kwargs.overwrite, "observations.mat", requested_case);
   if write_observation
      % Repeated non-overwrite imports keep the current observation bytes.
      targets = struct('format', 'subsurface_profile_bundle', ...
         'data', observations, 'metadata', obs_meta);
      targets = icemodel.verification.setup.stampArtifactMetadata(targets);
      save(fullfile(case_root, "observations.mat"), 'targets');
   end
   evaluation_file = fullfile(case_id, "observations.mat");
   colocation = caseColocation(evaluation_file, true, ...
      nearestPromice(site_location, kwargs), kwargs.radius_km);

   leg = struct();
   if kwargs.build_forcing && ~isempty(rcm_sources)
      [t1, t2] = icemodel.verification.setup.periodBounds(period);
      leg = icemodel.verification.setup.resolveLegWindows( ...
         rcm_sources, coverage, t1, t2);
   end
   state = struct('site_id', string(site.site_id), ...
      'case_id', string(site.case_id), 'site_name', string(site.site_name), ...
      'site_location', site_location, ...
      'point', [site.lat_wgs84, site.lon_wgs84], 'period', period, ...
      'evaluation_file_rel', char(evaluation_file), 'entry', struct(), ...
      'colocation', colocation, 'leg', leg, ...
      'comparison_variables', {comparison_variables}, ...
      'observation_variables', observationVariables(), ...
      'surface_zone', string(site.surface_zone), ...
      'eval_target', {string(site.eval_target)}, ...
      'permafrost_zone', string(site.permafrost_zone), ...
      'eval_source', "sumup_obs", 'native_timestep', "irregular", ...
      'notes', string(site.note), 'reuse_entry', false, 'dry_run', false);
end

function s = dryRunCase(site, period, kwargs)
   %DRYRUNCASE Build one metadata-only case state.
   period = casePeriod(site, period);
   site_location = caseLocation(site);
   colocation = caseColocation("", false, ...
      emptyNearestPromice(), kwargs.radius_km);
   comparison_variables = ["density", "subsurface_temperature", "smb"];
   s = struct('site_id', string(site.site_id), ...
      'case_id', string(site.case_id), 'site_name', string(site.site_name), ...
      'site_location', site_location, ...
      'point', [site.lat_wgs84, site.lon_wgs84], 'period', period, ...
      'evaluation_file_rel', '', 'entry', struct(), ...
      'colocation', colocation, 'leg', struct(), ...
      'comparison_variables', {comparison_variables}, ...
      'observation_variables', observationVariables(), ...
      'surface_zone', string(site.surface_zone), ...
      'eval_target', {string(site.eval_target)}, ...
      'permafrost_zone', string(site.permafrost_zone), ...
      'eval_source', "sumup_obs", 'native_timestep', "irregular", ...
      'notes', string(site.note), 'reuse_entry', false, 'dry_run', true);
end

function [period, used_site_period] = casePeriod(site, period)
   %CASEPERIOD Apply source-authored bounds to an otherwise blank request.
   used_site_period = isBlankPeriod(period) && isfield(site, 'period') ...
      && ~isBlankPeriod(site.period);
   if used_site_period
      period = site.period;
   end
end

function site_location = caseLocation(site)
   %CASELOCATION Build the canonical location record for a research site.
   [x3413, y3413] = projectedCoords(site);
   site_location = struct('lat_wgs84', site.lat_wgs84, ...
      'lon_wgs84', site.lon_wgs84, 'x_epsg3413', x3413, ...
      'y_epsg3413', y3413, 'elev_m', site.elev_m);
end

function colocation = caseColocation(evaluation_file, sumup_staged, ...
      nearest_promice, radius_km)
   %CASECOLOCATION Build stable native, SUMup, and PROMICE provenance.
   colocation = struct();
   colocation.research_site = struct('kind', 'research_site_anchor', ...
      'staged', false, 'met_files', strings(1, 0), ...
      'note', 'no native research_site station met staged');
   colocation.sumup = struct('kind', 'firn_profile_obs', ...
      'staged', sumup_staged, 'obs_file', char(evaluation_file), ...
      'selection_radius_km', radius_km);
   colocation.nearest_noncolocated_promice = nearest_promice;
end

function tf = isBlankPeriod(period)
   %ISBLANKPERIOD True when a period carries no usable bounds.
   tf = ~isstruct(period) || ~isfield(period, 'start') ...
      || ~isfield(period, 'end') ...
      || strlength(string(period.start)) == 0 ...
      || strlength(string(period.end)) == 0;
end

function period = observationPeriod(observations)
   %OBSERVATIONPERIOD Bound a research-site case by its staged observations.
   times = NaT(0, 1, 'TimeZone', 'UTC');
   if isfield(observations, 'density') && istable(observations.density) ...
         && ismember("datetime", string(observations.density.Properties.VariableNames))
      times = [times; observations.density.datetime(:)];
   end
   if isfield(observations, 'subsurface_temperature') ...
         && istimetable(observations.subsurface_temperature)
      times = [times; observations.subsurface_temperature.Time(:)];
   end
   if isfield(observations, 'smb') && istable(observations.smb)
      names = string(observations.smb.Properties.VariableNames);
      if all(ismember(["start_date", "end_date"], names))
         times = [times; observations.smb.start_date(:); observations.smb.end_date(:)];
      end
   end
   times = times(~isnat(times));
   period = struct('start', ...
      icemodel.verification.setup.formatManifestTime(min(times)), ...
      'end', icemodel.verification.setup.formatManifestTime(max(times)));
end

function [x3413, y3413] = projectedCoords(site)
   %PROJECTEDCOORDS Use metadata coords when present, otherwise compute them.
   x3413 = site.x_epsg3413;
   y3413 = site.y_epsg3413;
   if ~isfinite(x3413) || ~isfinite(y3413)
      proj = icemodel.forcing.helpers.psnProjection();
      [x3413, y3413] = projfwd(proj, site.lat_wgs84, site.lon_wgs84);
   end
end

function rec = nearestPromice(site_location, kwargs)
   %NEARESTPROMICE Record nearest staged CP1/SWC/JAR PROMICE provenance.
   rec = emptyNearestPromice();
   anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   if isempty(anchors)
      return
   end

   candidates = ["CP1", "SWC", "JAR"];
   anchor_site = upper(string({anchors.site}));
   anchor_source = upper(string({anchors.source_id}));
   anchor_case = upper(string({anchors.case_id}));
   keep = string({anchors.family}) == "promice" ...
      & (ismember(anchor_site, candidates) ...
      | ismember(anchor_source, candidates) ...
      | ismember(anchor_case, candidates));
   anchors = anchors(keep);
   if isempty(anchors)
      return
   end

   [~, anchor, distance_km] = icemodel.verification.setup.anchorColocation( ...
      site_location.x_epsg3413, site_location.y_epsg3413, ...
      anchors=anchors, threshold_km=realmax);
   rec.nearest_anchor = string(anchor.site);
   rec.nearest_source_id = string(anchor.source_id);
   rec.distance_km = distance_km;
   rec.is_colocated = false;
   rec.note = "nearest staged CP1/SWC/JAR PROMICE anchor; not treated as native met";
end

function rec = emptyNearestPromice()
   %EMPTYNEARESTPROMICE Return source-free PROMICE provenance metadata.
   rec = struct('candidate_sites', {{'CP1', 'SWC', 'JAR'}}, ...
      'nearest_anchor', "", 'nearest_source_id', "", ...
      'distance_km', NaN, 'threshold_km', 0, ...
      'is_colocated', false, ...
      'note', "no staged CP1/SWC/JAR promice anchors available");
end

function observation_variables = observationVariables()
   %OBSERVATIONVARIABLES Return staged observation metadata.
   observation_variables = icemodel.verification.setup.metadataStruct({ ...
      'density', 'SUMup density profile'
      'subsurface_temperature', 'SUMup subsurface temperature profile'
      'smb', 'SUMup surface mass balance records'});
end
