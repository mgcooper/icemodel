function manifest = importResearchSites(source_dir, kwargs)
   %IMPORTRESEARCHSITES Stage generic research-site firn targets.
   %
   %  manifest = icemodel.verification.setup.importResearchSites(source_dir)
   %  manifest = icemodel.verification.setup.importResearchSites(source_dir, ...
   %     site_ids="humphrey", dry_run=true)
   %
   % Minimal first implementation: Humphrey is represented as a generic
   % research_site anchor, with observations sourced from nearby SUMup records
   % via buildSumupObservations. It records no native station met. Optional RCM
   % forcing/Data legs are delegated to stageRcmForcing after observations stage.
   % When source_dir is blank, SUMup observations are read from
   % <repo>/data/verification/sumup and outputs stage under top-level
   % <repo>/data/{eval,input} unless roots are overridden.
   %
   % Role
   %  Setup/update tooling. Creates or refreshes research_site manifest entries
   %  and, when not a dry run, data-only observations.mat targets. Dry runs are
   %  metadata-only previews and do not resolve or read source/staging roots.
   %  forcing_sources selects RCM sources requested by the current call.
   %  Ordinary calls preserve omitted existing legs; overwrite_family=true
   %  deliberately replaces the whole family state.
   %  build_observations=false is a guarded non-dry fast path: requested cases
   %  must already exist in the target manifest, whose SUMup-derived observation
   %  entry is reused while selected RCM forcing is attached.
   %  Delegated RCM model met defaults to dt_out="15m"; pass dt_out="" for
   %  native model-met cadence. RCM Data/userdata defaults to hourly.
   %
   % See also: icemodel.verification.setup.researchSiteCatalog,
   %  icemodel.verification.setup.buildSumupObservations,
   %  icemodel.verification.setup.stageRcmForcing

   arguments
      source_dir (1, 1) string = ""
      kwargs.family (1, 1) string ...
         {mustBeMember(kwargs.family, "research_site")} = "research_site"
      kwargs.site_ids (1, :) string = "humphrey"
      kwargs.observation_source (1, 1) string ...
         {mustBeMember(kwargs.observation_source, "sumup")} = "sumup"
      kwargs.radius_km (1, 1) double {mustBePositive} = 7.5
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
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = true
      kwargs.dry_run (1, 1) logical = false
      kwargs.build_forcing (1, 1) logical = false
      kwargs.build_observations (1, 1) logical = true
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "15m"])} = "15m"
   end

   % Reject malformed windows before metadata/cache/root resolution.
   [window_start, window_end, window_enabled] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);

   sites = icemodel.verification.setup.researchSiteCatalog(kwargs.site_ids);
   dataset_family = kwargs.family;
   requested_ids = string({sites.case_id});

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

   % Record source provenance once so dry-run and persisted manifests match.
   source_version = "sumup-derived-research-site";
   retrieval_date = string(datetime('today'));

   % Build a metadata-only preview before resolving any source or staging root.
   % This keeps dry runs valid on a clean machine while preserving the complete
   % firn manifest-entry shape used by the non-dry importer.
   if kwargs.dry_run
      [state, alive, skipped] = ...
         icemodel.verification.setup.stageDatasetFamilyCases( ...
         1:numel(sites), emptyState(), ...
         @(~, n) dryRunResearchSite(sites(n), period, kwargs), ...
         skip_missing=kwargs.skip_missing, ...
         warning_id="icemodel:verification:importResearchSites:siteSkipped", ...
         label_callback=@(~, n) sites(n).case_id);
      manifest = icemodel.verification.setup.runDatasetFamilyDryRun( ...
         state, alive, dataset_family=dataset_family, ...
         requested_ids=requested_ids, skipped=skipped, ...
         source_version=source_version, retrieval_date=retrieval_date, ...
         entry_callback=@icemodel.verification.setup.stateCaseEntry);
      return
   end

   % Non-dry calls resolve the source and destination trees only after the
   % source-free preview path has returned.
   [evaluation_data_root, input_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   family_root = fullfile(evaluation_data_root, dataset_family);
   manifest_file = fullfile(family_root, "manifest.json");
   source_dir = icemodel.verification.setup.sumupCacheDir(source_dir);

   rcm_sources = kwargs.forcing_sources;
   coverage = struct();
   reuse_sources = strings(1, 0);
   if kwargs.build_forcing && ~isempty(rcm_sources)
      reuse_sources = rcm_sources;
      coverage = icemodel.verification.setup.promiceSourceCoverage( ...
         rcm_sources, struct('mar', kwargs.mar_dir, ...
         'merra', kwargs.merra_dir, 'racmo', kwargs.racmo_dir));
   end

   if ~kwargs.build_observations
      % A forcing-only refresh reuses the staged SUMup-derived observation entry
      % and does not require the SUMup cache again.
      [state, alive, skipped] = ...
         icemodel.verification.setup.reuseDatasetFamilyCases( ...
         manifest_file, requested_ids, emptyState(), ...
         forcing_sources=reuse_sources, coverage=coverage, ...
         startdate=kwargs.startdate, enddate=kwargs.enddate);
   else
      % Research-site observation builds are backed by SUMup, so validate the
      % cache before invoking any source-backed case callback.
      source_dir = icemodel.verification.setup.fetchSumup( ...
         cache_dir=source_dir, strict=~kwargs.skip_missing, ...
         silent=kwargs.skip_missing);

      % Stage each requested case into importer state.
      [state, alive, skipped] = ...
         icemodel.verification.setup.stageDatasetFamilyCases( ...
         1:numel(sites), emptyState(), ...
         @(~, n) stageResearchSite(sites(n), source_dir, family_root, ...
         window_start, window_end, period, kwargs, coverage, rcm_sources), ...
         skip_missing=kwargs.skip_missing, ...
         warning_id="icemodel:verification:importResearchSites:siteSkipped", ...
         label_callback=@(~, n) sites(n).case_id);
   end

   icemodel.helpers.ensureDirExists(family_root);
   leg_callback = [];
   if kwargs.build_forcing && ~isempty(rcm_sources)
      leg_callback = @(s, src) s.leg.(char(src));
   end
   [manifest, ~] = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family=dataset_family, ...
      manifest_file=manifest_file, requested_ids=requested_ids, ...
      skipped=skipped, ...
      source_version=source_version, retrieval_date=retrieval_date, ...
      overwrite_family=kwargs.overwrite_family, overwrite=kwargs.overwrite, ...
      entry_callback=@icemodel.verification.setup.stateCaseEntry, ...
      build_forcing=kwargs.build_forcing, forcing_sources=rcm_sources, ...
      leg_callback=leg_callback, met_outdir=fullfile(input_root, 'met'), ...
      userdata_outdir=fullfile(input_root, 'userdata'), ...
      mar_dir=kwargs.mar_dir, merra_dir=kwargs.merra_dir, ...
      racmo_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
      dt_out=kwargs.dt_out);
end

function s = emptyState()
   %EMPTYSTATE Prototype research-site staging state.
   s = struct('case_id', "", 'alias', "", 'point', [NaN NaN], ...
      'entry', emptyEntry(), 'period', struct('start', '', 'end', ''), ...
      'colocation', struct(), 'leg', struct(), 'reuse_entry', false, ...
      'dry_run', false);
end

function s = dryRunResearchSite(site, period, kwargs)
   %DRYRUNRESEARCHSITE Build one source-free research-site preview state.
   period = sitePeriod(site, period);
   site_location = siteLocation(site);
   colocation = researchSiteColocation("", false, ...
      emptyNearestPromice(), kwargs.radius_km);
   comparison_vars = ["density", "subsurface_temperature", "smb"];
   entry = siteEntry(site, site_location, period, "", comparison_vars, ...
      colocation, kwargs);
   s = researchSiteState(entry, struct(), true);
end

function s = stageResearchSite(site, source_dir, family_root, obs_start, ...
      obs_end, period, kwargs, coverage, rcm_sources)
   %STAGERESEARCHSITE Stage one source-backed research_site case.
   [period, used_site_period] = sitePeriod(site, period);
   if used_site_period
      % The source-authored default period also bounds SUMup selection.
      obs_start = period.start;
      obs_end = period.end;
   end
   site_location = siteLocation(site);
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
   comparison_vars = ...
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
   colocation = researchSiteColocation(evaluation_file, true, ...
      nearestPromice(site_location, kwargs), kwargs.radius_km);
   entry = siteEntry(site, site_location, period, evaluation_file, ...
      comparison_vars, colocation, kwargs);

   leg = struct();
   if kwargs.build_forcing && ~isempty(rcm_sources)
      [t1, t2] = icemodel.verification.setup.periodBounds(entry.period);
      leg = icemodel.verification.setup.resolveLegWindows( ...
         rcm_sources, coverage, t1, t2);
   end
   s = researchSiteState(entry, leg, false);
end

function s = researchSiteState(entry, leg, dry_run)
   %RESEARCHSITESTATE Wrap one manifest entry for shared importer orchestration.
   s = struct('case_id', string(entry.case_id), ...
      'alias', string(entry.case_id), ...
      'point', [entry.site_location.lat_wgs84, ...
      entry.site_location.lon_wgs84], ...
      'entry', entry, 'period', entry.period, ...
      'colocation', entry.colocation, 'leg', leg, 'reuse_entry', false, ...
      'dry_run', dry_run);
end

function entry = siteEntry(site, site_location, period, evaluation_file, ...
      comparison_vars, colocation, kwargs)
   %SITEENTRY Build one canonical research_site manifest entry.
   [forcing_sources, eval_sources] = sourceListsForSite(colocation, kwargs);
   obs_vars = observationVariables();
   values = { ...
      char(site.case_id)
      'firn_observational'
      char(site.site_id)
      char(site.site_name)
      char(site.surface_zone)
      cellstr(string(site.eval_target))
      char(site.permafrost_zone)
      site_location
      period
      char(evaluation_file)
      cellstr(forcing_sources)
      cellstr(eval_sources)
      cellstr(comparison_vars)
      obs_vars
      colocation
      'irregular'
      char(site.note)};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function [period, used_site_period] = sitePeriod(site, period)
   %SITEPERIOD Apply source-authored bounds to an otherwise blank request.
   used_site_period = isBlankPeriod(period) && isfield(site, 'period') ...
      && ~isBlankPeriod(site.period);
   if used_site_period
      period = site.period;
   end
end

function site_location = siteLocation(site)
   %SITELOCATION Build the canonical location record for a research site.
   [x3413, y3413] = projectedCoords(site);
   site_location = struct('lat_wgs84', site.lat_wgs84, ...
      'lon_wgs84', site.lon_wgs84, 'x_epsg3413', x3413, ...
      'y_epsg3413', y3413, 'elev_m', site.elev_m);
end

function colocation = researchSiteColocation(evaluation_file, sumup_staged, ...
      nearest_promice, radius_km)
   %RESEARCHSITECOLOCATION Build stable native, SUMup, and PROMICE provenance.
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

function [forcing_sources, eval_sources] = sourceListsForSite(colocation, kwargs)
   %SOURCELISTSFORSITE Derive source lists while preserving dry-run intent.
   if kwargs.dry_run
      forcing_sources = strings(0, 1);
      eval_sources = "sumup_obs";
   else
      [forcing_sources, eval_sources] = ...
         icemodel.verification.setup.colocationSourceLists(colocation);
   end
end

function entry = emptyEntry()
   %EMPTYENTRY Prototype research_site manifest entry.
   values = { ...
      ''
      'firn_observational'
      ''
      ''
      'unknown'
      {'firn'}
      'unknown'
      struct('lat_wgs84', NaN, 'lon_wgs84', NaN, ...
      'x_epsg3413', NaN, 'y_epsg3413', NaN, 'elev_m', NaN)
      struct('start', '', 'end', '')
      ''
      {}
      {'sumup_obs'}
      {}
      struct()
      struct()
      'irregular'
      ''};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
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

function obs_vars = observationVariables()
   %OBSERVATIONVARIABLES Manifest metadata for SUMup-derived research targets.
   obs_vars = icemodel.verification.setup.metadataStruct({ ...
      'density', 'SUMup density profile'
      'subsurface_temperature', 'SUMup subsurface temperature profile'
      'smb', 'SUMup surface mass balance records'});
end
