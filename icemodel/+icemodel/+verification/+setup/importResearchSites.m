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
   %
   % Role
   %  Setup/update tooling. Creates or refreshes research_site manifest entries
   %  and, when not a dry run, data-only observations.mat targets.
   %
   % See also: icemodel.verification.setup.researchSiteMetadata,
   %  icemodel.verification.setup.buildSumupObservations,
   %  icemodel.verification.setup.stageRcmForcing

   arguments
      source_dir (1, 1) string = ""
      kwargs.site_ids (1, :) string = "humphrey"
      kwargs.radius_km (1, 1) double {mustBePositive} = 7.5
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.dry_run (1, 1) logical = false
      kwargs.build_forcing (1, 1) logical = false
   end

   sites = icemodel.verification.setup.researchSiteMetadata(kwargs.site_ids);
   [eval_root, input_root] = resolveRoots(kwargs);
   family_root = fullfile(eval_root, "research_site");
   manifest_file = fullfile(family_root, "manifest.json");
   source_dir = resolveCacheDir(source_dir);

   % Resolve the optional observation/forcing window once so all staged sites use
   % the same explicit clamp, or all-available records when bounds are omitted.
   [obs_start, obs_end, period] = observationWindow(kwargs.startdate, kwargs.enddate);

   entries = repmat(emptyEntry(), 1, numel(sites));
   for k = 1:numel(sites)
      entries(k) = siteEntry(sites(k), source_dir, family_root, ...
         obs_start, obs_end, period, kwargs);
   end

   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "research_site", "", "", "sumup-derived-research-site", ...
      string(datetime('today')), entries);

   if kwargs.dry_run
      return
   end

   % Write the observation manifest first; optional RCM forcing updates it using
   % the existing manifest-mode stageRcmForcing path.
   icemodel.helpers.ensureDirExists(family_root);
   manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
      manifest_file, manifest, requested_ids=[entries.case_id], ...
      overwrite_family=kwargs.overwrite_family);

   if kwargs.build_forcing
      manifest = icemodel.verification.setup.stageRcmForcing( ...
         obs_manifest=manifest, manifest_file=manifest_file, ...
         met_outdir=fullfile(input_root, 'met'), ...
         userdata_outdir=fullfile(input_root, 'userdata'), ...
         mar_dir=kwargs.mar_dir, merra_dir=kwargs.merra_dir, ...
         racmo_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir);
   end
end

function entry = siteEntry(site, source_dir, family_root, obs_start, obs_end, ...
      period, kwargs)
   %SITEENTRY Build and optionally stage one research_site case entry.
   [x3413, y3413] = projectedCoords(site);
   case_id = string(site.case_id);
   site_location = struct('lat_wgs84', site.lat_wgs84, ...
      'lon_wgs84', site.lon_wgs84, 'x_epsg3413', x3413, ...
      'y_epsg3413', y3413, 'elev_m', site.elev_m);
   evaluation_file = "";
   comparison_vars = strings(0, 1);
   obs_vars = icemodel.verification.setup.metadataStruct({ ...
      'source', 'SUMup observations selected around research-site metadata'});

   if ~kwargs.dry_run
      case_root = fullfile(family_root, case_id);
      icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);
      [observations, obs_meta] = ...
         icemodel.verification.setup.buildSumupObservations( ...
         [site.lat_wgs84, site.lon_wgs84], source_dir=source_dir, ...
         radius_km=kwargs.radius_km, startdate=obs_start, enddate=obs_end);
      targets = struct('format', 'subsurface_profile_bundle', ...
         'data', observations, 'metadata', obs_meta);
      save(fullfile(case_root, "observations.mat"), 'targets');
      evaluation_file = fullfile(case_id, "observations.mat");
      comparison_vars = comparisonVariables(observations);
      obs_vars = observationVariables();
   end

   colocation = struct();
   colocation.research_site = struct('kind', 'research_site_anchor', ...
      'staged', false, 'met_files', strings(1, 0), ...
      'note', 'no native research_site station met staged');
   colocation.sumup = struct('kind', 'firn_profile_obs', ...
      'staged', ~kwargs.dry_run, 'obs_file', char(evaluation_file), ...
      'selection_radius_km', kwargs.radius_km);
   colocation.nearest_noncolocated_promice = ...
      nearestPromice(site_location, kwargs);

   values = { ...
      char(case_id)
      'firn_observational'
      char(site.site_id)
      char(site.site_name)
      char(site.surface_zone)
      cellstr(string(site.eval_target))
      char(site.permafrost_zone)
      site_location
      period
      char(evaluation_file)
      {}
      {'sumup_obs'}
      cellstr(comparison_vars)
      obs_vars
      colocation
      'irregular'
      char(site.note)};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function [eval_root, input_root] = resolveRoots(kwargs)
   %RESOLVEROOTS Match importPromiceSites output_root semantics.
   if kwargs.output_root ~= ""
      eval_root = fullfile(kwargs.output_root, 'eval');
      input_root = fullfile(kwargs.output_root, 'input');
      return
   end
   eval_root = icemodel.verification.helpers.evaluationDataRoot( ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   input_root = icemodel.verification.helpers.inputDataRoot( ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
end

function source_dir = resolveCacheDir(source_dir)
   %RESOLVECACHEDIR Resolve the SUMup cache dir used for research observations.
   if strlength(source_dir) == 0
      source_dir = string(fullfile(icemodel.getpath('data'), ...
         'verification', 'sumup'));
   end
end

function [obs_start, obs_end, period] = observationWindow(startdate, enddate)
   %OBSERVATIONWINDOW Validate paired bounds and build manifest period strings.
   has_start = ~strcmp(string(startdate), "");
   has_end = ~strcmp(string(enddate), "");
   if has_start ~= has_end
      error('icemodel:verification:importResearchSites:halfWindow', ...
         'startdate and enddate must be provided together')
   end
   if has_start
      obs_start = icemodel.verification.setup.ensureUtc(startdate);
      obs_end = icemodel.verification.setup.ensureUtc(enddate);
      period = struct('start', char(string(obs_start)), ...
         'end', char(string(obs_end)));
   else
      obs_start = "";
      obs_end = "";
      period = struct('start', '', 'end', '');
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
   rec = struct('candidate_sites', {{'CP1', 'SWC', 'JAR'}}, ...
      'nearest_anchor', "", 'nearest_source_id', "", ...
      'distance_km', NaN, 'note', "no staged CP1/SWC/JAR promice anchors available");
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
      anchors=anchors, threshold_km=eps);
   rec.nearest_anchor = string(anchor.site);
   rec.nearest_source_id = string(anchor.source_id);
   rec.distance_km = distance_km;
   rec.note = "nearest staged CP1/SWC/JAR PROMICE anchor; not treated as native met";
end

function vars = comparisonVariables(observations)
   %COMPARISONVARIABLES List SUMup observation axes that were actually staged.
   candidate = ["density"; "subsurface_temperature"; "smb"];
   present = false(numel(candidate), 1);
   for k = 1:numel(candidate)
      present(k) = isfield(observations, candidate(k)) ...
         && ~isempty(observations.(candidate(k)));
   end
   vars = candidate(present);
end

function obs_vars = observationVariables()
   %OBSERVATIONVARIABLES Manifest metadata for SUMup-derived research targets.
   obs_vars = icemodel.verification.setup.metadataStruct({ ...
      'density', 'SUMup density profile'
      'subsurface_temperature', 'SUMup subsurface temperature profile'
      'smb', 'SUMup surface mass balance records'});
end
