function manifest = importSumup(source_dir, kwargs)
   %IMPORTSUMUP Stage co-located SUMup firn evaluation cases.
   %
   %  manifest = icemodel.verification.setup.importSumup(source_dir)
   %  manifest = icemodel.verification.setup.importSumup(source_dir, ...
   %     points=[lat lon; ...], overwrite=true)
   %  manifest = icemodel.verification.setup.importSumup(source_dir, ...
   %     anchors=anchor_catalog)
   %
   %  Role
   %    Setup/update tooling. Creates or refreshes staged data under
   %    data/eval/sumup by default and is not part of normal verification runs.
   %
   %  Default roots
   %    source_dir="" reads <repo>/data/verification/sumup. With no output_root,
   %    observations go to <repo>/data/eval/sumup/<case_id>/observations.mat and
   %    RCM met/userdata go to <repo>/data/input/{met,userdata}/<source>/.
   %
   %  Stages SUMup firn observation cases under
   %  data/eval/sumup/<case_id>/ by default. Observation import is DECOUPLED
   %  from RCM forcing: for each selected SUMup point it
   %    - reads the nearest SUMup firn observations (density / subsurface
   %      temperature / SMB) via buildSumupObservations and stages them
   %      once as the bundled observations.mat eval target;
   %    - records whether the point is within a co-location threshold of a
   %      mixed family anchor (helpers.sumupColocation, default 7.5 km EPSG:3413)
   %      without rebuilding the anchor family's native files;
   %    - writes a FORCING-AGNOSTIC family manifest.json
   %      (case_type="firn_observational");
   %  then, when build_forcing is true (the convenience bundle), delegates the
   %  co-located gridded forcing to icemodel.verification.setup.stageRcmForcing,
   %  which stages MAR + MERRA met+Data and RACMO Data at each point (SUMup points
   %  carry no station met). The forcing lives in separate, runtime-discoverable
   %  met/userdata files recorded in the manifest colocation record by source id,
   %  never bundled. Setting build_forcing=false imports observations only; the
   %  RCM forcing can be built later, independently, on the staged manifest.
   %  forcing_sources selects RCM sources requested by the current call.
   %  Ordinary calls preserve omitted existing legs; overwrite_family=true
   %  deliberately replaces the whole family state.
   %
   %  Source: SUMup_2025 is read from the LOCAL verification cache
   %  (data/verification/sumup, NSIDC G02288); the files are committed there, so
   %  no download is required. fetchSumup verifies cache presence (and prints a
   %  retrieval banner only if it is genuinely absent). Provide points explicitly
   %  (SUMup is a point collection, not a curated site list); the default is the
   %  PROMICE anchor transect so the staged SUMup cases co-locate with the
   %  existing firn/promice bundles.
   %
   %  Inputs
   %    source_dir : string  SUMup cache dir (see fetchSumup). When blank, the
   %                 default data/verification/sumup is used.
   %
   %  Name-value
   %    points : Nx2 double  [lat lon] WGS84 points to stage. Default: the
   %             explicit anchors when anchors is provided, otherwise staged
   %             mixed-anchor coordinates from firn family manifests. The fully
   %             implicit default applies the audited SUMup 2025 canonical case
   %             map; explicit points/case_ids/anchors are never mapped.
   %    anchors : struct array  mixed anchor catalog with site/family/source_id
   %             and coordinates. When points is empty, anchors with finite
   %             lat/lon become the staging points. When points is provided,
   %             anchors are used only for nearest-anchor metadata.
   %    case_ids : string vector  case ids for each point (default sumup_NN).
   %    years : numeric vector  forcing-window years for the co-located MAR/MERRA
   %            legs (default 2012:2018, the RACMO subsurface span; RACMO itself
   %            is decoupled and always stages its full on-disk coverage). Used
   %            only when startdate/enddate are omitted.
   %    startdate / enddate : OPTIONAL comparison-window clamp; pass both or
   %            neither. The DEFAULT (omitted) is ALL AVAILABLE: no clamp, each
   %            point stages its full on-disk SUMup record and the manifest
   %            period records "" / "". There is no buried 2012-2015
   %            comparison window.
   %    radius_km : double (default 7.5)  SUMup point-selection radius.
   %    colocation_threshold_km : double (default 7.5)  PROMICE-anchor
   %            co-location threshold (see helpers.sumupColocation).
   %    mar_dir / merra_dir / racmo_dir / modis_dir : string source dirs for the
   %            co-located RCM builders (delegated to stageRcmForcing).
   %    dt_out : model-met output timestep (default "15m"). Pass "" only to
   %            preserve an explicit native model-met cadence. Data/userdata
   %            uses the shared hourly writer default.
   %    forcing_sources : string vector  RCM sources requested by this call when
   %            build_forcing=true. Defaults to
   %            icemodel.verification.namelists.rcmsources().
   %    build_forcing : logical (default false). When true the co-located
   %            MAR/MERRA met+Data and RACMO Data are staged after the SUMup
   %            observation import; when false ONLY the observations + manifest
   %            are written (forcing can be built later via stageRcmForcing).
   %    build_observations : logical (default true). When false, explicit
   %            case_ids must already exist in the target manifest. Their staged
   %            observation entries are reused while forcing_sources are attached
   %            without reading the SUMup cache or rediscovering points; reuse and
   %            fresh staging share the same manifest finalization path.
   %    output_root : base output root. When set, eval goes to
   %            <output_root>/eval and runtime files to <output_root>/input.
   %    evaluation_data_root / input_data_root / icemodel_config_casename :
   %            lower-level staging roots. DEFAULT (roots unset) stages to the
   %            gitignored top-level <repo>/data tree. Pass
   %            output_root/evaluation_data_root/input_data_root or a
   %            nonblank config casename only when intentionally targeting
   %            another tree.
   %    overwrite : logical (default false)  refresh a requested point's own
   %            staged case folder; other cases are never touched.
   %    overwrite_family : logical (default false)  force a FULL rewrite of the
   %            family manifest from the requested points alone. DEFAULT MERGE.
   %    skip_missing : logical (default true)  record data-gated points in
   %            manifest.skipped rather than aborting.
   %    dry_run : logical (default false)  return the manifest shape without
   %            writing observations, runtime files, or the family manifest.
   %
   %  Incremental staging (MERGE by default)
   %    Staging a point ADDS or UPDATES only that case in the sumup family
   %    manifest and PRESERVES every other committed case + files byte for byte
   %    through icemodel.verification.setup.writeFamilyManifestMerge.
   %    overwrite_family=true rebuilds the family root.
   %
   %  Colocation architecture: SUMup observations are staged first as the
   %    bundled data-only observations.mat eval target. Co-located RCM forcing
   %    legs (MAR + MERRA met+Data, RACMO Data) are staged as individual runtime
   %    files recorded in the manifest colocation record by source id. A failing
   %    RCM source degrades to a skipped forcing leg, never a skipped SUMup
   %    point.
   %
   %  Returns
   %    manifest : struct  Family manifest also written to manifest.json.
   %
   %  See also: icemodel.verification.setup.fetchSumup,
   %    icemodel.verification.setup.buildSumupObservations,
   %    icemodel.verification.setup.stageRcmForcing,
   %    icemodel.verification.setup.importPromiceSites,
   %    icemodel.verification.helpers.sumupColocation

   arguments
      source_dir (1, 1) string = ""
      kwargs.case_ids (1, :) string = strings(1, 0)
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources)} = ...
         icemodel.verification.namelists.rcmsources()
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.points double = []
      kwargs.anchors = []
      kwargs.years (1, :) double = 2012:2018
      kwargs.radius_km (1, 1) double {mustBePositive} = 7.5
      kwargs.colocation_threshold_km (1, 1) double {mustBePositive} = 7.5
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

   forcing_sources = ...
      icemodel.verification.setup.normalizeForcingSources( ...
      kwargs.forcing_sources, kwargs.build_forcing);
   kwargs.forcing_sources = forcing_sources;

   % Validate the optional clamp before any cache or staging side effect.
   [window_start, window_end, window_enabled] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);

   % Resolve the family identity and requested runtime source sets once.
   dataset_family = "sumup";
   rcm_sources = intersect(forcing_sources, ...
      icemodel.verification.namelists.rcmsources(), "stable");
   build_rcm_forcing = kwargs.build_forcing && ~isempty(rcm_sources);

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
   if ~kwargs.dry_run && kwargs.build_observations
      icemodel.helpers.ensureDirExists(family_root);
   end

   % Resolve the optional observation bounds before point discovery so the
   % forcing-only fast path needs only explicit staged case ids.
   if window_enabled
      obs_start = window_start;
      obs_end = window_end;
   else
      % Blank reader bounds retain the all-available SUMup observation policy.
      obs_start = "";
      obs_end = "";
   end

   % Build either manifest-reuse state or freshly staged point state, then send
   % both through the one shared dry-run/import finalization block below.
   leg_callback = [];
   reuse_existing = ~kwargs.dry_run && ~kwargs.build_observations;
   if reuse_existing
      requested_ids = lower(erase(reshape(kwargs.case_ids, 1, []), "_"));
      if isempty(requested_ids)
         error( ...
            'icemodel:verification:importSumup:fastPathRequiresCaseIds', ...
            ['build_observations=false requires explicit case_ids; point ' ...
            'discovery belongs to observation staging'])
      end
      reuse_sources = strings(1, 0);
      coverage = struct();
      f_start = NaT;
      f_end = NaT;
      if build_rcm_forcing
         reuse_sources = rcm_sources;
         [f_start, f_end] = forcingWindow( ...
            window_start, window_end, kwargs.years);
         coverage = icemodel.verification.setup.promiceSourceCoverage( ...
            rcm_sources, struct('mar', kwargs.mar_dir, ...
            'merra', kwargs.merra_dir, 'racmo', kwargs.racmo_dir));
      end
      % Reuse the staged case entry so an RCM-only attachment does not require
      % source caches.
      [state, alive, skipped] = ...
         icemodel.verification.setup.reuseDatasetFamilyCases( ...
         manifest_file, requested_ids, emptyState(), ...
         dataset_family=dataset_family, ...
         overwrite_family=kwargs.overwrite_family, ...
         forcing_sources=reuse_sources, coverage=coverage, ...
         startdate=kwargs.startdate, enddate=kwargs.enddate, ...
         forcing_startdate=f_start, forcing_enddate=f_end);
      for k = 1:numel(state)
         % Reused observations keep their scientific case ids while RCM files
         % use the same safe storage identity as a full SUMup import.
         state(k).storage_alias = ...
            icemodel.verification.setup.rcmStorageAlias( ...
            dataset_family, state(k).case_id);
      end
      if build_rcm_forcing
         % A reused all-available case can predate every requested RCM. Reject
         % only a disjoint leg before delegation; partial overlaps retain the
         % shared maximum source window for efficient artifact reuse.
         leg_callback = @(s, src) skipNonoverlappingLeg( ...
            s.leg.(char(src)), s.period);
      end
   else
      cases = kwargs.points;
      catalog_points = isempty(cases) && isempty(kwargs.case_ids);
      canonical_default = catalog_points && isempty(kwargs.anchors);
      if isempty(cases) && ~isempty(kwargs.anchors)
         cases = pointsFromAnchors(kwargs.anchors);
      elseif isempty(cases)
         cases = defaultAnchorPoints(evaluation_data_root);
      end

      if isempty(cases)
         error('icemodel:verification:importSumup:noPoints', ...
            ['no SUMup points to stage and no mixed-anchor coordinates ' ...
            'available; provide points=[lat lon; ...] or anchors=...'])
      end

      % Validate caches only when building observations.
      % Dry runs remain metadata-only; optional skips stay quiet while required
      % SUMup products print their retrieval guidance before failing.
      if kwargs.dry_run
         source_dir = icemodel.verification.setup.sumupCacheDir(source_dir);
         source_status = [];
      else
         [source_dir, source_status] = ...
            icemodel.verification.setup.fetchSumup( ...
            cache_dir=icemodel.verification.setup.sumupCacheDir(source_dir), ...
            strict=~kwargs.skip_missing, silent=kwargs.skip_missing);
      end

      proj = icemodel.forcing.helpers.psnProjection();
      requested_ids = sumupRequestedIds(cases, proj, kwargs);
      if catalog_points
         [cases, requested_ids] = ...
            deduplicateCatalogPoints( ...
            cases, requested_ids, canonical_default);
      end
      stage_callback = @(~, n) stageSumupPoint( ...
         cases(n, :), n, source_dir, family_root, proj, source_status, ...
         obs_start, obs_end, window_start, window_end, kwargs);

      % Stage each requested case into importer state.
      [state, alive, skipped] = ...
         icemodel.verification.setup.stageDatasetFamilyCases( ...
         1:numel(requested_ids), emptyState(), stage_callback, ...
         skip_missing=kwargs.skip_missing, ...
         warning_id="icemodel:verification:importSumup:caseSkipped", ...
         label_callback=@(~, n) requested_ids(n));

      % The delegated RCM pass receives source-local windows only after native
      % observation staging has completed.
      if ~kwargs.dry_run && build_rcm_forcing
         [f_start, f_end] = ...
            forcingWindow(window_start, window_end, kwargs.years);
         coverage = icemodel.verification.setup.promiceSourceCoverage( ...
            rcm_sources, struct('mar', kwargs.mar_dir, ...
            'merra', kwargs.merra_dir, 'racmo', kwargs.racmo_dir));
         leg = icemodel.verification.setup.resolveLegWindows( ...
            rcm_sources, coverage, f_start, f_end);
         % Resolve the global source inventory once, but reject point periods
         % with zero overlap before any expensive RCM builder runs. Partial
         % overlaps retain the shared maximum source window.
         leg_callback = @(s, src) skipNonoverlappingLeg( ...
            leg.(char(src)), s.period);
      end
   end

   % Record source provenance once so dry-run and persisted manifests match.
   source_doi = "10.18739/A2M61BR5M";
   source_url = "https://nsidc.org/data/g02288";
   source_version = "sumup2025";
   retrieval_date = string(datetime('today'));

   % SUMup alone finalizes comparison-window overlap while building each entry.
   entry_callback = @sumupCaseEntry;

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
      leg_callback=leg_callback, ...
      met_outdir=met_outdir, userdata_outdir=userdata_outdir, ...
      mar_dir=kwargs.mar_dir, merra_dir=kwargs.merra_dir, ...
      racmo_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
      method="nearest", dt_out=kwargs.dt_out, ...
      after_source_callback=@(st, live, src) ...
      icemodel.verification.setup.stageMarDensityProfiles( ...
      st, live, src, family_root=family_root, ...
      userdata_outdir=userdata_outdir, mar_dir=kwargs.mar_dir, ...
      overwrite=kwargs.overwrite));
end

function s = stageSumupPoint(point, n, source_dir, family_root, proj, ...
      source_status, obs_start, obs_end, window_start, window_end, kwargs)
   %STAGESUMUPPOINT Stage one SUMup point and return importer state.
   [x3413, y3413] = projfwd(proj, point(1), point(2));

   % Co-location records provenance only; it never rebuilds anchor files.
   [is_coloc, anchor, dist_km] = ...
      icemodel.verification.helpers.sumupColocation(x3413, y3413, ...
      threshold_km=kwargs.colocation_threshold_km, ...
      anchors=kwargs.anchors, output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);

   % Case ids and aliases stay compact because the family directory is already
   % encoded by the parent `sumup` folder.
   [case_id, case_dir, site_id, site_name] = ...
      resolveCaseId(kwargs.case_ids, n, is_coloc, anchor);

   site_location = struct( ...
      'lat_wgs84', point(1), 'lon_wgs84', point(2), ...
      'x_epsg3413', x3413, 'y_epsg3413', y3413, 'elev_m', NaN);

   obs_file = "observations.mat";
   if ~kwargs.dry_run
      % Stop a quiet, importer-level cache miss before the standalone builder's
      % strict fetch can print its user-facing retrieval banner.
      assertSumupSourceAvailable(source_status);
      % Observations are the primary target and are never blocked by forcing.
      [observations, obs_meta] = ...
         icemodel.verification.setup.buildSumupObservations(point, ...
         source_dir=source_dir, radius_km=kwargs.radius_km, ...
         startdate=obs_start, enddate=obs_end);
      if ~icemodel.verification.setup.hasObservationRecords(observations)
         error('icemodel:verification:importSumup:missingObservation', ...
            'no SUMup observations found for point %s', case_id)
      end
      period = observedPeriod(window_start, window_end, obs_start, ...
         obs_end, obs_meta);
      requested_case = struct('period', period, ...
         'site_location', site_location, 'artifact_metadata', obs_meta);
      case_root = fullfile(family_root, case_dir);
      write_observation = icemodel.verification.setup.prepareCaseRoot( ...
         case_root, kwargs.overwrite, obs_file, requested_case);
      if write_observation
         % Repeated non-overwrite imports keep the current observation bytes.
         targets = struct('format', 'subsurface_profile_bundle', ...
            'data', observations, 'metadata', obs_meta);
         targets = icemodel.verification.setup.stampArtifactMetadata(targets);
         save(fullfile(case_root, obs_file), 'targets');
      end
      comparison_variables = ...
         icemodel.verification.setup.sumupComparisonVariables(observations);
      observation_variables = ...
         sumupObservationVariables(observations, obs_meta);
   else
      comparison_variables = dryRunSumupComparisonVariables();
      observation_variables = ...
         dryRunSumupObservationVariables(kwargs.radius_km);
      period = struct('start', periodStr(window_start, obs_start), ...
         'end', periodStr(window_end, obs_end));
   end

   colocation = struct( ...
      'sumup', struct('kind', 'firn_profile_obs', ...
      'staged', ~kwargs.dry_run, ...
      'obs_file', char(fullfile(case_dir, obs_file)), ...
      'note', 'SUMup observation profiles (staged, not rebuilt).'));
   colocation.anchor = colocationRecord(is_coloc, anchor, dist_km, ...
      kwargs.colocation_threshold_km);

   [surface_zone, eval_target, permafrost_zone] = ...
      sumupZoneAndTarget(is_coloc, anchor);
   s = struct('case_id', string(case_id), 'storage_alias', ...
      icemodel.verification.setup.rcmStorageAlias("sumup", case_id), ...
      'site_id', string(site_id), ...
      'site_name', string(site_name), 'point', point, ...
      'site_location', site_location, ...
      'period', period, ...
      'evaluation_file_rel', char(fullfile(case_dir, obs_file)), ...
      'entry', struct(), 'colocation', colocation, 'leg', struct(), ...
      'comparison_variables', {comparison_variables}, ...
      'observation_variables', observation_variables, ...
      'surface_zone', string(surface_zone), ...
      'eval_target', {eval_target}, ...
      'permafrost_zone', string(permafrost_zone), ...
      'anchor_note', string(colocationNote(is_coloc, anchor)), ...
      'reuse_entry', false, 'dry_run', kwargs.dry_run);
end

function assertSumupSourceAvailable(status)
   %ASSERTSUMUPSOURCEAVAILABLE Throw a stable per-point missing-source error.
   if isempty(status) || all([status.present])
      return
   end
   error('icemodel:verification:importSumup:missingSource', ...
      'SUMup source cache is incomplete')
end

function requested_ids = sumupRequestedIds(points, proj, kwargs)
   %SUMUPREQUESTEDIDS Resolve merge-scope ids for every requested SUMup point.
   requested_ids = strings(1, size(points, 1));
   for n = 1:size(points, 1)
      point = points(n, :);
      try
         [x3413, y3413] = projfwd(proj, point(1), point(2));
         [is_coloc, anchor] = ...
            icemodel.verification.helpers.sumupColocation(x3413, y3413, ...
            threshold_km=kwargs.colocation_threshold_km, ...
            anchors=kwargs.anchors, output_root=kwargs.output_root, ...
            evaluation_data_root=kwargs.evaluation_data_root, ...
            icemodel_config_casename=kwargs.icemodel_config_casename);
         requested_ids(n) = resolveCaseId(kwargs.case_ids, n, is_coloc, ...
            anchor);
      catch
         requested_ids(n) = fallbackCaseId(kwargs.case_ids, n);
      end
   end
end

function case_id = fallbackCaseId(case_ids, n)
   %FALLBACKCASEID Return the explicit or numbered SUMup id for failure paths.
   if numel(case_ids) >= n && strlength(case_ids(n)) > 0
      case_id = lower(erase(case_ids(n), "_"));
   else
      case_id = sprintf("sumup%02d", n);
   end
end

function [t1, t2] = forcingWindow(window_start, window_end, years)
   %FORCINGWINDOW Resolve the co-located RCM forcing window for the SUMup legs.
   % The explicit comparison window wins; else the kwargs.years span; else
   % unbounded (NaT -> each source's full on-disk coverage in resolveLegWindows).
   if ~isnat(window_start) && ~isnat(window_end)
      t1 = window_start;
      t2 = window_end;
   elseif ~isempty(years)
      t1 = icemodel.verification.setup.ensureUtc( ...
         sprintf('%d-01-01', min(years)));
      t2 = icemodel.verification.setup.ensureUtc( ...
         sprintf('%d-12-31 23:00:00', max(years)));
   else
      t1 = NaT;
      t2 = NaT;
   end
end

%% Manifest assembly + RCM delegation
function s = emptyState()
   %EMPTYSTATE Prototype per-point staging state (preallocation seed).
   s = struct('case_id', "", 'storage_alias', "", ...
      'site_id', "", 'site_name', "", ...
      'point', [NaN NaN], 'site_location', struct(), ...
      'period', struct('start', '', 'end', ''), ...
      'evaluation_file_rel', '', ...
      'entry', struct(), 'colocation', struct(), 'leg', struct(), ...
      'comparison_variables', {strings(0, 1)}, ...
      'observation_variables', struct(), 'surface_zone', "", ...
      'eval_target', {strings(0, 1)}, 'permafrost_zone', "", ...
      'anchor_note', "", 'reuse_entry', false, ...
      'dry_run', false);
end

function entry = sumupCaseEntry(s)
   %SUMUPCASEENTRY Convert one SUMup state record to a manifest case entry.
   if s.reuse_entry
      % A forcing-only refresh keeps the staged observation contract verbatim
      % and updates only requested colocation legs plus derived source lists.
      % Apply the same overlap annotation as a normal observation build so a
      % fast refresh cannot advertise a different comparison contract.
      entry = s.entry;
      colocation = s.colocation;
      for source = reshape(string(fieldnames(s.leg)), 1, [])
         name = char(source);
         no_overlap = false;
         if isfield(colocation, name) && isfield(colocation.(name), 'reason')
            reason = lower(string(colocation.(name).reason));
            no_overlap = contains(reason, "no overlap") ...
               || contains(reason, "does not overlap");
         end
         % A repeat call must remain a no-op for a source that cannot overlap
         % this staged observation period. Restore any prior staged/diagnostic
         % leg verbatim, or omit the newly computed skip from a clean case.
         if isfield(colocation, name) ...
               && ~logical(colocation.(name).staged) && no_overlap
            if isfield(entry.colocation, name)
               colocation.(name) = entry.colocation.(name);
            else
               colocation = rmfield(colocation, name);
            end
         end
      end
      % Existing legs and newly requested forcing windows are already scoped per
      % leg. Do not relabel them from the observation period on a repeat call.
      entry.colocation = colocation;
      [forcing_sources, eval_sources] = ...
         icemodel.verification.setup.colocationSourceLists(entry.colocation);
      entry.forcing_sources = cellstr(forcing_sources);
      entry.eval_sources = cellstr(eval_sources);
      return
   end
   colocation = comparisonOverlapColocation(s.colocation, s.period);
   [forcing_sources, eval_sources] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);
   if s.dry_run && isempty(eval_sources)
      eval_sources = "sumup_obs";
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
      cellstr(forcing_sources)
      cellstr(eval_sources)
      cellstr(s.comparison_variables)
      s.observation_variables
      colocation
      'irregular'
      sprintf(['SUMup firn point%s; MAR/MERRA met + RACMO Data ' ...
      'co-located (forcing not bundled).'], char(s.anchor_note))};

   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);
end

%% Local helpers
function points = defaultAnchorPoints(evaluation_data_root)
   %DEFAULTANCHORPOINTS Mixed firn anchor [lat lon] from staged manifests.
   %
   % SUMup is a point collection with no curated default site list, so the
   % default points are the staged mixed-family anchors. Returns an empty 0x2
   % when no anchor manifests are present, so the caller can raise one explicit
   % no-points message rather than fabricating points.
   anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
      evaluation_data_root=evaluation_data_root, ...
      icemodel_config_casename="");
   points = pointsFromAnchors(anchors);
end

function points = pointsFromAnchors(anchors)
   %POINTSFROMANCHORS Resolve [lat lon] staging points from a mixed catalog.
   %
   % Only anchors with finite WGS84 coordinates become staging points. Anchors
   % without coordinates can still be used as metadata when explicit points are
   % supplied, but they cannot define point-selection geometry themselves.
   points = zeros(0, 2);
   if isempty(anchors) || ~isfield(anchors, 'lat_wgs84') ...
         || ~isfield(anchors, 'lon_wgs84')
      return
   end

   lat = [anchors.lat_wgs84];
   lon = [anchors.lon_wgs84];
   keep = isfinite(lat) & isfinite(lon);
   lat = lat(keep);
   lon = lon(keep);
   points = [lat.', lon.'];
   if isempty(points)
      return
   end

   % Multiple verification families can describe the same physical site. Default
   % SUMup staging needs one point per location, not one point per source alias.
   unique_rows = uniqueLocationRows(lat, lon);
   points = points(unique_rows, :);
end

function rows = uniqueLocationRows(lat, lon)
   %UNIQUELOCATIONROWS Keep the first anchor for each rounded WGS84 location.
   keys = strings(1, numel(lat));
   for k = 1:numel(lat)
      keys(k) = sprintf("%.3f:%.3f", lat(k), lon(k));
   end

   seen = strings(1, 0);
   rows = false(1, numel(keys));
   for k = 1:numel(keys)
      if ~ismember(keys(k), seen)
         rows(k) = true;
         seen(end + 1) = keys(k); %#ok<AGROW>
      end
   end
end

function [points, requested_ids] = deduplicateCatalogPoints( ...
      points, requested_ids, canonical_default)
   %DEDUPLICATECATALOGPOINTS Keep stable ids and canonicalize the 2025 default.
   %
   % Multiple staged source families can describe the same logical anchor with
   % slightly different coordinates. The default catalog path should stage one
   % SUMup case per resolved id; explicit points/case_ids still surface
   % duplicates through the manifest merge validator. MATLAB's stable indices
   % replace the former hand-maintained seen list without changing first-id wins.
   [~, first_rows] = unique(requested_ids, "stable");
   keep = false(1, numel(requested_ids));
   keep(first_rows) = true;
   points = points(keep, :);
   requested_ids = requested_ids(keep);

   if ~canonical_default
      return
   end

   % This fixed map is evidence for the immutable SUMup 2025 release, not a
   % dynamic proximity heuristic. A loser remains when its keeper is absent.
   losers = ["serb", "tasu", "thul", "thul2", "zaca", "zacl"];
   keepers = ["mit", "tasl", "thuu", "thuu", "zacu", "zacu"];
   drop = false(size(requested_ids));
   for k = 1:numel(losers)
      if any(requested_ids == keepers(k))
         drop = drop | requested_ids == losers(k);
      end
   end

   % The default mixed catalog already carries the identical externally owned
   % research_site/humphrey target, so it must not create a second SUMup case.
   drop = drop | requested_ids == "humphrey";

   % These staged anchors have no density, temperature, or SMB records within
   % the accepted 7.5 km radius in the immutable SUMup 2025 Greenland release.
   drop = drop | ismember(requested_ids, ["s23", "lynl", "lynt", "nukb"]);
   points = points(~drop, :);
   requested_ids = requested_ids(~drop);
end

function [case_id, case_dir, site_id, site_name] = ...
      resolveCaseId(case_ids, n, is_coloc, anchor)
   %RESOLVECASEID Resolve the n-th case id, dir, and site labels.
   %
   % The case_id and the on-disk case directory are the same compact id; the
   % family folder is already `sumup`, so NO `sumup` prefix is added. An
   % explicit case_ids(n) wins. Otherwise a KAN-co-located point inherits the
   % anchor's compact id (KAN_L -> kanl) and PROMICE-style site labels so SUMup
   % case_ids match the PROMICE convention; a non-co-located point falls back
   % to a numbered sumupNN id.
   if numel(case_ids) >= n && strlength(case_ids(n)) > 0
      case_id = lower(erase(case_ids(n), "_"));
      site_id = case_ids(n);
      site_name = case_ids(n);
   elseif is_coloc && ~isempty(anchor)
      site_id = string(anchor.site);                 % e.g. "KAN_L"
      case_id = lower(erase(site_id, "_"));           % e.g. "kanl"
      site_name = replace(site_id, "_", "-");         % e.g. "KAN-L"
   else
      case_id = sprintf("sumup%02d", n);
      site_id = case_id;
      site_name = case_id;
   end
   case_dir = case_id;
end

function rec = colocationRecord(is_coloc, anchor, dist_km, threshold_km)
   %COLOCATIONRECORD Build the JSON co-location provenance record.
   rec = struct('is_colocated', is_coloc, ...
      'threshold_km', threshold_km, ...
      'nearest_anchor', "", 'nearest_family', "", ...
      'nearest_source_id', "", 'distance_km', dist_km);
   if ~isempty(anchor)
      rec.nearest_anchor = string(anchor.site);
      if isfield(anchor, 'family')
         rec.nearest_family = string(anchor.family);
      end
      if isfield(anchor, 'source_id')
         rec.nearest_source_id = string(anchor.source_id);
      end
   end
end

function note = colocationNote(is_coloc, anchor)
   %COLOCATIONNOTE Short co-location phrase for the case note.
   if is_coloc && ~isempty(anchor)
      family = "anchor";
      if isfield(anchor, 'family') && strlength(string(anchor.family)) > 0
         family = string(anchor.family);
      end
      note = sprintf(' co-located with %s %s', family, anchor.site);
   else
      note = '';
   end
end

function [surface_zone, eval_target, permafrost_zone] = ...
      sumupZoneAndTarget(is_coloc, anchor)
   %SUMUPZONEANDTARGET Resolve a SUMup case surface_zone/eval_target/permafrost.
   %
   % When a SUMup point co-locates with a staged mixed-family anchor, inherit the
   % anchor's already-curated regime fields. Older/promice-only anchor structs
   % fall back to promiceSiteCatalog(anchor.site). Otherwise leave the fields empty
   % so the schema is satisfied without guessing a regime.
   surface_zone = "";
   eval_target = strings(0, 1);
   permafrost_zone = "";
   if is_coloc && ~isempty(anchor)
      if isfield(anchor, 'surface_zone')
         surface_zone = string(anchor.surface_zone);
      end
      if isfield(anchor, 'eval_target')
         eval_target = string(anchor.eval_target);
      end
      if isfield(anchor, 'permafrost_zone')
         permafrost_zone = string(anchor.permafrost_zone);
      end
      if strlength(surface_zone) > 0 && ~isempty(eval_target) ...
            && strlength(permafrost_zone) > 0
         return
      end
      try
         info = icemodel.verification.setup.promiceSiteCatalog(anchor.site);
         surface_zone = info.surface_zone;
         eval_target = string(info.eval_target);
         permafrost_zone = info.permafrost_zone;
      catch
         surface_zone = "";
         eval_target = strings(0, 1);
         permafrost_zone = "";
      end
   end
end

function vars = dryRunSumupComparisonVariables()
   %DRYRUNSUMUPCOMPARISONVARIABLES SUMup manifest axes used without raw files.
   vars = ["density"; "subsurface_temperature"; "smb"];
end

function observation_variables = dryRunSumupObservationVariables(radius_km)
   %DRYRUNSUMUPOBSERVATIONVARIABLES Observation metadata used without raw files.
   observation_variables = struct( ...
      'present_groups', {{'density', 'subsurface_temperature', 'smb'}}, ...
      'observation_source', "SUMup 2025 release (NSIDC G02288)", ...
      'selection_radius_km', radius_km, ...
      'density', 'firn/snow density profile (depth, density, error)', ...
      'subsurface_temperature', 'SUMup subsurface temperature T(z,t)', ...
      'smb', 'surface mass balance records (signed; accumulation or ablation)');
end

function observation_variables = ...
      sumupObservationVariables(observations, obs_meta)
   %SUMUPOBSERVATIONVARIABLES Metadata for SUMup groups actually staged.
   present = ...
      icemodel.verification.setup.sumupComparisonVariables(observations);
   observation_variables = struct( ...
      'present_groups', {cellstr(present)}, ...
      'observation_source', "SUMup 2025 release (NSIDC G02288)");
   if isfield(obs_meta, 'selection_radius_km')
      observation_variables.selection_radius_km = ...
         obs_meta.selection_radius_km;
   end
   if any(present == "density")
      observation_variables.density = ...
         'firn/snow density profile (depth, density, error)';
   end
   if any(present == "subsurface_temperature")
      observation_variables.subsurface_temperature = ...
         'SUMup subsurface temperature T(z,t)';
   end
   if any(present == "smb")
      observation_variables.smb = ...
         'surface mass balance records (signed; accumulation or ablation)';
   end
end

function s = periodStr(window_bound, obs_bound)
   %PERIODSTR Manifest period string: "" when the window is unbounded.
   if strcmp(string(obs_bound), "") || isnat(window_bound)
      s = '';
   else
      s = icemodel.verification.setup.formatManifestTime(window_bound);
   end
end

function period = observedPeriod(window_start, window_end, obs_start, ...
      obs_end, obs_meta)
   %OBSERVEDPERIOD Use explicit windows first, else actual SUMup coverage.
   if ~strcmp(string(obs_start), "") && ~isnat(window_start)
      period = struct('start', periodStr(window_start, obs_start), ...
         'end', periodStr(window_end, obs_end));
      return
   end

   period = struct('start', '', 'end', '');
   if ~isfield(obs_meta, 'observation_period_start') ...
         || ~isfield(obs_meta, 'observation_period_end')
      return
   end
   if isnat(obs_meta.observation_period_start) ...
         || isnat(obs_meta.observation_period_end)
      return
   end
   period.start = icemodel.verification.setup.formatManifestTime( ...
      obs_meta.observation_period_start);
   period.end = icemodel.verification.setup.formatManifestTime( ...
      obs_meta.observation_period_end);
end

function colocation = comparisonOverlapColocation(colocation, period)
   %COMPARISONOVERLAPCOLOCATION Clamp RCM legs to the comparable period.
   %
   % SUMup defaults to all available observations. RCM artifacts can still be
   % useful over the part of that observation period they actually overlap, so
   % manifests keep partial-overlap legs and record the usable clipped window.
   [case_start, case_end] = parsePeriod(period);
   if isnat(case_start) || isnat(case_end)
      return
   end

   for source = reshape(icemodel.verification.namelists.rcmsources(), 1, [])
      name = char(source);
      if ~isfield(colocation, name) || ~isstruct(colocation.(name))
         continue
      end
      leg = colocation.(name);
      if ~isfield(leg, 'staged') || ~logical(leg.staged) ...
            || ~isfield(leg, 'window')
         continue
      end
      [leg_start, leg_end] = parsePeriod(leg.window);
      if isnat(leg_start) || isnat(leg_end)
         continue
      end
      if leg_end < case_start || leg_start > case_end
         colocation.(name) = skippedComparisonLeg(leg, case_start, case_end);
      else
         colocation.(name).window = overlapWindow( ...
            leg_start, leg_end, case_start, case_end);
      end
   end
end

function leg = skipNonoverlappingLeg(leg, period)
   %SKIPNONOVERLAPPINGLEG Reject a zero-overlap RCM leg before staging.
   %
   % resolveLegWindows intentionally returns the source's all-available window
   % for an unbounded SUMup import. Each SUMup point has its own observation
   % period, however, so reject only zero-overlap points before an RCM builder
   % can write orphan files. Partially overlapping points retain the shared
   % maximum source window for efficient batching; the later colocation clamp
   % records each case's usable comparison window in its manifest.
   if ~logical(leg.staged)
      return
   end
   [case_start, case_end] = parsePeriod(period);
   if isnat(case_start) || isnat(case_end)
      return
   end
   if leg.end < case_start || leg.start > case_end
      leg.staged = false;
      leg.years = [];
      leg.start = NaT;
      leg.end = NaT;
      leg.reason = comparisonNoOverlapReason(case_start, case_end);
   end
end

function [t1, t2] = parsePeriod(period)
   %PARSEPERIOD Convert a manifest period/window struct to UTC datetimes.
   t1 = NaT;
   t2 = NaT;
   if ~isstruct(period) || ~isfield(period, 'start') || ~isfield(period, 'end')
      return
   end
   if strlength(string(period.start)) == 0 || strlength(string(period.end)) == 0
      return
   end
   t1 = icemodel.verification.setup.ensureUtc(string(period.start));
   t2 = icemodel.verification.setup.ensureUtc(string(period.end));
end

function window = overlapWindow(leg_start, leg_end, case_start, case_end)
   %OVERLAPWINDOW Return the bounded source/case overlap as manifest strings.
   window = struct( ...
      'start', icemodel.verification.setup.formatManifestTime( ...
      max(leg_start, case_start)), ...
      'end', icemodel.verification.setup.formatManifestTime( ...
      min(leg_end, case_end)));
end

function leg = skippedComparisonLeg(leg, case_start, case_end)
   %SKIPPEDCOMPARISONLEG Mark an RCM leg unusable for the SUMup case period.
   leg.staged = false;
   leg.reason = char(comparisonNoOverlapReason(case_start, case_end));
   for fieldname = ["met_files", "data_files", "model_output_files", ...
         "reference_file"]
      name = char(fieldname);
      if isfield(leg, name)
         leg = rmfield(leg, name);
      end
   end
end

function reason = comparisonNoOverlapReason(case_start, case_end)
   %COMPARISONNOOVERLAPREASON Describe an RCM/SUMup period mismatch once.
   reason = sprintf( ...
      'RCM window does not overlap SUMup observations %s to %s.', ...
      icemodel.verification.setup.formatManifestTime(case_start), ...
      icemodel.verification.setup.formatManifestTime(case_end));
end
