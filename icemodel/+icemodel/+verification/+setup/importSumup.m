function manifest = importSumup(source_dir, kwargs)
   %IMPORTSUMUP Stage co-located SUMup firn evaluation cases.
   %
   %  manifest = icemodel.verification.setup.importSumup(source_dir)
   %  manifest = icemodel.verification.setup.importSumup(source_dir, ...
   %     points=[lat lon; ...], overwrite=true)
   %  manifest = icemodel.verification.setup.importSumup(source_dir, ...
   %     anchors=anchor_catalog)
   %
   %  Stages SUMup firn observation cases under
   %  demo/data/eval/sumup/<case_id>/, mirroring importPromiceSites'
   %  structure. Observation import is DECOUPLED from RCM forcing: for each
   %  selected SUMup point it
   %    - reads the nearest SUMup firn observations (density / subsurface
   %      temperature / accumulation) via buildSumupObservations and stages them
   %      once as the bundled observations.mat eval target;
   %    - records whether the point is within a co-location threshold of a
   %      PROMICE anchor (helpers.sumupColocation, default 7.5 km EPSG:3413) -
   %      PROMICE forcing is NEVER rebuilt; the anchor's own staged met is reused
   %      via this provenance record;
   %    - writes a FORCING-AGNOSTIC family manifest.json
   %      (case_type="firn_observational");
   %  then, when build_forcing is true (the convenience bundle), delegates the
   %  co-located gridded forcing to icemodel.verification.setup.stageRcmForcing,
   %  which stages MAR + MERRA met+Data and RACMO Data at each point (SUMup points
   %  carry no station met). The forcing lives in separate, runtime-discoverable
   %  met/userdata files recorded in the manifest colocation record by source id,
   %  never bundled. Setting build_forcing=false imports observations only; the
   %  RCM forcing can be built later, independently, on the staged manifest.
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
   %             explicit anchors when anchors is provided, otherwise curated
   %             PROMICE anchor coordinates from the firn/promice manifest.
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
   %            period records "" / "" (mirrors the PROMICE all-available
   %            default). There is no buried 2012-2015 comparison window.
   %    radius_km : double (default 7.5)  SUMup point-selection radius.
   %    colocation_threshold_km : double (default 7.5)  PROMICE-anchor
   %            co-location threshold (see helpers.sumupColocation).
   %    mar_dir / merra_dir / racmo_dir / modis_dir : string source dirs for the
   %            co-located RCM builders (delegated to stageRcmForcing).
   %    build_forcing : logical (default false). When true the co-located
   %            MAR/MERRA met+Data and RACMO Data are staged after the SUMup
   %            observation import; when false ONLY the observations + manifest
   %            are written (forcing can be built later via stageRcmForcing).
   %    evaluation_data_root / input_data_root / icemodel_config_casename :
   %            staging roots (mirror importPromiceSites). DEFAULT (roots unset):
   %            resolves via icemodel_config_casename ("test"), which points at the
   %            COMMITTED <repo>/demo/data tree (icemodel.config) - NOT a research
   %            root. To stage the gitignored research set pass
   %            evaluation_data_root/input_data_root at <repo>/data explicitly.
   %    overwrite : logical (default false)  refresh a requested point's own
   %            staged case folder; other cases are never touched.
   %    overwrite_family : logical (default false)  force a FULL rewrite of the
   %            family manifest from the requested points alone. DEFAULT MERGE.
   %    skip_missing : logical (default true)  record data-gated points in
   %            manifest.skipped rather than aborting.
   %
   %  Incremental staging (MERGE by default)
   %    Staging a point ADDS or UPDATES only that case in the sumup family
   %    manifest and PRESERVES every other committed case + files byte for byte
   %    (icemodel.verification.setup.writeFamilyManifestMerge, shared with
   %    importPromiceSites). overwrite_family=true rebuilds the family root.
   %
   %  Returns
   %    manifest : struct  family manifest also written to manifest.json.
   %
   %  Role
   %    Setup/update tooling. Creates or refreshes staged data under
   %    demo/data/eval/sumup and is not part of normal verification runs.
   %
   %  Colocation architecture (LOCKED: forcing-agnostic eval): the SUMup
   %    observation parser (buildSumupObservations) reads the real 2025 Greenland
   %    files end-to-end and stages them as the bundled data-only observations.mat
   %    eval target. The co-located RCM forcing legs (MAR + MERRA met+Data, RACMO
   %    Data) are the forcing/reference side - staged through stageRcmForcing as
   %    INDIVIDUAL files recorded in the manifest colocation record by source id,
   %    never bundled. This mirrors importPromiceSites: the eval target is
   %    bundled, the forcing is not. The SUMup observations are the PRIMARY target
   %    and are staged first; a failing RCM source degrades to a skipped forcing
   %    leg, never a skipped SUMup point.
   %
   % See also: icemodel.verification.setup.fetchSumup,
   %  icemodel.verification.setup.buildSumupObservations,
   %  icemodel.verification.setup.stageRcmForcing,
   %  icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.helpers.sumupColocation

   arguments
      source_dir (1, 1) string = ""
      kwargs.points double = defaultAnchorPoints()
      kwargs.anchors = []
      kwargs.case_ids (1, :) string = strings(1, 0)
      kwargs.years (1, :) double = 2012:2018
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.radius_km (1, 1) double {mustBePositive} = 7.5
      kwargs.colocation_threshold_km (1, 1) double {mustBePositive} = 7.5
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = true
      kwargs.build_forcing (1, 1) logical = false
   end

   points = kwargs.points;
   if isempty(points) && ~isempty(kwargs.anchors)
      points = pointsFromAnchors(kwargs.anchors);
   end

   if isempty(points)
      error('icemodel:verification:importSumup:noPoints', ...
         ['no SUMup points to stage and no mixed-anchor coordinates ' ...
         'available; provide points=[lat lon; ...] or anchors=...'])
   end
   n_points = size(points, 1);

   % Verify the SUMup cache up front (single source of truth for presence).
   % strict=true errors with the retrieval banner when the cache is empty, so
   % nothing is staged from a missing cache.
   source_dir = icemodel.verification.setup.fetchSumup( ...
      cache_dir=resolveCacheDir(source_dir), strict=~kwargs.skip_missing);

   % Resolve the comparison window (paired or both-blank).
   has_start = ~strcmp(string(kwargs.startdate), "");
   has_end = ~strcmp(string(kwargs.enddate), "");
   if has_start ~= has_end
      error('icemodel:verification:importSumup:halfWindow', ...
         'startdate and enddate must be provided together')
   end
   if has_start && has_end
      window_start = icemodel.verification.setup.ensureUtc(kwargs.startdate);
      window_end = icemodel.verification.setup.ensureUtc(kwargs.enddate);
      obs_start = window_start;
      obs_end = window_end;
   else
      % DEFAULT (omitted) is ALL AVAILABLE: no comparison clamp, so each SUMup
      % point stages its full on-disk observation record (mirrors the PROMICE
      % all-available-years default). The obs/forcing builders treat ""
      % (and NaT) as unbounded; the manifest period records "" / "" too.
      window_start = NaT;
      window_end = NaT;
      obs_start = "";
      obs_end = "";
   end

   % Resolve the sumup family root + the input layout.
   dataset_family = "sumup";
   evaluation_data_root = icemodel.verification.helpers.evaluationDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   family_root = fullfile(evaluation_data_root, dataset_family);
   icemodel.helpers.ensureDirExists(family_root);
   manifest_file = fullfile(family_root, "manifest.json");

   input_root = icemodel.verification.helpers.inputDataRoot( ...
      "input_data_root", kwargs.input_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   met_outdir = fullfile(input_root, 'met');
   userdata_outdir = fullfile(input_root, 'userdata');

   proj = icemodel.forcing.helpers.psnProjection();

   % Per-point state; a point yields at most one staged case OR one skip and
   % exactly one requested id.
   state = repmat(emptyState(), 1, n_points);
   alive = false(1, n_points);
   skipped = repmat(struct('site', "", 'reason', ""), 1, n_points);
   n_skipped = 0;
   requested_ids = strings(1, n_points);
   n_requested = 0;

   % --- Pass 1: stage the SUMup observations per point (the PRIMARY target). ---
   % Observation import is DECOUPLED from RCM forcing: the obs are staged and the
   % manifest written BEFORE any forcing build (Pass 2), so a killed/aborted
   % forcing run never destroys the imported observations.
   for n = 1:n_points
      point = points(n, :);
      case_id = sprintf("point%02d", n);   % fallback id for the catch path

      try
         [x3413, y3413] = projfwd(proj, point(1), point(2));

         % Co-location with a PROMICE anchor (records provenance only).
         [is_coloc, anchor, dist_km] = ...
            icemodel.verification.helpers.sumupColocation(x3413, y3413, ...
            threshold_km=kwargs.colocation_threshold_km, ...
            anchors=kwargs.anchors, ...
            evaluation_data_root=kwargs.evaluation_data_root, ...
            icemodel_config_casename=kwargs.icemodel_config_casename);

         % Resolve the case id / on-disk dir / site labels. case_id and the
         % case dir are the SAME compact id (no redundant family prefix - the
         % family folder is already `sumup`); a KAN-co-located case inherits
         % the anchor's compact id so the SUMup case_id matches the PROMICE
         % convention (kanl/kanm/kanu, site_id KAN_L/M/U).
         [case_id, alias, site_id, site_name] = ...
            resolveCaseId(kwargs.case_ids, n, is_coloc, anchor);
         n_requested = n_requested + 1;
         requested_ids(n_requested) = string(case_id);

         site_location = struct( ...
            'lat_wgs84', point(1), 'lon_wgs84', point(2), ...
            'x_epsg3413', x3413, 'y_epsg3413', y3413, 'elev_m', NaN);

         case_root = fullfile(family_root, alias);
         icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);

         % SUMup observations (targets), staged once as observations.mat (NOT
         % rebuilt later). The observations are the PRIMARY target and are never
         % blocked by a forcing failure (forcing is the SECONDARY Pass-2 step).
         [observations, obs_meta] = ...
            icemodel.verification.setup.buildSumupObservations(point, ...
            source_dir=source_dir, radius_km=kwargs.radius_km, ...
            startdate=obs_start, enddate=obs_end);
         targets = struct('format', 'subsurface_profile_bundle', ...
            'data', observations, 'metadata', obs_meta);
         obs_file = "observations.mat";
         save(fullfile(case_root, obs_file), 'targets');

         % Colocation = the SUMup obs leg + the PROMICE-anchor provenance. The
         % co-located RCM forcing legs (MAR/MERRA met + RACMO Data) are added in
         % Pass 2 by stageRcmForcing; PROMICE forcing is NEVER rebuilt here (the
         % anchor's own staged met is reused via this provenance record).
         colocation = struct( ...
            'sumup', struct('kind', 'firn_profile_obs', 'staged', true, ...
               'obs_file', char(fullfile(alias, obs_file)), ...
               'note', 'SUMup observation profiles (staged, not rebuilt).'));
         colocation.anchor = colocationRecord(is_coloc, anchor, dist_km, ...
            kwargs.colocation_threshold_km);

         comparison_vars = sumupComparisonVariables(observations);
         obs_vars = icemodel.verification.setup.metadataStruct({ ...
            'density', 'firn/snow density profile (depth, density, error)'
            'subsurface_temperature', 'SUMup subsurface temperature T(z,t)'
            'smb', 'surface mass balance records (signed; accumulation or ablation)'});
         [zone, target, pfz] = sumupZoneAndTarget(is_coloc, anchor);

         state(n) = struct('case_id', string(case_id), ...
            'alias', string(alias), 'site_id', string(site_id), ...
            'site_name', string(site_name), 'point', point, ...
            'site_location', site_location, ...
            'period', struct('start', periodStr(window_start, obs_start), ...
               'end', periodStr(window_end, obs_end)), ...
            'obs_file_rel', char(fullfile(alias, obs_file)), ...
            'colocation', colocation, ...
            'comparison_vars', {comparison_vars}, 'obs_vars', obs_vars, ...
            'zone', string(zone), 'target', {target}, 'pfz', string(pfz), ...
            'anchor_note', string(colocationNote(is_coloc, anchor)));
         alive(n) = true;

      catch err
         if ~kwargs.skip_missing
            rethrow(err)
         end
         n_skipped = n_skipped + 1;
         skipped(n_skipped) = struct('site', case_id, ...
            'reason', string(err.message));
         warning('icemodel:verification:importSumup:pointSkipped', ...
            'skipping %s: %s', case_id, err.message);
      end
   end

   requested_ids = requested_ids(1:n_requested);

   % --- Pass 2: write the obs manifest first, then delegate the
   % co-located RCM forcing/Data (MAR/MERRA met + RACMO Data) to stageRcmForcing.
   manifest = assembleAndWrite(state, alive, skipped, n_skipped, ...
      dataset_family, manifest_file, requested_ids, kwargs);

   if kwargs.build_forcing
      % Forcing window: the explicit comparison window when given, else the
      % kwargs.years span (preserves the years kwarg's meaning), else unbounded
      % (each source's full on-disk coverage).
      [f_start, f_end] = forcingWindow(window_start, window_end, kwargs.years);
      % Persist the manifest after each RCM source so a partial forcing run
      % keeps the completed sources' legs.
      persist = @(st) assembleAndWrite(st, alive, skipped, n_skipped, ...
         dataset_family, manifest_file, requested_ids, kwargs);
      state = stageColocatedForcing(state, alive, f_start, f_end, ...
         met_outdir, userdata_outdir, kwargs, persist);
      manifest = persist(state);
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
   s = struct('case_id', "", 'alias', "", 'site_id', "", 'site_name', "", ...
      'point', [NaN NaN], 'site_location', struct(), ...
      'period', struct('start', '', 'end', ''), 'obs_file_rel', '', ...
      'colocation', struct(), 'comparison_vars', {strings(0, 1)}, ...
      'obs_vars', struct(), 'zone', "", 'target', {strings(0, 1)}, ...
      'pfz', "", 'anchor_note', "");
end

function manifest = assembleAndWrite(state, alive, skipped, n_skipped, ...
      dataset_family, manifest_file, requested_ids, kwargs)
   %ASSEMBLEANDWRITE Build the forcing-agnostic SUMup manifest + MERGE-write it.
   %
   % Called twice: once with observations-only colocation (before
   % any forcing build) and once after the RCM legs are merged in.
   case_entries = cell(1, nnz(alive));
   n_cases = 0;
   for n = 1:numel(state)
      if ~alive(n)
         continue
      end
      s = state(n);
      [forcing_sources, eval_sources] = ...
         icemodel.verification.setup.colocationSourceLists(s.colocation);

      case_values = { ...
         char(s.case_id)
         'firn_observational'
         char(s.site_id)
         char(s.site_name)
         char(s.zone)
         cellstr(s.target)
         char(s.pfz)
         s.site_location
         s.period
         s.obs_file_rel
         cellstr(forcing_sources)
         cellstr(eval_sources)
         cellstr(s.comparison_vars)
         s.obs_vars
         s.colocation
         'irregular'
         sprintf(['SUMup firn point%s; MAR/MERRA met + RACMO Data ' ...
         'co-located (forcing not bundled).'], char(s.anchor_note))};

      n_cases = n_cases + 1;
      case_entries{n_cases} = ...
         icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);
   end
   case_entries = case_entries(1:n_cases);

   % Family manifest. Provenance points at the SUMup release; the co-located
   % MAR/MERRA/RACMO per-model provenance lives in each builder.
   source_doi = "10.18739/A2M61BR5M";
   source_url = "https://nsidc.org/data/g02288";
   source_version = "sumup2024[mar+merra+racmo colocated]";
   retrieval_date = string(datetime('today'));

   if isempty(case_entries)
      cases = struct([]);
   else
      cases = vertcat(case_entries{:});
   end
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      dataset_family, source_doi, source_url, source_version, ...
      retrieval_date, cases);
   manifest.skipped = skipped(1:n_skipped);

   % MERGE by default: add/update only the requested points' cases, preserving
   % every other committed SUMup case + files untouched (shared helper with
   % importPromiceSites). overwrite_family forces a full rewrite.
   manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
      manifest_file, manifest, requested_ids=requested_ids, ...
      overwrite_family=kwargs.overwrite_family);
end

function state = stageColocatedForcing(state, alive, window_start, window_end, ...
      met_outdir, userdata_outdir, kwargs, persist)
   %STAGECOLOCATEDFORCING Delegate the co-located MAR/MERRA/RACMO legs.
   %
   % SUMup points carry no station met, so PROMICE is NOT a forcing leg here -
   % only the gridded RCMs are staged. Each source is staged in its OWN
   % stageRcmForcing call, the legs merged into each point's state, and `persist`
   % MERGE-writes the manifest after each source with a per-source progress
   % line. stageRcmForcing writes MAR/MERRA met+Data and RACMO Data,
   % degrading a failing source's legs to skip-with-reason without losing the
   % others. The forcing window follows the comparison window (NaT = all-available
   % = each source's full on-disk coverage), resolved fail-early per source.
   rcm_models = ["mar", "merra", "racmo"];
   alive_idx = find(alive);
   if isempty(alive_idx)
      return
   end

   coverage = icemodel.verification.setup.promiceSourceCoverage(rcm_models, ...
      struct('mar', kwargs.mar_dir, 'merra', kwargs.merra_dir, ...
      'racmo', kwargs.racmo_dir));
   leg = icemodel.verification.setup.resolveLegWindows(rcm_models, coverage, ...
      window_start, window_end);

   points = vertcat(state(alive_idx).point);
   for src = rcm_models
      fprintf('[staging] %s forcing for %d SUMup point(s)...\n', ...
         upper(char(src)), numel(alive_idx));

      legspec = repmat(legProto(src), 1, numel(alive_idx));
      for j = 1:numel(alive_idx)
         legspec(j).alias = state(alive_idx(j)).alias;
         legspec(j).(char(src)) = leg.(char(src));
      end

      colocation = icemodel.verification.setup.stageRcmForcing(points, ...
         legspec=legspec, models=src, ...
         met_outdir=met_outdir, userdata_outdir=userdata_outdir, ...
         mar_dir=kwargs.mar_dir, merra_dir=kwargs.merra_dir, ...
         racmo_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
         method="nearest");

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

%% Local helpers
function cache_dir = resolveCacheDir(source_dir)
   %RESOLVECACHEDIR Resolve the SUMup cache dir, defaulting to the standard one.
   if strlength(source_dir) > 0
      cache_dir = source_dir;
   else
      cache_dir = string(fullfile(icemodel.getpath('data'), ...
         'verification', 'sumup'));
   end
end

function points = defaultAnchorPoints()
   %DEFAULTANCHORPOINTS PROMICE anchor [lat lon] from the promice manifest.
   %
   % SUMup is a point collection with no curated default site list, so the
   % default points are the committed PROMICE anchors - this stages SUMup
   % cases that co-locate with the existing promice bundles. Returns an
   % empty 0x2 when no promice manifest is present (the caller then errors with
   % an explicit no-points message rather than fabricating points).
   points = zeros(0, 2);
   eval_root = icemodel.verification.helpers.evaluationDataRoot( ...
      "icemodel_config_casename", "test");
   manifest_file = fullfile(eval_root, "promice", "manifest.json");
   if ~isfile(manifest_file)
      return
   end
   m = jsondecode(fileread(manifest_file));
   % The case count is known from the decoded manifest, so size the output up
   % front and fill it (no growing in the loop).
   points = zeros(numel(m.cases), 2);
   for n = 1:numel(m.cases)
      loc = m.cases(n).site_location;
      points(n, :) = [loc.lat_wgs84, loc.lon_wgs84];
   end
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
   points = [lat(keep).', lon(keep).'];
end

function [case_id, alias, site_id, site_name] = ...
      resolveCaseId(case_ids, n, is_coloc, anchor)
   %RESOLVECASEID Resolve the n-th case id, dir, and site labels.
   %
   % The case_id and the on-disk case dir (alias) are the SAME compact id; the
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
   alias = case_id;
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

function [zone, target, pfz] = sumupZoneAndTarget(is_coloc, anchor)
   %SUMUPZONEANDTARGET Resolve a SUMup case surface_zone/eval_target/permafrost.
   %
   % When a SUMup point co-locates with a curated PROMICE anchor, the
   % glaciological zone, capability descriptor, and permafrost zone are inherited
   % from promicesiteinfo(anchor.site) (the single source of truth). Otherwise all
   % are left empty ("" / string(0,1)) so the schema is satisfied without guessing
   % a regime.
   zone = "";
   target = strings(0, 1);
   pfz = "";
   if is_coloc && ~isempty(anchor)
      try
         info = icemodel.verification.helpers.promicesiteinfo(anchor.site);
         zone = info.surface_zone;
         target = string(info.eval_target);
         pfz = info.permafrost_zone;
      catch
         zone = "";
         target = strings(0, 1);
         pfz = "";
      end
   end
end

function vars = sumupComparisonVariables(observations)
   %SUMUPCOMPARISONVARIABLES Comparison axes present in the SUMup obs bundle.
   candidate = ["density"; "subsurface_temperature"; "smb"];
   present = false(numel(candidate), 1);
   for k = 1:numel(candidate)
      present(k) = isfield(observations, candidate(k)) ...
         && ~isempty(observations.(candidate(k)));
   end
   vars = candidate(present);
end

function s = periodStr(window_bound, obs_bound)
   %PERIODSTR Manifest period string: "" when the window is unbounded.
   if strcmp(string(obs_bound), "") || isnat(window_bound)
      s = '';
   else
      s = char(string(window_bound));
   end
end
