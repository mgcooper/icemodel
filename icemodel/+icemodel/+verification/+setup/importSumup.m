function manifest = importSumup(source_dir, kwargs)
   %IMPORTSUMUP Stage co-located SUMup firn evaluation cases.
   %
   %  manifest = icemodel.verification.setup.importSumup(source_dir)
   %  manifest = icemodel.verification.setup.importSumup(source_dir, ...
   %     points=[lat lon; ...], overwrite=true)
   %
   %  Stages SUMup firn observation cases under
   %  demo/data/eval/sumup/<case_id>/, mirroring importPromiceSites'
   %  structure. For each selected SUMup point it:
   %    - reads the nearest SUMup firn observations (density / subsurface
   %      temperature / accumulation) via buildSumupObservations;
   %    - builds the co-located MAR point met + RACMO point Data via
   %      buildSumupForcing (SUMup points carry no station met);
   %    - records whether the point is within a co-location threshold of a
   %      PROMICE anchor (helpers.sumupColocation, default 7.5 km EPSG:3413);
   %    - writes a METADATA-ONLY family manifest.json
   %      (case_type="firn_observational"). The SUMup observation profiles are
   %      staged once as a per-case observations.mat obs bundle (not rebuilt);
   %      the co-located MAR met + RACMO Data are staged as individual files and
   %      recorded in the manifest colocation record by source id - no bundled
   %      evaluation.mat/reference.mat colocation copy is written.
   %
   %  Data-gated: SUMup is access-gated (NASA Earthdata, NSIDC G02288). With
   %  no populated cache, fetchSumup(strict=true) errors with the retrieval
   %  banner and nothing is staged. Provide points explicitly (SUMup is a
   %  point collection, not a curated site list); the default is the PROMICE
   %  anchor transect so the staged SUMup cases co-locate with the existing
   %  firn/promice bundles.
   %
   %  Inputs
   %    source_dir : string  SUMup cache dir (see fetchSumup). When blank, the
   %                 default data/verification/sumup is used.
   %
   %  Name-value
   %    points : Nx2 double  [lat lon] WGS84 points to stage. Default: the
   %             curated PROMICE anchor coordinates from the firn/promice
   %             manifest (so SUMup cases co-locate with the anchors).
   %    case_ids : string vector  case ids for each point (default sumup_NN).
   %    years : numeric vector  forcing years (default 2012:2015, the RACMO-
   %            bound firn window for the co-located MAR/RACMO legs).
   %    startdate / enddate : OPTIONAL comparison-window clamp; pass both or
   %            neither. The DEFAULT (omitted) is ALL AVAILABLE: no clamp, each
   %            point stages its full on-disk SUMup record and the manifest
   %            period records "" / "" (mirrors the PROMICE all-available
   %            default). There is no buried 2012-2015 comparison window.
   %    radius_km : double (default 7.5)  SUMup point-selection radius.
   %    colocation_threshold_km : double (default 7.5)  PROMICE-anchor
   %            co-location threshold (see helpers.sumupColocation).
   %    mar_dir / racmo_dir / modis_dir : string source dirs for the builders.
   %    evaluation_data_root / input_data_root / icemodel_config_casename :
   %            staging roots (mirror importPromiceSites). DEFAULT is SAFE: the
   %            roots resolve to the configured per-case RESEARCH root, NOT the
   %            committed demo/data/eval tree, unless a root pointing at
   %            demo/data is explicitly passed.
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
   %  Colocation architecture (LOCKED: metadata-only): the SUMup observation
   %    parser (buildSumupObservations) reads the real 2025 Greenland files
   %    end-to-end, and the co-located RCM forcing legs (buildSumupForcing: MAR
   %    met + RACMO Data) ARE staged through this driver - as INDIVIDUAL files
   %    recorded in the manifest colocation record by source id, NOT as a bundled
   %    per-case evaluation.mat/reference.mat. This mirrors importPromiceSites:
   %    colocation is metadata, not a bundled data copy. The SUMup observations
   %    are the PRIMARY target and are staged first; a throwing buildSumupForcing
   %    degrades to a skipped forcing leg, never a skipped SUMup point.
   %
   % See also: icemodel.verification.setup.fetchSumup,
   %  icemodel.verification.setup.buildSumupObservations,
   %  icemodel.verification.setup.buildSumupForcing,
   %  icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.helpers.sumupColocation

   arguments
      source_dir (1, 1) string = ""
      kwargs.points double = defaultAnchorPoints()
      kwargs.case_ids (1, :) string = strings(1, 0)
      kwargs.years (1, :) double = 2012:2015
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.radius_km (1, 1) double {mustBePositive} = 7.5
      kwargs.colocation_threshold_km (1, 1) double {mustBePositive} = 7.5
      kwargs.mar_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = true
   end

   if isempty(kwargs.points)
      error('icemodel:verification:importSumup:noPoints', ...
         ['no SUMup points to stage and no PROMICE anchor coordinates ' ...
         'available; provide points=[lat lon; ...]'])
   end
   points = kwargs.points;
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
      % RR1 all-available-years default). The obs/forcing builders treat ""
      % (and NaT) as unbounded; the manifest period records "" / "" too.
      window_start = NaT;
      window_end = NaT;
      obs_start = "";
      obs_end = "";
   end
   years = kwargs.years;

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

   % Every point yields at most one case OR one skip and exactly one requested
   % id, so preallocate to n_points and trim with the running counts (no growing
   % in the loop).
   case_entries = cell(1, n_points);
   n_cases = 0;
   skipped = repmat(struct('site', "", 'reason', ""), 1, n_points);
   n_skipped = 0;
   requested_ids = strings(1, n_points);
   n_requested = 0;

   for n = 1:n_points
      point = points(n, :);

      try
         [x3413, y3413] = projfwd(proj, point(1), point(2));

         % Co-location with a PROMICE anchor (records provenance only).
         [is_coloc, anchor, dist_km] = ...
            icemodel.verification.helpers.sumupColocation(x3413, y3413, ...
            threshold_km=kwargs.colocation_threshold_km, ...
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

         % --- SUMup observations (targets). ---
         [observations, obs_meta] = ...
            icemodel.verification.setup.buildSumupObservations(point, ...
            source_dir=source_dir, radius_km=kwargs.radius_km, ...
            startdate=obs_start, enddate=obs_end);

         % SUMup observation profiles are staged as a per-case obs bundle (NOT
         % rebuilt here; staged_kind=observations_only). The observations are
         % the PRIMARY target, so they are staged BEFORE the co-located forcing
         % legs and never blocked by a forcing failure (see the leg guard below).
         targets = struct('format', 'subsurface_profile_bundle', ...
            'data', observations, 'metadata', obs_meta);
         obs_file = "observations.mat";
         save(fullfile(case_root, obs_file), 'targets');

         colocation = struct( ...
            'sumup', struct('kind', 'firn_profile_obs', ...
               'obs_file', char(fullfile(alias, obs_file)), ...
               'note', 'SUMup observation profiles (staged, not rebuilt).'));

         % --- Co-located MAR met + RACMO Data (colocation METADATA). ---
         % The co-located RCM legs are SECONDARY: a throwing buildSumupForcing
         % (e.g. a source dir that vanished, or a read error) degrades to a
         % SKIPPED forcing leg, not a skipped SUMup point - the SUMup obs are the
         % primary target and have already been staged above. This mirrors the
         % per-leg guard the PROMICE importer uses (47f5a31).
         forcing_sources = strings(0, 1);
         eval_sources = "sumup_obs";
         try
            [forcing, ~] = ...
               icemodel.verification.setup.buildSumupForcing(point, years, ...
               mar_dir=kwargs.mar_dir, racmo_dir=kwargs.racmo_dir, ...
               modis_dir=kwargs.modis_dir, ...
               window_start=obs_start, window_end=obs_end);

            mar_files = icemodel.forcing.helpers.writemet( ...
               forcing.mar_met, alias, "mar", outdir=met_outdir, naming="window");
            racmo_files = icemodel.forcing.helpers.writeuserdata( ...
               forcing.racmo_data, alias, "racmo", outdir=userdata_outdir);

            colocation.mar = struct('kind', 'point_met', ...
               'met_files', icemodel.verification.setup.relpaths(mar_files, met_outdir), ...
               'sample_method', 'nearest');
            colocation.racmo = struct('kind', 'point_data_smb_eval', ...
               'data_files', icemodel.verification.setup.relpaths(racmo_files, userdata_outdir), ...
               'sample_method', 'nearest', ...
               'note', 'SMB/eval Data only; RACMO is not a met source.');
            forcing_sources = "mar";
            eval_sources = ["sumup_obs", "racmo"];
         catch leg_err
            % Record the co-located forcing as a skipped leg; obs still staged.
            colocation.mar = skippedLeg('point_met', leg_err.message);
            colocation.racmo = skippedLeg('point_data_smb_eval', leg_err.message);
            warning('icemodel:verification:importSumup:forcingSkipped', ...
               'forcing legs skipped for %s: %s', case_id, leg_err.message);
         end

         colocation.anchor = colocationRecord(is_coloc, anchor, dist_km, ...
            kwargs.colocation_threshold_km);

         comparison_vars = sumupComparisonVariables(observations);
         obs_vars = icemodel.verification.setup.metadataStruct({ ...
            'density', 'firn/snow density profile (depth, density, error)'
            'subsurface_temperature', 'SUMup subsurface temperature T(z,t)'
            'accumulation', 'SMB / accumulation records'});

         [zone, target, pfz] = sumupZoneAndTarget(is_coloc, anchor);

         case_values = { ...
            char(case_id)
            'firn_observational'
            char(site_id)
            char(site_name)
            char(zone)
            cellstr(target)
            char(pfz)
            site_location
            struct('start', periodStr(window_start, obs_start), ...
            'end', periodStr(window_end, obs_end))
            char(fullfile(alias, obs_file))
            cellstr(forcing_sources)
            cellstr(eval_sources)
            cellstr(comparison_vars)
            obs_vars
            colocation
            'irregular'
            sprintf(['SUMup firn point%s; MAR met + RACMO Data ' ...
            'co-located (metadata-only).'], colocationNote(is_coloc, anchor))};

         n_cases = n_cases + 1;
         case_entries{n_cases} = ...
            icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);

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

   % Trim the preallocated buffers to the counts actually produced.
   case_entries = case_entries(1:n_cases);
   skipped = skipped(1:n_skipped);
   requested_ids = requested_ids(1:n_requested);

   % Family manifest. Provenance points at the SUMup release; the co-located
   % MAR/RACMO per-model provenance lives in each builder.
   source_doi = "10.18739/A2M61BR5M";
   source_url = "https://nsidc.org/data/g02288";
   source_version = "sumup2024[mar+racmo colocated]";
   retrieval_date = string(datetime('today'));

   if isempty(case_entries)
      cases = struct([]);
   else
      cases = vertcat(case_entries{:});
   end
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      dataset_family, source_doi, source_url, source_version, ...
      retrieval_date, cases);
   manifest.skipped = skipped;

   % MERGE by default: add/update only the requested points' cases, preserving
   % every other committed SUMup case + files untouched (shared helper with
   % importPromiceSites). overwrite_family forces a full rewrite.
   manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
      manifest_file, manifest, requested_ids=requested_ids, ...
      overwrite_family=kwargs.overwrite_family);
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
      'nearest_anchor', "", 'distance_km', dist_km);
   if is_coloc && ~isempty(anchor)
      rec.nearest_anchor = string(anchor.site);
   end
end

function leg = skippedLeg(kind, reason)
   %SKIPPEDLEG Manifest entry for a co-located forcing leg that failed to build.
   %
   % Mirrors importPromiceSites' skipped-leg record so a throwing
   % buildSumupForcing degrades to a skipped leg (obs still staged), not a
   % skipped SUMup point.
   leg = struct('kind', kind, 'staged', false, ...
      'reason', char(string(reason)));
end

function note = colocationNote(is_coloc, anchor)
   %COLOCATIONNOTE Short co-location phrase for the case note.
   if is_coloc && ~isempty(anchor)
      note = sprintf(' co-located with PROMICE %s', anchor.site);
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
   candidate = ["density"; "subsurface_temperature"; "accumulation"];
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
