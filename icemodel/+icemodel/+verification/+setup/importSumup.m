function manifest = importSumup(source_dir, kwargs)
   %IMPORTSUMUP Stage co-located SUMup firn evaluation cases.
   %
   %  manifest = icemodel.verification.setup.importSumup(source_dir)
   %  manifest = icemodel.verification.setup.importSumup(source_dir, ...
   %     points=[lat lon; ...], overwrite=true)
   %
   %  Stages SUMup firn observation cases under
   %  demo/data/eval/firn/sumup/<case_id>/, mirroring importPromiceSites'
   %  structure. For each selected SUMup point it:
   %    - reads the nearest SUMup firn observations (density / subsurface
   %      temperature / accumulation) via buildSumupObservations;
   %    - builds the co-located MAR point met + RACMO point Data via
   %      buildSumupForcing (SUMup points carry no station met);
   %    - records whether the point is within a co-location threshold of a
   %      PROMICE anchor (helpers.sumupColocation, default 7.5 km EPSG:3413);
   %    - writes evaluation.mat (SUMup observations as targets) + reference.mat
   %      (co-located RACMO Data) + a family manifest.json with
   %      case_type="firn_observational".
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
   %                 default data/verification/firn/sumup is used.
   %
   %  Name-value
   %    points : Nx2 double  [lat lon] WGS84 points to stage. Default: the
   %             curated PROMICE anchor coordinates from the firn/promice
   %             manifest (so SUMup cases co-locate with the anchors).
   %    case_ids : string vector  case ids for each point (default sumup_NN).
   %    years : numeric vector  forcing years (default 2012:2015, the RACMO-
   %            bound firn window shared with importPromiceSites).
   %    startdate / enddate : datetime or "" comparison-window clamp.
   %    radius_km : double (default 7.5)  SUMup point-selection radius.
   %    colocation_threshold_km : double (default 7.5)  PROMICE-anchor
   %            co-location threshold (see helpers.sumupColocation).
   %    mar_dir / racmo_dir / modis_dir : string source dirs for the builders.
   %    evaluation_data_root / input_data_root / icemodel_config_casename :
   %            staging roots (mirror importPromiceSites).
   %    overwrite : logical (default false)
   %    skip_missing : logical (default true)  record data-gated points in
   %            manifest.skipped rather than aborting.
   %
   %  Returns
   %    manifest : struct  family manifest also written to manifest.json.
   %
   %  Role
   %    Setup/update tooling. Creates or refreshes staged data under
   %    demo/data/eval/firn/sumup and is not part of normal verification runs.
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
   else
      window_start = icemodel.verification.setup.ensureUtc("2012-01-01");
      window_end = icemodel.verification.setup.ensureUtc("2015-12-31 23:00:00");
   end
   years = kwargs.years;

   % Resolve the firn/sumup family root + the input layout.
   dataset_family = "sumup";
   firn_data_root = icemodel.verification.helpers.firnDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   family_root = fullfile(firn_data_root, dataset_family);
   icemodel.helpers.ensureDirExists(family_root);
   manifest_file = fullfile(family_root, "manifest.json");

   input_root = icemodel.verification.helpers.inputDataRoot( ...
      "input_data_root", kwargs.input_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   met_outdir = fullfile(input_root, 'met');
   userdata_outdir = fullfile(input_root, 'userdata');

   proj = icemodel.forcing.helpers.psnProjection();

   case_entries = {};
   skipped = struct('site', {}, 'reason', {});

   for n = 1:n_points
      point = points(n, :);
      [case_id, alias] = resolveCaseId(kwargs.case_ids, n);

      try
         [x3413, y3413] = projfwd(proj, point(1), point(2));

         % Co-location with a PROMICE anchor (records provenance only).
         [is_coloc, anchor, dist_km] = ...
            icemodel.verification.helpers.sumupColocation(x3413, y3413, ...
            threshold_km=kwargs.colocation_threshold_km, ...
            evaluation_data_root=kwargs.evaluation_data_root, ...
            icemodel_config_casename=kwargs.icemodel_config_casename);

         site_location = struct( ...
            'lat_wgs84', point(1), 'lon_wgs84', point(2), ...
            'x_epsg3413', x3413, 'y_epsg3413', y3413, 'elev_m', NaN);

         case_root = fullfile(family_root, alias);
         icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);

         % --- SUMup observations (targets). ---
         [observations, obs_meta] = ...
            icemodel.verification.setup.buildSumupObservations(point, ...
            source_dir=source_dir, radius_km=kwargs.radius_km, ...
            startdate=window_start, enddate=window_end);

         % --- Co-located MAR met + RACMO Data. ---
         [forcing, ~] = ...
            icemodel.verification.setup.buildSumupForcing(point, years, ...
            mar_dir=kwargs.mar_dir, racmo_dir=kwargs.racmo_dir, ...
            modis_dir=kwargs.modis_dir, ...
            window_start=window_start, window_end=window_end);

         mar_files = icemodel.forcing.helpers.writemet( ...
            forcing.mar_met, alias, "mar", outdir=met_outdir, naming="window");
         racmo_files = icemodel.forcing.helpers.writeuserdata( ...
            forcing.racmo_data, alias, "racmo", outdir=userdata_outdir);

         colocated = struct( ...
            'mar', struct('kind', 'point_met', ...
               'met_files', relpaths(mar_files, met_outdir), ...
               'sample_method', 'nearest'), ...
            'racmo', struct('kind', 'point_data_smb_eval', ...
               'data_files', relpaths(racmo_files, userdata_outdir), ...
               'sample_method', 'nearest', ...
               'note', 'SMB/eval Data only; RACMO is not a met source.'), ...
            'colocation', colocationRecord(is_coloc, anchor, dist_km, ...
               kwargs.colocation_threshold_km));

         % --- Evaluation (SUMup obs) + reference (RACMO Data) artifacts. ---
         targets = struct('format', 'firn_profile_bundle', ...
            'data', observations, 'metadata', obs_meta);
         reference = struct('format', 'timeseries', ...
            'data', forcing.racmo_data, ...
            'metadata', icemodel.verification.setup.metadataStruct({ ...
            'reference_kind', 'colocated_rcm'
            'reference_source', 'RACMO2.3p3 FGRN11 (point extraction)'
            'note', 'SMB/eval Data only; RACMO carries no met channels.'}));

         save(fullfile(case_root, "evaluation.mat"), 'targets');
         save(fullfile(case_root, "reference.mat"), 'reference');

         comparison_vars = sumupComparisonVariables(observations);
         obs_vars = icemodel.verification.setup.metadataStruct({ ...
            'density', 'firn/snow density profile (depth, density, error)'
            'subsurface_temperature', 'SUMup subsurface temperature T(z,t)'
            'accumulation', 'SMB / accumulation records'});

         case_values = { ...
            char(alias)
            'firn_observational'
            char(case_id)
            char(case_id)
            site_location
            char(fullfile(alias, "evaluation.mat"))
            char(fullfile(alias, "reference.mat"))
            'irregular'
            struct('start', char(string(window_start)), ...
            'end', char(string(window_end)))
            cellstr(comparison_vars)
            obs_vars
            colocated
            sprintf(['SUMup firn point%s; MAR met + RACMO Data ' ...
            'co-located bundle.'], colocationNote(is_coloc, anchor))};

         case_entries{end+1} = ...
            icemodel.verification.setup.makeFirnCaseManifestEntry(case_values); %#ok<AGROW>

      catch err
         if ~kwargs.skip_missing
            rethrow(err)
         end
         skipped(end+1) = struct('site', case_id, ...
            'reason', string(err.message)); %#ok<AGROW>
         warning('icemodel:verification:importSumup:pointSkipped', ...
            'skipping %s: %s', case_id, err.message);
      end
   end

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

   icemodel.verification.setup.writeManifest(manifest_file, manifest);
end

%% Local helpers
function cache_dir = resolveCacheDir(source_dir)
   %RESOLVECACHEDIR Resolve the SUMup cache dir, defaulting to the standard one.
   if strlength(source_dir) > 0
      cache_dir = source_dir;
   else
      cache_dir = string(fullfile(icemodel.getpath('data'), ...
         'verification', 'firn', 'sumup'));
   end
end

function points = defaultAnchorPoints()
   %DEFAULTANCHORPOINTS PROMICE anchor [lat lon] from the firn/promice manifest.
   %
   % SUMup is a point collection with no curated default site list, so the
   % default points are the committed PROMICE anchors - this stages SUMup
   % cases that co-locate with the existing firn/promice bundles. Returns an
   % empty 0x2 when no promice manifest is present (the caller then errors with
   % an explicit no-points message rather than fabricating points).
   points = zeros(0, 2);
   firn_root = icemodel.verification.helpers.firnDataRoot( ...
      "icemodel_config_casename", "test");
   manifest_file = fullfile(firn_root, "promice", "manifest.json");
   if exist(manifest_file, 'file') ~= 2
      return
   end
   m = jsondecode(fileread(manifest_file));
   for n = 1:numel(m.cases)
      loc = m.cases(n).site_location;
      points(end + 1, :) = [loc.lat_wgs84, loc.lon_wgs84]; %#ok<AGROW>
   end
end

function [case_id, alias] = resolveCaseId(case_ids, n)
   %RESOLVECASEID Resolve the n-th case id, defaulting to sumup_NN.
   if numel(case_ids) >= n && strlength(case_ids(n)) > 0
      case_id = case_ids(n);
   else
      case_id = sprintf("sumup_%02d", n);
   end
   alias = lower(erase(case_id, "_"));
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

function note = colocationNote(is_coloc, anchor)
   %COLOCATIONNOTE Short co-location phrase for the case note.
   if is_coloc && ~isempty(anchor)
      note = sprintf(' co-located with PROMICE %s', anchor.site);
   else
      note = '';
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

function rel = relpaths(filenames, base)
   %RELPATHS Reduce absolute staged paths to base-relative names for JSON.
   filenames = string(filenames);
   base = string(base);
   rel = erase(filenames, base + filesep);
   rel = reshape(rel, 1, []);
end
