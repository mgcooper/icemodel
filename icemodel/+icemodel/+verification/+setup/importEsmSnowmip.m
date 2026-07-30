function manifest = importEsmSnowmip(source_dir, kwargs)
   %IMPORTESMSNOWMIP Stage ESM-SnowMIP site fixtures for the verification suite.
   %
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir)
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir, ...
   %     overwrite=true, case_ids=["cdp","sod"])
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir, ...
   %     case_ids="cdp", startdate="2005-10-01", enddate="2006-07-01")
   %
   %  Role
   %    Setup/update tooling. Creates or refreshes staged data under the resolved
   %    data/eval/esm_snowmip tree and is not part of normal verification runs.
   %
   %  Default roots
   %    source_dir="" reads the configured ESM-SnowMIP cache. With no output_root,
   %    cases go to <repo>/data/eval/esm_snowmip/<case_id>/ and met files go to
   %    <repo>/data/input/met/esm_snowmip/.
   %
   %  Stages the requested ESM-SnowMIP site cases under the resolved
   %  data/eval/esm_snowmip/<sitename>/ tree. All 10 reference sites are
   %  supported; the per-site forcing and observation artifacts are produced by
   %  the reusable builders buildEsmSnowmipForcing / buildEsmSnowmipObservations
   %  so the same conversion path is used for staging and for any future
   %  on-the-fly icemodel run.
   %  This importer intentionally has no forcing_sources/build_observations split:
   %  each site's source forcing and observations share one requested window and
   %  are staged as one atomic ESM-SnowMIP case conversion. There is no optional
   %  RCM attachment path to refresh independently.
   %
   %  Window resolution
   %    - With no startdate / enddate, each site stages the full source record
   %      available in the ESM-SnowMIP forcing and observation files.
   %    - With explicit startdate and enddate, that single window is
   %      applied to every site listed in case_ids. Per-site staging
   %      is the natural unit of the importer; multi-window staging
   %      should be done by repeated calls.
   %
   %  Inputs
   %    source_dir : string
   %        Directory containing met_<kind>_<sitename>_*.nc and
   %        obs_insitu_<sitename>_*.nc files.
   %        See icemodel.verification.setup.fetchEsmSnowmip.
   %
   %  Name-value
   %    case_ids : string vector. Default is all 10 canonical ESM-SnowMIP sites.
   %        Pass an explicit list to stage a subset.
   %    startdate, enddate : datetime / string. Optional explicit forcing and
   %        observation window; pass both or neither. With both omitted, each
   %        site stages its full source record.
   %    output_root : string. Convenience root for eval and input artifacts.
   %        When set, eval artifacts go to <output_root>/eval and runtime files
   %        go to <output_root>/input.
   %    evaluation_data_root, input_data_root, icemodel_config_casename : lower-
   %        level root resolution when output_root is unset. With all roots
   %        unset, the importer targets the gitignored top-level <repo>/data tree.
   %    dt_out : model-met output timestep (default "15m"); pass "" to retain
   %        native model-met cadence.
   %    overwrite : logical (default false). Refresh requested case artifacts;
   %        other cases are never touched.
   %    overwrite_family : logical (default false). Replace the family manifest
   %        from the requested cases alone. The default is merge.
   %    skip_missing : logical (default false). Record source-local skips in
   %        manifest.skipped instead of aborting. The default import is strict,
   %        so absent source caches fail with retrieval instructions before any
   %        manifest merge can occur.
   %    dry_run : logical (default false). Return the manifest shape without
   %        reading source caches, writing artifacts, or merging the manifest.
   %
   %  Incremental staging (MERGE by default)
   %    Staging one case adds or updates only that case in the family manifest
   %    and preserves every other committed case and file. Re-staging the same
   %    case updates exactly its entry. Set overwrite_family=true only to
   %    deliberately rebuild the family root.
   %
   %  Returns
   %    manifest : struct  Family manifest also written to manifest.json.
   %
   %  Source guarantees
   %    Layout / file-presence guarantees come from
   %    icemodel.verification.setup.fetchEsmSnowmip; downstream from that
   %    point this importer can write the resolved output files directly
   %    without repeating per-file existence checks.
   %
   %  See also: icemodel.verification.setup.fetchEsmSnowmip,
   %    icemodel.verification.setup.buildEsmSnowmipForcing,
   %    icemodel.verification.setup.buildEsmSnowmipObservations,
   %    icemodel.verification.namelists.snowmipsite,
   %    icemodel.verification.setup.esmSnowmipSiteCatalog,
   %    icemodel.verification.helpers.default_smoke_window

   arguments
      source_dir (1, 1) string = ""
      kwargs.case_ids (1, :) string ...
         {icemodel.verification.validators.mustBeSnowmipSite} = ...
         icemodel.verification.namelists.snowmipsite()
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "15m"])} = "15m"
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = false
      kwargs.dry_run (1, 1) logical = false
   end

   cases = reshape(kwargs.case_ids, 1, []);
   requested_ids = lower(cases);

   % Validate the optional clamp before any cache or staging side effect.
   [window_start, window_end, window_enabled] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);

   % Validate source caches only for real staging. Dry runs remain metadata-only.
   if kwargs.dry_run
      source_dir = resolveDryRunSourceDir(source_dir);
   else
      source_dir = icemodel.verification.setup.fetchEsmSnowmip( ...
         cache_dir=string(source_dir), stations=cases, ...
         strict=~kwargs.skip_missing, silent=kwargs.skip_missing);
   end

   % Name the family once so dry-run and persisted paths use one identifier.
   dataset_family = "esm_snowmip";

   % Resolve output roots and paths before staging artifacts.
   % ESM-SnowMIP writes both evaluation cases and atomic input/met artifacts.
   [evaluation_data_root, input_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   [family_root, manifest_file] = ...
      icemodel.verification.setup.datasetFamilyStagingPaths( ...
      evaluation_data_root, input_root, dataset_family);
   if ~kwargs.dry_run
      icemodel.helpers.ensureDirExists(family_root);
   end

   stage_callback = @(~, n) stageCase( ...
      cases(n), source_dir, family_root, input_root, ...
      window_enabled, window_start, window_end, kwargs, dataset_family);

   % Stage each requested case into importer state.
   [state, alive, skipped] = ...
      icemodel.verification.setup.stageDatasetFamilyCases( ...
      1:numel(requested_ids), emptyState(), stage_callback, ...
      skip_missing=kwargs.skip_missing, ...
      warning_id="icemodel:verification:importEsmSnowmip:caseSkipped", ...
      label_callback=@(~, n) requested_ids(n));

   % Record source provenance once so dry-run and persisted manifests match.
   source_doi = "10.1594/PANGAEA.897575";
   source_url = "https://doi.org/10.1594/PANGAEA.897575";
   source_version = "ESM-SnowMIP_all.zip";
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

   % Persist each atomic case without separately attachable source legs.
   [manifest, ~] = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family=dataset_family, ...
      manifest_file=manifest_file, requested_ids=requested_ids, ...
      skipped=skipped, ...
      source_doi=source_doi, source_url=source_url, ...
      source_version=source_version, retrieval_date=retrieval_date, ...
      overwrite_family=kwargs.overwrite_family, overwrite=kwargs.overwrite, ...
      entry_callback=entry_callback);
end

%% Local helpers
function s = emptyState()
   %EMPTYSTATE Prototype dataset-family staging state.
   s = struct('case_id', "", 'entry', struct());
end

function s = stageCase(sitename, source_dir, family_root, input_root, ...
      window_enabled, window_start, window_end, kwargs, dataset_family)
   %STAGECASE Stage one requested case and return importer state.
   info = icemodel.verification.setup.esmSnowmipSiteCatalog(sitename);

   % Resolve the per-site case root and the staged obs-bundle path.
   case_root = fullfile(family_root, sitename);
   observations_output_file = fullfile(case_root, "observations.mat");

   if ~kwargs.dry_run
      % Build forcing and observations through the reusable builders. Explicit
      % kwargs request one shared window; otherwise the full source records set
      % the manifest period from their staged forcing/observation time axes.
      if window_enabled
         [forcing_tt, ~] = ...
            icemodel.verification.setup.buildEsmSnowmipForcing( ...
            sitename, source_dir=source_dir, startdate=window_start, ...
            enddate=window_end);
         [obs_tt, obs_meta] = ...
            icemodel.verification.setup.buildEsmSnowmipObservations( ...
            sitename, source_dir=source_dir, startdate=window_start, ...
            enddate=window_end);
      else
         [forcing_tt, ~] = ...
            icemodel.verification.setup.buildEsmSnowmipForcing( ...
            sitename, source_dir=source_dir);
         [obs_tt, obs_meta] = ...
            icemodel.verification.setup.buildEsmSnowmipObservations( ...
            sitename, source_dir=source_dir);
         [window_start, window_end] = stagedTimeBounds(forcing_tt, obs_tt);
      end

      % Create or clear the output root only after source reads succeeded, so
      % skip-missing imports do not delete existing cases or leave empty roots.
      requested_case = struct('period', struct( ...
         'start', icemodel.verification.setup.formatManifestTime(window_start), ...
         'end', icemodel.verification.setup.formatManifestTime(window_end)), ...
         'artifact_metadata', obs_meta);
      write_observation = icemodel.verification.setup.prepareCaseRoot( ...
         case_root, kwargs.overwrite, "observations.mat", requested_case);

      % Stage forcing under the standard icemodel input layout so loadmet consumes it
      % without verification-only branches.
      writeMetFiles(forcing_tt, sitename, input_root, kwargs.dt_out, ...
         kwargs.overwrite, dataset_family);

      if write_observation
         % Repeated non-overwrite imports keep the current observation bytes.
         targets = struct('format', 'timeseries', 'data', obs_tt, ...
            'metadata', obs_meta);
         targets = icemodel.verification.setup.stampArtifactMetadata(targets);
         save(observations_output_file, 'targets');
      end
      comparison_variables = obsComparisonVariables(obs_tt);
      observation_variables = icemodel.verification.setup.metadataStruct({ ...
         'snow_depth_source', obs_meta.snow_depth_source
         'swe_source', obs_meta.swe_source
         'soil_depths_m', obs_meta.soil_depths_m
         'obs_file', char(fullfile(sitename, "observations.mat"))});
   else
      [window_start, window_end] = ...
         icemodel.verification.helpers.default_smoke_window(sitename);
      comparison_variables = dryRunComparisonVariables();
      observation_variables = dryRunObservationVariables(sitename);
   end

   % Forcing-agnostic schema: evaluation_file references observations.mat and
   % reference_file is empty because the old smoke reference was redundant.
   % native_timestep records the staged model-met cadence used by standard
   % runtime filename resolution, not the raw hourly source cadence.
   staged_timestep = kwargs.dt_out;
   if staged_timestep == ""
      staged_timestep = "hourly";
   end
   case_values = { ...
      char(sitename)
      'esm_site'
      char(sitename)
      char(info.long_name)
      'land'
      {'seasonal_snow'}
      char(casePermafrostZone(sitename))
      char(fullfile(sitename, "observations.mat"))
      ''
      char(staged_timestep)
      struct('start', icemodel.verification.setup.formatManifestTime(window_start), ...
      'end', icemodel.verification.setup.formatManifestTime(window_end))
      cellstr(comparison_variables)
      observation_variables
      char(caseNote(info))};

   s = struct('case_id', sitename, ...
      'entry', icemodel.verification.setup.makeCaseManifestEntry(case_values));
end

%% ESM-SnowMIP source helpers
function [window_start, window_end] = stagedTimeBounds(forcing_tt, obs_tt)
   %STAGEDTIMEBOUNDS Return the full staged period from nonempty time axes.
   times = [forcing_tt.Time(:); obs_tt.Time(:)];
   times = times(~isnat(times));
   assert(~isempty(times), 'ESM-SnowMIP builders returned no timestamps')
   window_start = min(times);
   window_end = max(times);
end

function permafrost_zone = casePermafrostZone(sitename)
   %CASEPERMAFROSTZONE Obu et al. (2019) extent per SnowMIP site.
   %
   % Hard-coded results of a point-in-polygon test of the Obu et al. (2019) ESA
   % GlobPermafrost / UiO PEX permafrost-zone map at each ESM-SnowMIP site
   % (test/interactive/site_classification/classify_site_facies.m). All ten sites are off-ice land
   % surfaces. Sites outside any permafrost polygon -> "none". Vocabulary:
   % icemodel.verification.namelists.permafrostzone. Replaces the v1 Brown et al.
   % (1997) source.
   switch lower(string(sitename))
      case "sod"   % Sodankyla, boreal Lapland
         permafrost_zone = "continuous";
      case "snb"   % Senator Beck, Colorado alpine
         permafrost_zone = "discontinuous";
      case {"swa", "wfj"}  % Swamp Angel (CO), Weissfluhjoch (Alps)
         permafrost_zone = "sporadic";
      case "ojp"   % Old Jack Pine, boreal Saskatchewan
         permafrost_zone = "isolated";
      otherwise    % cdp, oas, obs, rme, sap: outside Obu permafrost polygons
         permafrost_zone = "none";
   end
end

function vars = obsComparisonVariables(obs_tt)
   %OBSCOMPARISONVARIABLES Pick comparison variables from staged obs.
   %
   % Returns a string column with the canonical ordering:
   %   snow_depth_m, swe_kg_m2, surface_temp_C, soil_temp_<k>_C
   %
   % A canonical variable is included when its column exists in the
   % staged obs timetable. The obs builder
   % (buildEsmSnowmipObservations) decides which canonical columns to
   % stage based on upstream NetCDF channel availability, so a simple
   % presence check is sufficient here. Sparseness within the staged
   % window is preserved on the comparison axis (plotcase renders
   % sparse markers); never-observed variables are absent and never
   % appear as comparison rows.
   present = string(obs_tt.Properties.VariableNames);
   canonical = ["snow_depth_m"; "swe_kg_m2"; "surface_temp_C"];
   soil = reshape(present(startsWith(present, "soil_temp_")), [], 1);
   vars = [canonical(ismember(canonical, present)); sort(soil)];
end

function vars = dryRunComparisonVariables()
   %DRYRUNCOMPARISONVARIABLES Manifest axes used without raw files.
   vars = ["snow_depth_m"; "swe_kg_m2"; "surface_temp_C"];
end

function observation_variables = dryRunObservationVariables(sitename)
   %DRYRUNOBSERVATIONVARIABLES Observation metadata used without raw files.
   observation_variables = icemodel.verification.setup.metadataStruct({ ...
      'snow_depth_source', 'ESM-SnowMIP snow-depth observation channels'
      'swe_source', 'ESM-SnowMIP SWE observation channels'
      'soil_depths_m', []
      'obs_file', char(fullfile(sitename, "observations.mat"))});
end

function source_dir = resolveDryRunSourceDir(source_dir)
   %RESOLVEDRYRUNSOURCEDIR Resolve the nominal source root without I/O.
   source_dir = icemodel.verification.setup.resolveFetchCacheDir( ...
      source_dir, fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'esm_snowmip'));
end

function note = caseNote(info)
   %CASENOTE Return the staged case description.
   note = sprintf( ...
      'ESM-SnowMIP %s reference site (%s, %s).', ...
      info.long_name, info.location, info.sitename);
end

function writeMetFiles(forcing_tt, sitename, input_root, dt_out, overwrite, ...
      dataset_family)
   %WRITEMETFILES Stage one multi-year met file under input/met/.
   %
   % Delegates met-file naming, validation, and saving to the shared
   % icemodel.forcing.helpers.writemet (window form, 15-minute by default), so
   % configureRun + loadmet resolve the file without verification-only branches.
   % Existing files are additive no-ops unless overwrite=true.

   icemodel.forcing.helpers.writemet(forcing_tt, sitename, dataset_family, ...
      outdir=fullfile(input_root, 'met'), naming="window", dt_out=dt_out, ...
      overwrite=overwrite);
end
