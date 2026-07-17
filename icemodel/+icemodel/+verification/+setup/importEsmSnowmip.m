function manifest = importEsmSnowmip(source_dir, kwargs)
   %IMPORTESMSNOWMIP Stage ESM-SnowMIP site fixtures for the verification suite.
   %
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir)
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir, ...
   %     overwrite=true, case_ids=["cdp","sod"])
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir, ...
   %     case_ids="cdp", startdate="2005-10-01", enddate="2006-07-01")
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
   %    output_root : string
   %        Base output root. When set, eval artifacts go to output_root/eval
   %        and met files go to output_root/input.
   %    evaluation_data_root : string (default config-derived)
   %        Base evaluation-data root to stage into.
   %    icemodel_config_casename : string (default "")
   %        Config casename when evaluation_data_root is blank.
   %    case_ids : string vector (default all 10 ESM-SnowMIP sites)
   %        Site case ids to stage.
   %    startdate : datetime / string ("" default)
   %        Optional explicit window start. When omitted, each site stages
   %        the full source record.
   %    enddate : datetime / string ("" default)
   %        Optional explicit window end. Required when startdate is
   %        provided.
   %    dt_out : string (default "15m")
   %        Model-met output timestep. Pass "" only to preserve the explicit
   %        native hourly cadence.
   %    overwrite : logical (default false)
   %        Refresh requested case folders when true.
   %    overwrite_family : logical (default false)
   %        Replace the family manifest instead of merging requested cases.
   %    skip_missing : logical (default false)
   %        When true, missing per-case source files are recorded as skipped cases.
   %        The default import path is strict so absent source caches fail with
   %        retrieval instructions before any manifest merge can occur.
   %    dry_run : logical (default false)
   %        Return a one-year smoke-window manifest preview without reading
   %        source files or writing eval/input artifacts.
   %
   %  Returns
   %    manifest : struct  Family manifest also written to manifest.json.
   %
   %  Role
   %    Setup/update tooling. This function creates or refreshes staged data
   %    under the resolved data/eval/esm_snowmip tree and is not part of normal verification
   %    runs.
   %    Layout / file-presence guarantees come from
   %    icemodel.verification.setup.fetchEsmSnowmip; downstream from that
   %    point this importer can write the resolved output files directly
   %    without repeating per-file existence checks.
   %
   % See also: icemodel.verification.setup.fetchEsmSnowmip,
   %  icemodel.verification.setup.buildEsmSnowmipForcing,
   %  icemodel.verification.setup.buildEsmSnowmipObservations,
   %  icemodel.verification.namelists.snowmipsite,
   %  icemodel.verification.setup.esmSnowmipSiteCatalog,
   %  icemodel.verification.helpers.default_smoke_window

   arguments
      source_dir (1, 1) string = ""
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.case_ids (1, :) string ...
         {icemodel.verification.validators.mustBeSnowmipSite} = ...
         icemodel.verification.namelists.snowmipsite()
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "15m"])} = "15m"
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = false
      kwargs.dry_run (1, 1) logical = false
   end

   case_ids = reshape(kwargs.case_ids, 1, []);

   % Normalize the optional shared site window before source/cache work.
   [window_start, window_end, window_enabled] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);

   % Confirm the source-cache layout before real staging. Dry runs are
   % metadata-only contract previews, so they resolve the default root without
   % creating cache folders or reading upstream NetCDF files.
   if kwargs.dry_run
      source_dir = resolveDryRunSourceDir(source_dir);
   else
      source_dir = icemodel.verification.setup.fetchEsmSnowmip( ...
         cache_dir=string(source_dir), sitenames=case_ids, ...
         strict=~kwargs.skip_missing, silent=kwargs.skip_missing);
   end

   % Name the source family and runnable cases once. dataset_family is the
   % staged source folder/manifest family; case_ids are the site cases inside
   % that family.
   dataset_family = "esm_snowmip";

   % Resolve the paired eval/input roots once, so staged met files land beside
   % the eval tree that declares them.
   [evaluation_data_root, input_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   family_root = fullfile(evaluation_data_root, dataset_family);
   if ~kwargs.dry_run
      icemodel.helpers.ensureDirExists(family_root);
   end

   % Resolve the path to the dataset family manifest
   %   <evaluation_data_root>/esm_snowmip/manifest.json
   manifest_file = fullfile(family_root, "manifest.json");

   % Stage each requested case into importer state.
   [state, alive, skipped] = ...
      icemodel.verification.setup.stageDatasetFamilyCases( ...
      1:numel(case_ids), emptyState(), ...
      @(~, n) stageEsmCase(case_ids(n), source_dir, family_root, ...
      input_root, window_enabled, window_start, window_end, kwargs), ...
      skip_missing=kwargs.skip_missing, ...
      warning_id="icemodel:verification:importEsmSnowmip:caseSkipped", ...
      label_callback=@(~, n) case_ids(n));

   % Record source provenance once so dry-run and persisted manifests match.
   source_doi = "10.1594/PANGAEA.897575";
   source_url = "https://doi.org/10.1594/PANGAEA.897575";
   source_version = "ESM-SnowMIP_all.zip";
   retrieval_date = string(datetime('today'));

   if kwargs.dry_run
      manifest = icemodel.verification.setup.runDatasetFamilyDryRun( ...
         state, alive, dataset_family=dataset_family, ...
         requested_ids=case_ids, skipped=skipped, ...
         source_doi=source_doi, source_url=source_url, ...
         source_version=source_version, retrieval_date=retrieval_date, ...
         entry_callback=@caseEntry);
      return
   end

   [manifest, ~] = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family=dataset_family, manifest_file=manifest_file, ...
      requested_ids=case_ids, skipped=skipped, source_doi=source_doi, ...
      source_url=source_url, source_version=source_version, retrieval_date=retrieval_date, ...
      overwrite_family=kwargs.overwrite_family, overwrite=kwargs.overwrite, ...
      entry_callback=@caseEntry);
end

function s = emptyState()
   %EMPTYSTATE Prototype ESM-SnowMIP staging state.
   s = struct('case_id', "", 'entry', struct());
end

function s = stageEsmCase(sitename, source_dir, family_root, input_root, ...
      window_enabled, window_start, window_end, kwargs)
   %STAGEESMCASE Stage one ESM-SnowMIP case and return importer state.
   info = icemodel.verification.setup.esmSnowmipSiteCatalog(sitename);

   % Resolve the staged time window. Explicit kwargs request one shared window;
   % otherwise the builders read the full source record and the manifest period
   % is taken from the staged forcing/observation time axes.
   % Resolve the per-site case root and the staged obs-bundle path.
   case_root = fullfile(family_root, sitename);
   observations_output_file = fullfile(case_root, "observations.mat");

   if ~kwargs.dry_run
      % Build forcing and observations through the reusable builders.
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
         kwargs.overwrite);

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
      comparison_variables = dryRunEsmComparisonVariables();
      observation_variables = dryRunEsmObservationVariables(sitename);
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
      char(snowmipPermafrostZone(sitename))
      char(fullfile(sitename, "observations.mat"))
      ''
      char(staged_timestep)
      struct('start', icemodel.verification.setup.formatManifestTime(window_start), ...
      'end', icemodel.verification.setup.formatManifestTime(window_end))
      cellstr(comparison_variables)
      observation_variables
      char(siteCaseNote(info))};

   s = struct('case_id', sitename, ...
      'entry', icemodel.verification.setup.makeCaseManifestEntry(case_values));
end

function entry = caseEntry(s)
   %CASEENTRY Return the staged ESM-SnowMIP case entry.
   entry = s.entry;
end

%% Local helpers
function [window_start, window_end] = stagedTimeBounds(forcing_tt, obs_tt)
   %STAGEDTIMEBOUNDS Return the full staged period from nonempty time axes.
   times = [forcing_tt.Time(:); obs_tt.Time(:)];
   times = times(~isnat(times));
   assert(~isempty(times), 'ESM-SnowMIP builders returned no timestamps')
   window_start = min(times);
   window_end = max(times);
end

function pfz = snowmipPermafrostZone(sitename)
   %SNOWMIPPERMAFROSTZONE Obu et al. (2019) permafrost extent per SnowMIP site.
   %
   % Hard-coded results of a point-in-polygon test of the Obu et al. (2019) ESA
   % GlobPermafrost / UiO PEX permafrost-zone map at each ESM-SnowMIP site
   % (test/interactive/site_classification/classify_site_facies.m). All ten sites are off-ice land
   % surfaces. Sites outside any permafrost polygon -> "none". Vocabulary:
   % icemodel.verification.namelists.permafrostzone. Replaces the v1 Brown et al.
   % (1997) source.
   switch lower(string(sitename))
      case "sod"   % Sodankyla, boreal Lapland
         pfz = "continuous";
      case "snb"   % Senator Beck, Colorado alpine
         pfz = "discontinuous";
      case {"swa", "wfj"}  % Swamp Angel (CO), Weissfluhjoch (Alps)
         pfz = "sporadic";
      case "ojp"   % Old Jack Pine, boreal Saskatchewan
         pfz = "isolated";
      otherwise    % cdp, oas, obs, rme, sap: outside Obu permafrost polygons
         pfz = "none";
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

function vars = dryRunEsmComparisonVariables()
   %DRYRUNESMCOMPARISONVARIABLES Source-light ESM manifest comparison axes.
   vars = ["snow_depth_m"; "swe_kg_m2"; "surface_temp_C"];
end

function obs_vars = dryRunEsmObservationVariables(sitename)
   %DRYRUNESMOBSERVATIONVARIABLES Metadata for source-light ESM dry runs.
   obs_vars = icemodel.verification.setup.metadataStruct({ ...
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

function note = siteCaseNote(info)
   %SITECASENOTE Short manifest note describing the staged case.
   note = sprintf( ...
      'ESM-SnowMIP %s reference site (%s, %s).', ...
      info.long_name, info.location, info.sitename);
end

function writeMetFiles(forcing_tt, sitename, input_root, dt_out, overwrite)
   %WRITEMETFILES Stage one multi-year met file under input/met/.
   %
   % Delegates met-file naming, validation, and saving to the shared
   % icemodel.forcing.helpers.writemet (window form, 15-minute by default), so
   % configureRun + loadmet resolve the file without verification-only branches.
   % Existing files are additive no-ops unless overwrite=true.

   icemodel.forcing.helpers.writemet(forcing_tt, sitename, "esm_snowmip", ...
      outdir=fullfile(input_root, 'met'), naming="window", dt_out=dt_out, ...
      overwrite=overwrite);
end
