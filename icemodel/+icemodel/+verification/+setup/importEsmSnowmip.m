function manifest = importEsmSnowmip(source_dir, kwargs)
   %IMPORTESMSNOWMIP Stage ESM-SnowMIP site fixtures for the verification suite.
   %
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir)
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir, ...
   %     overwrite=true, case_ids=["cdp","sod"])
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir, ...
   %     case_ids="cdp", startdate="2005-10-01", enddate="2006-07-01")
   %
   %  Stages the requested ESM-SnowMIP site cases under
   %  demo/data/eval/snow/esm_snowmip/<sitename>/. All 10 reference sites are
   %  supported; the per-site forcing and observation artifacts are produced by
   %  the reusable builders buildEsmSnowmipForcing / buildEsmSnowmipObservations
   %  so the same conversion path is used for staging and for any future
   %  on-the-fly icemodel run.
   %
   %  Window resolution
   %    - With no startdate / enddate, each site uses its
   %      icemodel.verification.helpers.default_smoke_window (one snow
   %      water year). This is the canonical default-staged fixture.
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
   %    evaluation_data_root : string (default config-derived)
   %        Base evaluation-data root to stage into.
   %    icemodel_config_casename : string (default "test")
   %        Config casename when evaluation_data_root is blank.
   %    case_ids : string vector (default all 10 ESM-SnowMIP sites)
   %        Site case ids to stage.
   %    startdate : datetime / string ("" default)
   %        Optional explicit window start. When omitted, each site
   %        uses default_smoke_window.
   %    enddate : datetime / string ("" default)
   %        Optional explicit window end. Required when startdate is
   %        provided.
   %    overwrite : logical (default false)
   %        Refresh staged setup artifacts when true.
   %
   %  Returns
   %    manifest : struct  Family manifest also written to manifest.json.
   %
   %  Role
   %    Setup/update tooling. This function creates or refreshes staged data
   %    under demo/data/eval/snow and is not part of normal verification runs.
   %    Layout / file-presence guarantees come from
   %    icemodel.verification.setup.fetchEsmSnowmip; downstream from that
   %    point this importer can write the resolved output files directly
   %    without repeating per-file existence checks.
   %
   % See also: icemodel.verification.setup.fetchEsmSnowmip,
   %  icemodel.verification.setup.buildEsmSnowmipForcing,
   %  icemodel.verification.setup.buildEsmSnowmipObservations,
   %  icemodel.verification.namelists.snowmipsite,
   %  icemodel.verification.namelists.snowmipcatalog,
   %  icemodel.verification.helpers.default_smoke_window

   arguments
      source_dir (1, :) string
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.case_ids (1, :) string = ...
         icemodel.verification.namelists.snowmipsite()
      kwargs.startdate = ""
      kwargs.enddate   = ""
      kwargs.overwrite (1, 1) logical = false
   end

   % Validate the requested case ids against the canonical sitename namelist.
   case_ids = reshape(kwargs.case_ids, 1, []);
   bad = setdiff(case_ids, icemodel.verification.namelists.snowmipsite());
   if ~isempty(bad)
      error('unsupported ESM-SnowMIP case ids: %s', strjoin(bad, ', '));
   end

   % Explicit windows must be paired. Either both bounds are given (the
   % staged window is shared by every requested site), or both are blank
   % (each site falls back to its default_smoke_window). Mixing is an
   % error rather than a silent half-window.
   has_startdate = ~strcmp(string(kwargs.startdate), "");
   has_enddate   = ~strcmp(string(kwargs.enddate),   "");
   if has_startdate ~= has_enddate
      error(['icemodel:verification:importEsmSnowmip:halfWindow ' ...
         'startdate and enddate must be provided together']);
   end
   use_explicit_window = has_startdate && has_enddate;

   % Confirm the source-cache layout is complete before the per-site loop.
   % fetch is the single source of truth for "are the upstream files
   % present and readable?", so the importer does not repeat that check.
   icemodel.verification.setup.fetchEsmSnowmip( ...
      cache_dir=string(source_dir), sitenames=case_ids);

   % Name the source family and runnable cases once. dataset_family is the
   % staged source folder/manifest family; case_ids are the site cases inside
   % that family.
   dataset_family = "esm_snowmip";

   % Resolve the path to the dataset family sub-folder
   %   <snow_data_root>/esm_snowmip
   snow_data_root = icemodel.verification.helpers.snowDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   family_root = fullfile(snow_data_root, dataset_family);
   icemodel.helpers.ensureDirExists(family_root);

   % Resolve the path to the dataset family manifest
   %   <snow_data_root>/esm_snowmip/manifest.json
   manifest_file = fullfile(family_root, "manifest.json");

   % Stage each requested case.
   case_entries = cell(numel(case_ids), 1);
   for i = 1:numel(case_ids)
      sitename = case_ids(i);
      info = icemodel.verification.namelists.snowmipcatalog(sitename);

      % Resolve the staged time window. Explicit kwargs override the
      % per-site default smoke window.
      if use_explicit_window
         window_start = icemodel.verification.setup.ensureUtc(kwargs.startdate);
         window_end   = icemodel.verification.setup.ensureUtc(kwargs.enddate);
      else
         [window_start, window_end] = ...
            icemodel.verification.helpers.default_smoke_window(sitename);
      end

      % Resolve the per-site case root and the three staged-artifact paths.
      case_root = fullfile(family_root, sitename);
      forcing_output_file    = fullfile(case_root, "forcing.mat");
      evaluation_output_file = fullfile(case_root, "evaluation.mat");
      reference_output_file  = fullfile(case_root, "reference.mat");

      % prepareCaseRoot owns the overwrite guard. From here down, this
      % importer can write the resolved output files directly without
      % repeating per-file existence checks.
      icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);

      % Build forcing and observations through the reusable builders. Both
      % builders share NetCDF / time-window / channel-selection logic via
      % the local readers in this setup namespace.
      [forcing_tt, forcing_meta] = ...
         icemodel.verification.setup.buildEsmSnowmipForcing( ...
         sitename, source_dir=source_dir, ...
         startdate=window_start, enddate=window_end);

      [obs_tt, obs_meta] = ...
         icemodel.verification.setup.buildEsmSnowmipObservations( ...
         sitename, source_dir=source_dir, ...
         startdate=window_start, enddate=window_end);

      % Pack the artifacts in the common verification schema. forcing
      % and observation timetables are wrapped in struct envelopes
      % carrying provenance metadata so downstream consumers do not
      % need to know which builder produced them.
      forcing = struct( ...
         'format',   'timeseries', ...
         'data',     forcing_tt, ...
         'metadata', forcing_meta);

      targets = struct( ...
         'format',   'timeseries', ...
         'data',     obs_tt, ...
         'metadata', obs_meta);

      % Until production snow physics exists, the smoke reference
      % duplicates the staged observations so comparecase / plotcase
      % can run on a fresh clone (synthetic-snow hook still applies
      % when run_snow_verification_suite('run_icemodel', true)).
      reference = struct( ...
         'format',   'timeseries', ...
         'data',     obs_tt, ...
         'metadata', icemodel.verification.setup.metadataStruct({ ...
         'reference_kind', 'smoke_observed_reference'
         'notes', ['Observed series duplicated as the initial smoke ' ...
         'reference so the staged harness can run before a ' ...
         'snow model exists.']}));

      save(forcing_output_file, 'forcing');
      save(evaluation_output_file, 'targets');
      save(reference_output_file, 'reference');

      % Choose comparison variables from what's actually present in the
      % staged obs timetable. The obs builder decides which canonical
      % columns to stage based on upstream channel availability, so a
      % simple presence check is sufficient here.
      comparison_variables = obsComparisonVariables(obs_tt);

      % Provenance for the staged observations. The richer obs_meta
      % carries snow_depth / SWE source channels and soil-temperature
      % depths; these are reduced to a metadataStruct so the manifest
      % entry stays JSON-friendly.
      observation_variables = ...
         icemodel.verification.setup.metadataStruct({ ...
         'snow_depth_source', obs_meta.snow_depth_source
         'swe_source',        obs_meta.swe_source
         'soil_depths_m',     obs_meta.soil_depths_m});

      case_values = { ...
         char(sitename)
         'esm_site'
         char(sitename)
         char(info.long_name)
         char(fullfile(sitename, "forcing.mat"))
         char(fullfile(sitename, "evaluation.mat"))
         char(fullfile(sitename, "reference.mat"))
         'hourly'
         struct('start', char(string(window_start)), ...
         'end',   char(string(window_end)))
         cellstr(comparison_variables)
         observation_variables
         char(siteCaseNote(info))};

      case_entries{i} = icemodel.verification.setup.makeCaseManifestEntry( ...
         case_values);
   end

   % Provenance.
   source_doi = "10.1594/PANGAEA.897575";
   source_url = "https://doi.org/10.1594/PANGAEA.897575";
   source_version = "ESM-SnowMIP_all.zip";
   retrieval_date = string(datetime('today'));

   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      dataset_family, source_doi, source_url, source_version, retrieval_date, ...
      vertcat(case_entries{:}));
   icemodel.verification.setup.writeManifest(manifest_file, manifest);
end

% =====================================================================
% Local helpers
% =====================================================================

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

function note = siteCaseNote(info)
   %SITECASENOTE Short manifest note describing the staged case.
   note = sprintf( ...
      'ESM-SnowMIP %s reference site (%s, %s).', ...
      info.long_name, info.location, info.sitename);
end
