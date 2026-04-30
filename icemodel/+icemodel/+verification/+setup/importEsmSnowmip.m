function manifest = importEsmSnowmip(source_dir, kwargs)
   %IMPORTESMSNOWMIP Stage ESM-SnowMIP site fixtures for the verification suite.
   %
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir)
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir, ...
   %     overwrite=true, tier="full", case_ids=["cdp","sod"])
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir, ...
   %     tier="custom", startdate="2005-10-01", enddate="2006-07-01")
   %
   %  Stages the requested ESM-SnowMIP site cases under
   %  demo/data/eval/snow/esm_snowmip/<sitename>/. All 10 reference sites are
   %  supported; the per-site forcing and observation artifacts are produced by
   %  the reusable builders buildEsmSnowmipForcing / buildEsmSnowmipObservations
   %  so the same conversion path is used for staging and for any future
   %  on-the-fly icemodel run.
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
   %    tier : "smoke" | "full" | "custom" (default "smoke")
   %        Window selector. "smoke" stages a single representative
   %        snow year per site (cheap, agent feedback loop). "full"
   %        stages every full snow year in the site's coverage.
   %        "custom" requires startdate / enddate kwargs.
   %    startdate : datetime / string ("" default)
   %        Custom-tier window start. Only used when tier="custom".
   %    enddate : datetime / string ("" default)
   %        Custom-tier window end.
   %    overwrite : logical (default false)
   %        Refresh staged setup artifacts when true.
   %
   %  Returns
   %    manifest : struct  Family manifest also written to manifest.json.
   %
   % See also: icemodel.verification.setup.fetchEsmSnowmip,
   %  icemodel.verification.setup.buildEsmSnowmipForcing,
   %  icemodel.verification.setup.buildEsmSnowmipObservations,
   %  icemodel.verification.namelists.snowmipsite

   arguments
      source_dir (1, :) string
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.case_ids (1, :) string = ...
         icemodel.verification.namelists.caseid("esm_snowmip")
      kwargs.tier (1, 1) string {mustBeMember(kwargs.tier, ...
         ["smoke", "full", "custom"])} = "smoke"
      kwargs.startdate = ""
      kwargs.enddate   = ""
      kwargs.overwrite (1, 1) logical = false
   end

   % Validate the requested case ids against the canonical site catalog.
   valid_ids = icemodel.verification.namelists.caseid("esm_snowmip");
   case_ids = reshape(kwargs.case_ids, 1, []);
   bad = setdiff(case_ids, valid_ids);
   if ~isempty(bad)
      error('unsupported ESM-SnowMIP case ids: %s', strjoin(bad, ', '));
   end

   % Custom-tier requires explicit start / end dates. The constraint is
   % checked here so site iteration cannot accidentally produce
   % zero-row staged artifacts when one date is missing.
   if kwargs.tier == "custom"
      if strcmp(string(kwargs.startdate), "") || strcmp(string(kwargs.enddate), "")
         error('tier=custom requires explicit startdate and enddate kwargs');
      end
   end

   % Resolve and prepare the dataset-family root.
   dataset_family = "esm_snowmip";
   snow_data_root = icemodel.verification.helpers.snowDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   family_root = fullfile(snow_data_root, dataset_family);
   icemodel.helpers.ensureDirExists(family_root);
   manifest_file = fullfile(family_root, "manifest.json");

   % Stage each requested case.
   case_entries = cell(numel(case_ids), 1);
   for i = 1:numel(case_ids)
      sitename = case_ids(i);
      info = icemodel.verification.namelists.snowmipsite(sitename);

      % Resolve the requested time window per the selected tier.
      [window_start, window_end] = resolveTierWindow( ...
         kwargs.tier, info, kwargs.startdate, kwargs.enddate);

      case_root = fullfile(family_root, sitename);
      forcing_output_file    = fullfile(case_root, "forcing.mat");
      evaluation_output_file = fullfile(case_root, "evaluation.mat");
      reference_output_file  = fullfile(case_root, "reference.mat");

      icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);

      % Build forcing and observations through the reusable builders.
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

      % Choose comparison variables from what's actually present in
      % the obs timetable (some sites lack soil temperature / surface
      % temperature). This avoids comparecase emitting "not_applicable"
      % rows for variables that were never observed at the site.
      comparison_variables = obsComparisonVariables(obs_tt);

      observation_variables = ...
         icemodel.verification.setup.metadataStruct({ ...
         'snow_depth_source', obs_meta.snow_depth_source
         'swe_source',        obs_meta.swe_source
         'soil_depths_m',     obs_meta.soil_depths_m});

      case_values = { ...
         char(sitename)
         'esm_site'
         char(kwargs.tier)
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
         char(siteCaseNote(info, kwargs.tier))};

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

function [window_start, window_end] = resolveTierWindow( ...
      tier, info, custom_start, custom_end)
   %RESOLVETIERWINDOW Decide the staged time window per tier.

   % insitu_window in the namelist is [start_year, end_year] inclusive.
   insitu_start = info.insitu_window(1);
   insitu_end = info.insitu_window(2);

   switch tier
      case "smoke"
         % Single representative snow-season year. Pick the second
         % full year of the site coverage so the window starts after
         % any spin-up gaps in the upstream files.
         smoke_year = insitu_start + 1;
         window_start = datetime(smoke_year, ...
            10, 1, 0, 0, 0, 'TimeZone', 'UTC');
         window_end = datetime(smoke_year + 1, ...
            9, 30, 23, 0, 0, 'TimeZone', 'UTC');

      case "full"
         % Full available range as advertised by the namelist.
         window_start = datetime(insitu_start, ...
            10, 1,  0, 0, 0, 'TimeZone', 'UTC');
         window_end   = datetime(insitu_end,  ...
            9, 30, 23, 0, 0, 'TimeZone', 'UTC');

      case "custom"
         window_start = ensureUtc(custom_start);
         window_end = ensureUtc(custom_end);
   end
end

function vars = obsComparisonVariables(obs_tt)
   %OBSCOMPARISONVARIABLES Pick comparison variables from staged obs.
   %
   % Returns a string column with the canonical ordering: snow_depth_m,
   % swe_kg_m2, surface_temp_C, soil_temp_<k>_C. Variables absent from
   % the staged timetable are dropped so comparecase does not emit
   % "not_applicable" rows for never-observed variables.
   present = string(obs_tt.Properties.VariableNames);
   canonical = ["snow_depth_m"; "swe_kg_m2"; "surface_temp_C"];
   soil = present(startsWith(present, "soil_temp_"));
   ordered = [canonical(ismember(canonical, present)); sort(soil(:))];
   vars = ordered;
end

function note = siteCaseNote(info, tier)
   %SITECASENOTE Short manifest note describing the staged case.
   note = sprintf( ...
      'ESM-SnowMIP %s reference site (%s, %s). Tier=%s.', ...
      info.long_name, info.location, info.sitename, tier);
end

function t = ensureUtc(t)
   if isstring(t) || ischar(t)
      t = datetime(t, 'TimeZone', 'UTC');
   elseif isdatetime(t) && isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   end
end
