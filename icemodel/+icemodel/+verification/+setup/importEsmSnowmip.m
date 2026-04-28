function manifest = importEsmSnowmip(source_dir, kwargs)
   %IMPORTESMSNOWMIP Stage curated ESM-SnowMIP smoke fixtures.
   %
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir)
   %  manifest = icemodel.verification.setup.importEsmSnowmip(source_dir, ...
   %     overwrite=true)
   %
   % Inputs
   %  source_dir                 Folder containing ESM-SnowMIP NetCDF files.
   %  evaluation_data_root       Base evaluation-data root to stage into.
   %  icemodel_config_casename   Config casename used when evaluation_data_root
   %                             is blank. Defaults to the repo test/demo data.
   %  case_ids                   ESM-SnowMIP case ids to stage.
   %  overwrite                  Refresh setup artifacts when true; protect
   %                             existing staged data when false.
   %
   % Outputs
   %  manifest   Family manifest struct also written to manifest.json.
   %
   % Role
   %  Setup/update tooling. This function creates or refreshes staged data
   %  under demo/data/eval/snow and is not part of normal verification runs.
   %
   % See also icemodel.verification.setup.importLaughTests

   arguments
      source_dir (1, :) string
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.case_ids (1, :) string = ...
         icemodel.verification.namelists.caseid("esm_snowmip")
      kwargs.overwrite (1, 1) logical = false
   end

   % Name the source family and runnable cases once. dataset_family is the
   % staged source folder/manifest family; case_ids are the site cases inside
   % that family.
   dataset_family = "esm_snowmip";
   case_ids = reshape(kwargs.case_ids, 1, []);

   % Convert case ids into source-file/window specs. This keeps the import loop
   % data-driven without hiding site-specific provenance. Validate before any
   % setup directories are created so bad options fail without side effects.
   specs = selectSpecs(defaultSpecs(), case_ids);
   case_entries = cell(numel(specs), 1);

   % Resolve the path to the snow evaluation data folder
   snow_data_root = icemodel.verification.helpers.snowDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);

   % Resolve the path to the dataset family sub-folder
   %   <snow_data_root>/esm_snowmip
   family_root = fullfile(snow_data_root, dataset_family);
   icemodel.helpers.ensureDirExists(family_root);

   % Resolve the path to the dataset family manifest
   %   <snow_data_root>/esm_snowmip/manifest.json
   manifest_file = fullfile(family_root, "manifest.json");

   % Stage each selected site using the same write order as importLaughTests:
   % prepare root, build artifacts, write artifacts, then build manifest entry.
   for i = 1:numel(specs)
      spec = specs(i);

      % Resolve the path to the case id sub-folder
      %   <snow_data_root>/esm_snowmip/<case_id>
      case_root = fullfile(family_root, spec.case_id);

      % Resolve the exact staged paths produced for this case:
      %   <snow_data_root>/esm_snowmip/<case_id>/forcing.mat
      %   <snow_data_root>/esm_snowmip/<case_id>/evaluation.mat
      %   <snow_data_root>/esm_snowmip/<case_id>/reference.mat
      forcing_output_file = fullfile(case_root, "forcing.mat");
      evaluation_output_file = fullfile(case_root, "evaluation.mat");
      reference_output_file = fullfile(case_root, "reference.mat");

      % prepareCaseRoot owns the overwrite guard. From here down, this importer
      % can write the resolved output files directly without repeating per-file
      % existence checks.
      icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);

      [forcing, targets, reference, case_values] = buildEsmArtifacts( ...
         source_dir, spec);
      save(forcing_output_file, 'forcing');
      save(evaluation_output_file, 'targets');
      save(reference_output_file, 'reference');

      case_entries{i} = icemodel.verification.setup.makeCaseManifestEntry( ...
         case_values);
   end

   % Keep data provenance values explicit.
   source_doi = "10.1594/PANGAEA.897575";
   source_url = "https://doi.org/10.1594/PANGAEA.897575";
   source_version = "ESM-SnowMIP_all.zip (2019-01-23)";
   retrieval_date = string(datetime('today'));

   % Build and write the family manifest after all selected case artifacts exist
   % so the manifest never points at a partially staged family.
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      dataset_family, source_doi, source_url, source_version, retrieval_date, ...
      vertcat(case_entries{:}));
   icemodel.verification.setup.writeManifest(manifest_file, manifest);
end

function specs = defaultSpecs()
   %DEFAULTSPECS Return curated ESM-SnowMIP smoke-case definitions.

   % Each spec is the site-specific setup contract: source files, comparison
   % window, display name, and variables retained for smoke verification. This
   % list stays local because these fields are owned only by this importer and
   % are not part of the public manifest schema.
   specs = struct( ...
      'case_id', {"cdp", "wfj"}, ...
      'site_id', {"cdp", "wfj"}, ...
      'site_name', {"Col de Porte", "Weissfluhjoch"}, ...
      'met_file', {"met_insitu_cdp_1994_2014.nc", "met_insitu_wfj_1996_2016.nc"}, ...
      'obs_file', {"obs_insitu_cdp_1994_2014.nc", "obs_insitu_wfj_1996_2016.nc"}, ...
      'window_start', {datetime(2005, 10, 1, 0, 0, 0, 'TimeZone', 'UTC'), ...
      datetime(2000, 10, 1, 0, 0, 0, 'TimeZone', 'UTC')}, ...
      'window_end', {datetime(2006, 7, 1, 0, 0, 0, 'TimeZone', 'UTC'), ...
      datetime(2001, 7, 1, 0, 0, 0, 'TimeZone', 'UTC')}, ...
      'comparison_variables', {["snow_depth_m", "swe_kg_m2", "surface_temp_C", ...
      "soil_temp_1_C", "soil_temp_2_C", "soil_temp_3_C"], ...
      ["snow_depth_m", "swe_kg_m2", "surface_temp_C", "soil_temp_1_C"]});
end

function specs = selectSpecs(all_specs, case_ids)
   %SELECTSPECS Filter curated specs by requested case id.

   % Validate against the shared namelist so importer options and workflow
   % filters use the same canonical case catalog.
   valid = icemodel.verification.namelists.caseid("esm_snowmip");
   invalid = setdiff(reshape(case_ids, [], 1), valid);
   if ~isempty(invalid)
      error('unsupported ESM-SnowMIP verification cases: %s', ...
         strjoin(invalid, ', '))
   end

   % Preserve the default spec order while retaining only requested cases.
   keep = ismember([all_specs.case_id], case_ids);
   specs = all_specs(keep);
   if isempty(specs)
      error('no ESM-SnowMIP specs selected for %s', strjoin(case_ids, ', '))
   end
end

function [forcing, targets, reference, case_values] = buildEsmArtifacts( ...
      source_dir, spec)
   %BUILDESMARTIFACTS Normalize one ESM-SnowMIP site window.

   % Resolve the exact ESM-SnowMIP source files consumed by this site:
   %   <source_dir>/<spec.met_file>
   %   <source_dir>/<spec.obs_file>
   metfile = fullfile(source_dir, spec.met_file);
   obsfile = fullfile(source_dir, spec.obs_file);

   % Resolve the requested smoke-season time window before reading variables so
   % every forcing and target channel is subset with the same index.
   Time = readNetcdfTime(metfile, 'time');
   idx = Time >= spec.window_start & Time <= spec.window_end;
   Time = Time(idx);

   % Load the source forcing vectors once, then subset the chosen smoke season.
   tair = double(ncread(metfile, 'Tair'));
   swdown = double(ncread(metfile, 'SWdown'));
   lwdown = double(ncread(metfile, 'LWdown'));
   wind = double(ncread(metfile, 'Wind'));
   psfc = double(ncread(metfile, 'Psurf'));
   qair = double(ncread(metfile, 'Qair'));
   rainf = double(ncread(metfile, 'Rainf'));
   snowf = double(ncread(metfile, 'Snowf'));

   % Load observation channels before deriving forcing fields because observed
   % snow depth and albedo both influence the staged forcing timetable.
   snd_auto = readObs(obsfile, 'snd_auto');
   snd_man = readObs(obsfile, 'snd_man');
   swe_auto = optionalObs(obsfile, 'snw_auto', numel(tair));
   swe_man = optionalObs(obsfile, 'snw_man', numel(tair));
   ts = readObs(obsfile, 'ts');
   albs = readObs(obsfile, 'albs');
   tsl = readObs(obsfile, 'tsl');
   sdepth = double(ncread(obsfile, 'sdepth'));

   % Merge primary/fallback observations and derive icemodel-native forcing
   % variables for the selected window.
   snow_depth = preferPrimary(snd_auto(idx), snd_man(idx));
   swe = preferPrimary(swe_auto(idx), swe_man(idx));
   albedo = buildAlbedo(albs(idx), snow_depth);
   rh = icemodel.vapor.relative_humidity_from_specific_humidity( ...
      qair(idx), psfc(idx), tair(idx));
   ppt = (rainf(idx) + snowf(idx)) ./ 1000;

   % Create a forcing timetable using icemodel-native variable names.
   forcing_tt = timetable(Time, tair(idx), swdown(idx), lwdown(idx), albedo, ...
      wind(idx), rh, psfc(idx), ppt, snow_depth, ...
      'VariableNames', {'tair', 'swd', 'lwd', 'albedo', 'wspd', 'rh', ...
      'psfc', 'ppt', 'snow_depth'});

   % Package into a struct.
   forcing = struct( ...
      'format', 'timeseries', ...
      'data', forcing_tt, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'albedo_policy', 'observed albedo with snow/ground fallback'}));

   % Assemble observed targets. Soil temperatures are expanded into named
   % columns so comparison variables can select depths explicitly.
   observations = buildObservationTimetable(Time, snow_depth, swe, ts(idx), ...
      tsl(:, idx)');

   targets = struct( ...
      'format', 'timeseries', ...
      'data', observations, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'soil_depths_m', sdepth(:)
      'snow_depth_source', 'snd_auto with snd_man fallback'
      'swe_source', 'snw_auto with snw_man fallback when available'}));

   % Until the snow model exists, the smoke reference duplicates observations so
   % comparecase and plotcase can run in a fresh clone.
   reference = struct( ...
      'format', 'timeseries', ...
      'data', observations, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'reference_kind', 'smoke_observed_reference'
      'notes', ['Observed series duplicated as the initial smoke reference ' ...
      'so the staged harness can run before a snow model exists.']}));

   % Return the same canonical case-value contract used by importLaughTests.
   % The second value is case_type: this is a field-site case, not a
   % dataset-family id.
   observation_variables = icemodel.verification.setup.metadataStruct({ ...
      'snow_depth_m', 'snd_auto|snd_man'
      'swe_kg_m2', 'snw_auto|snw_man'
      'surface_temp_C', 'ts'
      'soil_depths_m', sdepth(:)});
   case_values = { ...
      spec.case_id
      'esm_site'
      'smoke'
      spec.site_id
      spec.site_name
      fullfile(spec.case_id, 'forcing.mat')
      fullfile(spec.case_id, 'evaluation.mat')
      fullfile(spec.case_id, 'reference.mat')
      'hourly'
      struct('start', string(spec.window_start), ...
      'end', string(spec.window_end))
      spec.comparison_variables
      observation_variables
      ['Forcing albedo is interpolated from observed albedo and falls back ' ...
      'to snow/ground constants where observations are missing. The smoke ' ...
      'reference duplicates the staged observations by design.']};
end

function observations = buildObservationTimetable(Time, snow_depth, swe, ...
      surface_temp, soil_temp)
   %BUILDOBSERVATIONTIMETABLE Assemble the staged target timetable.

   % Generate one soil-temperature variable per available observation depth.
   n_soil = size(soil_temp, 2);
   soil_names = "soil_temp_" + string(1:n_soil)' + "_C";
   variable_names = ["snow_depth_m"; "swe_kg_m2"; "surface_temp_C"; soil_names];
   variable_values = [{snow_depth, swe, surface_temp}, num2cell(soil_temp, 1)];

   % Keep all observations in one timetable so plotting and comparison code can
   % remain dataset-family agnostic.
   observations = timetable(Time, variable_values{:}, ...
      'VariableNames', cellstr(variable_names));
end

function time = readNetcdfTime(pathname, varname)
   %READNETCDFTIME Read an ESM-SnowMIP time coordinate into UTC datetime.

   % ESM-SnowMIP forcing files use hours since a UTC reference timestamp.
   raw = double(ncread(pathname, varname));
   units = string(ncreadatt(pathname, varname, 'units'));
   tref = datetime(extractAfter(units, 'hours since '), ...
      'InputFormat', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
   time = tref + hours(raw);
   time = time(:);
end

function values = readObs(pathname, varname)
   %READOBS Read one observed NetCDF variable and replace fill values with NaN.

   values = double(ncread(pathname, varname));
   values(values <= -900) = NaN;
end

function values = optionalObs(pathname, varname, ntime)
   %OPTIONALOBS Read an optional observed variable or return NaNs.

   % Some sites lack automatic SWE channels. Returning NaNs preserves the
   % target schema while allowing manual observations to fill gaps.
   info = ncinfo(pathname);
   if any(string({info.Variables.Name}) == varname)
      values = readObs(pathname, varname);
   else
      values = nan(ntime, 1);
   end
end

function merged = preferPrimary(primary, fallback)
   %PREFERPRIMARY Merge primary/fallback observation vectors.

   % Primary automated observations are retained unless missing; manual
   % observations fill only those gaps.
   merged = primary;
   replace = ~isfinite(merged) & isfinite(fallback);
   merged(replace) = fallback(replace);
end

function albedo = buildAlbedo(raw_albedo, snow_depth)
   %BUILDALBEDO Build a continuous albedo series for staged forcing files.

   % Start from observed albedo and remove out-of-range fill values.
   albedo = raw_albedo;
   albedo(albedo < 0 | albedo > 1) = NaN;
   albedo = fillmissing(albedo, 'linear', 'EndValues', 'nearest');

   % Fill remaining gaps with simple snow/ground constants. This keeps forcing
   % continuous while documenting that this is setup normalization, not a snow
   % albedo model.
   snow_mask = isfinite(snow_depth) & snow_depth > 0.05;
   albedo(~isfinite(albedo) & snow_mask) = 0.85;
   albedo(~isfinite(albedo)) = 0.45;
   albedo = min(max(albedo, 0.05), 0.95);
end
