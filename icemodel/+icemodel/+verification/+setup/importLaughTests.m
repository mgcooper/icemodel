function manifest = importLaughTests(laugh_tests_source_dir, kwargs)
   %IMPORTLAUGHTESTS Stage selected Laugh-Tests synthetic snow benchmarks.
   %
   %  manifest = icemodel.verification.setup.importLaughTests(source_dir)
   %  manifest = icemodel.verification.setup.importLaughTests(source_dir, ...
   %     overwrite=true)
   %
   % Inputs
   %  laugh_tests_source_dir     Root of a local Laugh-Tests checkout.
   %  evaluation_data_root       Base evaluation-data root to stage into.
   %  icemodel_config_casename   Config casename used when evaluation_data_root
   %                             is blank. Defaults to the repo test/demo data.
   %  case_id                    Laugh-Tests case to stage. Currently only
   %                             "colbeck1976" is supported.
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
   % See also icemodel.verification.setup.importEsmSnowmip

   arguments
      laugh_tests_source_dir (1, :) string
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.case_id (1, 1) string = ...
         icemodel.verification.namelists.caseid("laugh_tests")
      kwargs.overwrite (1, 1) logical = false
   end

   % Name the source family and runnable case once. dataset_family is the
   % staged source folder/manifest family; case_id is the benchmark case inside
   % that family.
   dataset_family = "laugh_tests";
   case_id = kwargs.case_id;

   % Resolve the path to the snow evaluation data folder
   snow_data_root = icemodel.verification.helpers.snowDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);

   % Resolve the path to the dataset family sub-folder
   family_root = fullfile(snow_data_root, dataset_family);

   % Resolve the path to the case id folder sub-folder
   %   <snow_data_root>/laugh_tests/colbeck1976
   case_root = fullfile(family_root, case_id);

   % Resolve the exact staged paths produced by this importer. Keeping these
   % paths explicit makes it clear what the setup pass writes:
   %   <snow_data_root>/laugh_tests/manifest.json
   %   <snow_data_root>/laugh_tests/colbeck1976/forcing.mat
   %   <snow_data_root>/laugh_tests/colbeck1976/evaluation.mat
   %   <snow_data_root>/laugh_tests/colbeck1976/reference.mat
   manifest_file = fullfile(family_root, "manifest.json");
   forcing_output_file = fullfile(case_root, "forcing.mat");
   evaluation_output_file = fullfile(case_root, "evaluation.mat");
   reference_output_file = fullfile(case_root, "reference.mat");

   % Choose the case-specific artifact builder before preparing output folders.
   % This keeps unsupported case ids from creating staged-data directories.
   switch case_id
      case "colbeck1976"
         buildArtifacts = @buildColbeckArtifacts;
      otherwise
         valid_cases = icemodel.verification.namelists.caseid(dataset_family);
         error('unsupported Laugh-Tests verification case %s. Valid cases: %s', ...
            case_id, strjoin(valid_cases, ', '))
   end

   % prepareCaseRoot owns the overwrite guard. From here down, this importer
   % can write the resolved output files directly without repeating per-file
   % existence checks.
   icemodel.verification.setup.prepareCaseRoot(case_root, kwargs.overwrite);

   % Build case-specific artifacts. The builder choice above is intentionally
   % explicit so future Laugh-Tests cases can document their own source paths.
   [forcing, targets, reference, case_values] = buildArtifacts( ...
      laugh_tests_source_dir);

   % Write the three case artifacts named in the manifest entry.
   save(forcing_output_file, 'forcing');
   save(evaluation_output_file, 'targets');
   save(reference_output_file, 'reference');

   % Keep provenance values explicit instead of burying constants in the
   % makeFamilyManifest call.
   source_doi = "";
   source_url = "https://github.com/KyleKlenk/Laugh-Tests";
   source_version = "m2_mac_Sept23 validation bundle";
   retrieval_date = string(datetime('today'));

   % Build and write the family manifest after case artifacts exist so the JSON
   % never points at files that failed to stage.
   case_entry = icemodel.verification.setup.makeCaseManifestEntry(case_values);
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      dataset_family, source_doi, source_url, source_version, retrieval_date, ...
      case_entry);
   icemodel.verification.setup.writeManifest(manifest_file, manifest);
end

function [forcing, targets, reference, case_values] = buildColbeckArtifacts( ...
      laugh_tests_source_dir)
   %BUILDCOLBECKARTIFACTS Normalize the Colbeck process benchmark.

   % Resolve the exact Laugh-Tests source paths consumed by this case:
   %   <laugh_tests_source_dir>/test_cases/input_data/colbeck1976/colbeck1976_forcing.nc
   %   <laugh_tests_source_dir>/validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp1_G1-1_timestep.nc
   %   <laugh_tests_source_dir>/validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp2_G1-1_timestep.nc
   %   <laugh_tests_source_dir>/validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp3_G1-1_timestep.nc
   forcing_file = fullfile(laugh_tests_source_dir, 'test_cases', 'input_data', ...
      'colbeck1976', 'colbeck1976_forcing.nc');
   validation_root = fullfile(laugh_tests_source_dir, 'validation_data', ...
      'm2_mac_Sept23', 'colbeck1976');

   % Read forcing channels and map them onto icemodel-native meteorological
   % variable names. Precipitation is converted from mm s-1 to m s-1.
   Time = readLaughTime(forcing_file, 'time');
   tair = readOneVariable(forcing_file, 'airtemp');
   psfc = readOneVariable(forcing_file, 'airpres');
   qair = readOneVariable(forcing_file, 'spechum');
   wspd = readOneVariable(forcing_file, 'windspd');
   swd = readOneVariable(forcing_file, 'SWRadAtm');
   lwd = readOneVariable(forcing_file, 'LWRadAtm');
   ppt = readOneVariable(forcing_file, 'pptrate') ./ 1000;
   rh = icemodel.vapor.relative_humidity_from_specific_humidity(qair, psfc, tair);

   % Laugh-Tests does not provide these as forcing channels. Fixed albedo and
   % missing snow depth are explicit setup placeholders, not model assumptions.
   albedo = 0.85 + zeros(size(tair));
   snow_depth = nan(size(tair));

   % Store forcing as a timetable artifact so future model adapters can consume
   % the same layout as ESM-SnowMIP site cases.
   forcing_tt = timetable(Time, tair, swd, lwd, albedo, wspd, rh, psfc, ...
      ppt, snow_depth, 'VariableNames', {'tair', 'swd', 'lwd', 'albedo', ...
      'wspd', 'rh', 'psfc', 'ppt', 'snow_depth'});
   forcing = struct( ...
      'format', 'timeseries', ...
      'data', forcing_tt, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'albedo_policy', 'constant 0.85 placeholder'}));

   % Build one experiment timetable for each frozen SUMMA validation output.
   exp_ids = ["exp1"; "exp2"; "exp3"];
   experiment_rows = cell(numel(exp_ids), 1);
   for i = 1:numel(exp_ids)
      exp_id = exp_ids(i);
      experiment_file = fullfile(validation_root, ...
         sprintf('colbeck1976-%s_G1-1_timestep.nc', exp_id));
      experiment_rows{i} = buildColbeckExperiment(experiment_file);
   end
   experiments = cell2struct(experiment_rows, cellstr(exp_ids), 1);

   % Document the intentional diagnostic scope. These frozen SUMMA outputs are
   % useful process targets, but they are not an analytical wetting-front target.
   targets_note = ['Wetting-front depth is intentionally omitted from the ' ...
      'staged benchmark until a true analytical front target is imported.'];
   case_note_text = ['The staged target is derived from the frozen SUMMA ' ...
      'validation outputs bundled in Laugh-Tests, not from the full ' ...
      'upstream repo layout. A true wetting-front target is intentionally ' ...
      'not staged yet because these frozen SUMMA outputs already contain ' ...
      'liquid water throughout the snow column from the first saved step, ' ...
      'so a threshold-based wetting-front diagnostic is not trustworthy.'];

   % Targets and references intentionally share the same experiment timetables
   % for the smoke lane; future model output can replace only the candidate.
   targets = struct( ...
      'format', 'experiment_bundle', ...
      'experiments', experiments, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'notes', targets_note}));

   reference = struct( ...
      'format', 'experiment_bundle', ...
      'experiments', experiments, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'reference_kind', 'frozen_summa_validation_output'
      'validation_bundle', 'm2_mac_Sept23'}));

   % Build the same canonical manifest value list used by ESM-SnowMIP import.
   % The second value is case_type: this is a synthetic process benchmark, not
   % a site-validation case and not a dataset-family id.
   observation_variables = icemodel.verification.setup.metadataStruct({ ...
      'snow_liquid_water_storage_m', ...
      'derived from mLayerVolFracLiq * mLayerDepth in snow layers'
      'bottom_outflow_mps', 'scalarRainPlusMelt'});
   case_values = { ...
      'colbeck1976'
      'synthetic_process'
      'smoke'
      'colbeck1976'
      'Colbeck 1976 synthetic snow infiltration benchmark'
      fullfile('colbeck1976', 'forcing.mat')
      fullfile('colbeck1976', 'evaluation.mat')
      fullfile('colbeck1976', 'reference.mat')
      '1 minute output / sub-hour forcing'
      struct('start', '1990-01-01 00:01:00', 'end', '1990-01-01 10:00:00')
      {'snow_liquid_water_storage_m', 'bottom_outflow_mps'}
      observation_variables
      case_note_text};
end

function tt = buildColbeckExperiment(pathname)
   %BUILDCOLBECKEXPERIMENT Derive scalar diagnostics from one Colbeck run.
   %
   % The raw SUMMA output contains full layer profiles. The staged smoke
   % reference keeps only two scalar diagnostics that are useful to an early
   % snow-model implementation: total snow liquid water storage and the bottom
   % outflow proxy.

   % Read raw profile and scalar variables from one frozen SUMMA output file.
   time = readLaughTime(pathname, 'time');
   liq = squeeze(double(ncread(pathname, 'mLayerVolFracLiq')));
   depth = squeeze(double(ncread(pathname, 'mLayerDepth')));
   iface = squeeze(double(ncread(pathname, 'iLayerHeight')));
   outflow = squeeze(double(ncread(pathname, 'scalarRainPlusMelt')));

   % Integrate liquid-water storage only through snow layers at each output
   % time. Soil layers are excluded by the interface-height mask.
   ntime = numel(time);
   snow_liquid_water_storage_m = zeros(ntime, 1);
   for k = 1:ntime
      layer_depth = depth(:, k);
      layer_liq = liq(:, k);
      interfaces = iface(:, k);

      mid_height = 0.5 * (interfaces(1:end-1) + interfaces(2:end));
      snow_mask = layer_depth > 0 & mid_height > 0;
      if ~any(snow_mask)
         continue
      end

      snow_liquid_water_storage_m(k) = ...
         sum(layer_liq(snow_mask) .* layer_depth(snow_mask));
   end

   % Return the experiment as a timetable so comparecase and plotcase can share
   % the same generic experiment-bundle code.
   Time = time(:);
   tt = timetable(Time, snow_liquid_water_storage_m, outflow(:), ...
      'VariableNames', {'snow_liquid_water_storage_m', 'bottom_outflow_mps'});
end

function out = readOneVariable(pathname, varname)
   %READONEVARIABLE Read one NetCDF variable as a column double vector.

   out = squeeze(double(ncread(pathname, varname)));
   out = out(:);
end

function time = readLaughTime(pathname, varname)
   %READLAUGHTIME Read a Laugh-Tests time coordinate into UTC datetime.

   % Most Laugh-Tests files use seconds since a timestamp that includes a
   % trailing timezone component, so parse the reference timestamp explicitly.
   time = double(ncread(pathname, varname));
   unit = string(ncreadatt(pathname, varname, 'units'));
   tref = datetime(extractBefore(extractAfter(unit, 'since '), ' -'), ...
      'InputFormat', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
   time = tref + seconds(time);
   time = time(:);
end
