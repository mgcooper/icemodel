function [forcing, targets, reference] = ...
      buildLaughTestsArtifacts(source_dir, case_id)
   %BUILDLAUGHTESTSARTIFACTS Build one Laugh-Tests evaluation/reference bundle.
   %
   %  [forcing, targets, reference] = ...
   %     icemodel.verification.setup.buildLaughTestsArtifacts(source_dir)
   %
   % Raw Laugh-Tests reads and source-specific normalization live here so
   % importLaughTests remains orchestration-only. The returned structs are not
   % written; the importer owns atomic evaluation/reference persistence.

   arguments
      source_dir (1, 1) string
      case_id (1, 1) string = "colbeck1976"
   end

   % Dispatch only truly case-specific construction at this builder boundary.
   switch case_id
      case "colbeck1976"
         [forcing, targets, reference] = buildColbeckArtifacts(source_dir);
      otherwise
         valid_cases = icemodel.verification.namelists.caseid("laugh_tests");
         error('icemodel:verification:buildLaughTestsArtifacts:unsupportedCase', ...
            'unsupported Laugh-Tests case %s. Valid cases: %s', ...
            case_id, strjoin(valid_cases, ', '))
   end
end

function [forcing, targets, reference] = buildColbeckArtifacts(source_dir)
   %BUILDCOLBECKARTIFACTS Normalize the Colbeck process benchmark.

   % Resolve the exact forcing and frozen validation paths consumed by the
   % canonical Colbeck case.
   forcing_file = fullfile(source_dir, 'test_cases', 'input_data', ...
      'colbeck1976', 'colbeck1976_forcing.nc');
   validation_root = fullfile(source_dir, 'validation_data', ...
      'm2_mac_Sept23', 'colbeck1976');

   % Map raw forcing channels to icemodel names. Precipitation changes from
   % millimetres per second to metres per second.
   Time = readLaughTime(forcing_file, 'time');
   tair = readOneVariable(forcing_file, 'airtemp');
   psfc = readOneVariable(forcing_file, 'airpres');
   qair = readOneVariable(forcing_file, 'spechum');
   wspd = readOneVariable(forcing_file, 'windspd');
   swd = readOneVariable(forcing_file, 'SWRadAtm');
   lwd = readOneVariable(forcing_file, 'LWRadAtm');
   ppt = readOneVariable(forcing_file, 'pptrate') ./ 1000;
   rh = icemodel.vapor.relative_humidity_from_specific_humidity( ...
      qair, psfc, tair);

   % Laugh-Tests omits albedo and snow depth. The Colbeck SWRadAtm channel is
   % zero, so the upstream SUMMA albedoMax value is a traceable placeholder.
   albedo = 0.84 + zeros(size(tair));
   snow_depth = nan(size(tair));
   forcing_tt = timetable(Time, tair, swd, lwd, albedo, wspd, rh, psfc, ...
      ppt, snow_depth, 'VariableNames', {'tair', 'swd', 'lwd', 'albedo', ...
      'wspd', 'rh', 'psfc', 'ppt', 'snow_depth'});
   forcing = struct( ...
      'format', 'timeseries', ...
      'data', forcing_tt, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'albedo_policy', ['constant 0.84 from Colbeck SUMMA albedoMax; ' ...
      'SWRadAtm is zero so albedo is not an active process forcing']
      'snow_depth_policy', ...
      'NaN placeholder because Laugh-Tests provides no snow-depth forcing channel'}));

   % Build one experiment timetable for each frozen SUMMA validation output.
   exp_ids = ["exp1"; "exp2"; "exp3"];
   experiment_rows = cell(numel(exp_ids), 1);
   for n = 1:numel(exp_ids)
      experiment_file = fullfile(validation_root, sprintf( ...
         'colbeck1976-%s_G1-1_timestep.nc', exp_ids(n)));
      experiment_rows{n} = buildColbeckExperiment(experiment_file);
   end
   experiments = cell2struct(experiment_rows, cellstr(exp_ids), 1);

   % Build the analytical Clark 2017 bundle in the same experiment schema.
   def = icemodel.verification.colbeck.caseDefinition();
   analytical_rows = cell(numel(exp_ids), 1);
   for n = 1:numel(exp_ids)
      analytical_rows{n} = analyticalExperimentTimetable(exp_ids(n), def);
   end
   analytical_experiments = cell2struct(analytical_rows, ...
      cellstr(exp_ids), 1);

   % Return the two evaluation targets and frozen reference bundle.
   targets = struct( ...
      'numerical_summa', struct( ...
      'format', 'experiment_bundle', ...
      'experiments', experiments, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'source', 'frozen_summa_validation_output'
      'validation_bundle', 'm2_mac_Sept23'})), ...
      'analytical_clark2017', struct( ...
      'format', 'experiment_bundle', ...
      'experiments', analytical_experiments, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'source', 'icemodel.verification.colbeck.analyticalSolution'
      'method', 'Clark 2017 wetting-front / kinematic-wave'})));
   reference = struct( ...
      'format', 'experiment_bundle', ...
      'experiments', experiments, ...
      'metadata', icemodel.verification.setup.metadataStruct({ ...
      'reference_kind', 'frozen_summa_validation_output'
      'validation_bundle', 'm2_mac_Sept23'}));
end

function tt = analyticalExperimentTimetable(experiment_name, def)
   %ANALYTICALEXPERIMENTTIMETABLE Wrap analyticalSolution for one experiment.
   sol = icemodel.verification.colbeck.analyticalSolution( ...
      experiment_name, def);
   tt = timetable(sol.time_datetime, ...
      sol.snow_liquid_water_storage_m, sol.bottom_outflow_mps, ...
      'VariableNames', {'snow_liquid_water_storage_m', 'bottom_outflow_mps'});
end

function tt = buildColbeckExperiment(pathname)
   %BUILDCOLBECKEXPERIMENT Derive scalar diagnostics from one Colbeck run.
   %
   % Frozen SUMMA output carries layer profiles. The staged reference keeps
   % total snow liquid-water storage and the bottom-outflow proxy.

   % Read raw profile and scalar variables from one validation file.
   time = readLaughTime(pathname, 'time');
   liq = squeeze(double(ncread(pathname, 'mLayerVolFracLiq')));
   depth = squeeze(double(ncread(pathname, 'mLayerDepth')));
   iface = squeeze(double(ncread(pathname, 'iLayerHeight')));
   outflow = squeeze(double(ncread(pathname, 'scalarRainPlusMelt')));

   % Integrate liquid storage only through snow layers above the soil surface.
   ntime = numel(time);
   snow_liquid_water_storage_m = zeros(ntime, 1);
   for n = 1:ntime
      layer_depth = depth(:, n);
      layer_liq = liq(:, n);
      interfaces = iface(:, n);
      mid_height = 0.5 * (interfaces(1:end-1) + interfaces(2:end));
      snow_mask = layer_depth > 0 & mid_height < 0;
      if ~any(snow_mask)
         continue
      end
      snow_liquid_water_storage_m(n) = ...
         sum(layer_liq(snow_mask) .* layer_depth(snow_mask));
   end

   % Return the generic experiment timetable consumed by comparison drivers.
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
   time = double(ncread(pathname, varname));
   unit = string(ncreadatt(pathname, varname, 'units'));
   tref = datetime(extractBefore(extractAfter(unit, 'since '), ' -'), ...
      'InputFormat', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
   time = tref + seconds(time);
   time = time(:);
end
