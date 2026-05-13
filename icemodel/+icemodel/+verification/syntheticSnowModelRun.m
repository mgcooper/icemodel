function [ice1, ice2, opts] = syntheticSnowModelRun(opts)
   %SYNTHETICSNOWMODELRUN Return snow-model-like icemodel outputs.
   %
   %  [ice1, ice2, opts] = icemodel.verification.syntheticSnowModelRun(opts)
   %
   % This verification-only model stand-in lets the suite exercise the
   % icemodel(opts) call path and the ICE1/ICE2 adapter before production
   % snow physics is implemented. It reads the same staged forcing and
   % target artifacts used by the verification suite, perturbs target
   % snow variables by the offsets / scales in opts.verification_*, and
   % stores the result as icemodel output fields.
   %
   % WARNING - this is NOT a real snow-model run
   % -------------------------------------------
   % The output is a deliberate +5 % / +0.25 K / +2 cm perturbation of
   % the staged targets, NOT an icemodel-physics output. Use this path
   % only to validate the verification adapter / runner end-to-end.
   % When run_snow_verification_suite('run_icemodel', true) prints
   % a +5 % storage bias against staged targets, that bias IS the
   % synthetic perturbation, not a model error.
   %
   % Retirement of this hook is tracked under icemodel-tk6.7. The
   % Colbeck verification path uses real icemodel.column.infiltration
   % outputs and does not route through this stand-in.

   manifest = opts.verification_case_manifest;

   % Stand-in candidates are derived directly from the staged targets;
   % no forcing read is needed.
   targets = icemodel.verification.helpers.loadArtifact( ...
      manifest.evaluation_path, "targets");

   % Multi-source schema fallback: when an evaluation.mat carries
   % multiple target sources keyed at the top level (e.g. colbeck1976
   % carries numerical_summa and analytical_clark2017), pick the
   % default numerical_summa source for the synthetic perturbation
   % path. compareSolutions handles the multi-source comparison.
   if ~isfield(targets, 'format') && isfield(targets, 'numerical_summa')
      targets = targets.numerical_summa;
   end

   switch targets.format
      case "timeseries"
         [ice1, ice2] = timeseriesCandidateOutput(targets.data, opts);
      case "experiment_bundle"
         [ice1, ice2] = experimentBundleCandidateOutput( ...
            targets.experiments, opts);
      otherwise
         error('unsupported verification target format: %s', targets.format)
   end

   opts.verification_candidate_source = "synthetic_snow_model_run";
end

function [ice1, ice2] = timeseriesCandidateOutput(targets, opts)
   %TIMESERIESCANDIDATEOUTPUT Build ICE1/ICE2 for ESM-SnowMIP site cases.

   data = targets;
   variable_names = string(data.Properties.VariableNames);

   if any(variable_names == "snow_depth_m")
      data.snow_depth_m = nonnegativeFinite( ...
         data.snow_depth_m + opts.verification_snow_depth_offset_m);
   end
   if any(variable_names == "swe_kg_m2")
      data.swe_kg_m2 = nonnegativeFinite( ...
         data.swe_kg_m2 .* opts.verification_swe_scale);
   end
   if any(variable_names == "surface_temp_C")
      data.surface_temp_C = data.surface_temp_C ...
         + opts.verification_surface_temp_offset_C;
   end

   ice1 = struct();
   ice1.Time = data.Properties.RowTimes;
   ice1.verification_format = "timeseries";

   % Store both verification-native names and core-adjacent names where they
   % exist today. The adapter prefers native snow names and falls back to Tsfc.
   for name = variable_names(:)'
      ice1.(name) = data.(name);
   end
   if any(variable_names == "surface_temp_C")
      Tf = icemodel.physicalConstant('Tf');
      ice1.Tsfc = data.surface_temp_C + Tf;
   end
   if any(variable_names == "snow_depth_m")
      ice1.snow_depth = data.snow_depth_m;
   end
   if all(ismember(["snow_depth_m", "swe_kg_m2"], variable_names))
      ice1.snow_density_kg_m3 = data.swe_kg_m2 ./ data.snow_depth_m;
      ice1.snow_density_kg_m3(data.snow_depth_m <= 0) = NaN;
   end

   ice2 = struct();
   ice2.verification_format = "timeseries";
end

function [ice1, ice2] = experimentBundleCandidateOutput(experiments, opts)
   %EXPERIMENTBUNDLECANDIDATEOUTPUT Build ICE1/ICE2 for process benchmarks.

   experiment_names = fieldnames(experiments);
   candidate_experiments = experiments;

   for n = 1:numel(experiment_names)
      name = experiment_names{n};
      data = candidate_experiments.(name);
      variable_names = string(data.Properties.VariableNames);

      if any(variable_names == "snow_liquid_water_storage_m")
         data.snow_liquid_water_storage_m = nonnegativeFinite( ...
            data.snow_liquid_water_storage_m ...
            .* opts.verification_liquid_water_scale);
      end

      candidate_experiments.(name) = data;
   end

   ice1 = struct( ...
      "verification_format", "experiment_bundle", ...
      "verification_experiments", candidate_experiments);
   ice2 = struct("verification_format", "experiment_bundle");
end

function values = nonnegativeFinite(values)
   %NONNEGATIVEFINITE Clamp finite negative values while preserving NaNs.
   values(isfinite(values) & values < 0) = 0;
end
