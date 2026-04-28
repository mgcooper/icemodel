function [ice1, ice2, opts] = syntheticSnowModelRun(opts)
   %SYNTHETICSNOWMODELRUN Return snow-model-like icemodel outputs.
   %
   %  [ice1, ice2, opts] = icemodel.verification.syntheticSnowModelRun(opts)
   %
   % This verification-only model stand-in lets the suite exercise the real
   % icemodel(opts) call and ICE1/ICE2 adapter before production snow physics is
   % implemented. It reads the same staged forcing and target artifacts used by
   % the verification suite, perturbs target snow variables, and stores the
   % result as icemodel output fields.

   manifest = opts.verification_case_manifest;

   % Load the forcing artifact so the verification run is tied to the case's
   % staged inputs now, even though this temporary stand-in uses targets to
   % synthesize model-like outputs.
   forcing = loadVerificationForcing(manifest.forcing_path);
   targets = icemodel.verification.helpers.loadArtifact( ...
      manifest.evaluation_path, "targets");

   switch targets.format
      case "timeseries"
         [ice1, ice2] = timeseriesCandidateOutput(targets.data, forcing, opts);
      case "experiment_bundle"
         [ice1, ice2] = experimentBundleCandidateOutput( ...
            targets.experiments, forcing, opts);
      otherwise
         error('unsupported verification target format: %s', targets.format)
   end

   opts.verification_candidate_source = "synthetic_snow_model_run";
end

function forcing = loadVerificationForcing(pathname)
   %LOADVERIFICATIONFORCING Normalize staged forcing MAT-file shapes.

   contents = load(pathname);
   names = string(fieldnames(contents));

   if any(names == "met")
      forcing = struct("format", "timeseries", "data", contents.met);
   elseif any(names == "forcing")
      forcing = contents.forcing;
   else
      error('verification forcing file does not contain met or forcing: %s', ...
         pathname)
   end
end

function [ice1, ice2] = timeseriesCandidateOutput(targets, forcing, opts)
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
   ice1.verification_forcing_time = forcing.data.Properties.RowTimes;
   ice1.verification_format = "timeseries";

   % Store both verification-native names and core-adjacent names where they
   % exist today. The adapter prefers native snow names and falls back to Tsfc.
   for i = 1:numel(variable_names)
      name = variable_names(i);
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

function [ice1, ice2] = experimentBundleCandidateOutput(experiments, ...
      forcing, opts)
   %EXPERIMENTBUNDLECANDIDATEOUTPUT Build ICE1/ICE2 for process benchmarks.

   experiment_names = fieldnames(experiments);
   candidate_experiments = experiments;

   for i = 1:numel(experiment_names)
      name = experiment_names{i};
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
      "verification_experiments", candidate_experiments, ...
      "verification_forcing_time", forcing.data.Properties.RowTimes);
   ice2 = struct("verification_format", "experiment_bundle");
end

function values = nonnegativeFinite(values)
   %NONNEGATIVEFINITE Clamp finite negative values while preserving NaNs.

   values(isfinite(values) & values < 0) = 0;
end
