function candidate = runIcemodelSnowCandidate(case_manifest, kwargs)
   %RUNICEMODELSNOWCANDIDATE Run icemodel and return a verification candidate.
   %
   %  candidate = icemodel.verification.runIcemodelSnowCandidate(case_manifest)
   %
   % This is the bridge agents should use while developing the snow model. It
   % calls the normal icemodel entry point, receives ICE1/ICE2, then converts
   % those outputs into the same candidate bundle consumed by comparecase.
   %
   % WARNING - synthetic-candidate path
   % ----------------------------------
   % Until production snow physics exists, icemodel is run with an explicit
   % verification synthetic-snow hook (icemodel.verification.syntheticSnowModelRun)
   % that perturbs the staged targets. The hook applies these hard-coded
   % offsets / scales by default:
   %
   %     snow_depth_offset_m   = 0.02   (+2 cm)
   %     swe_scale             = 1.05   (+5 %)
   %     surface_temp_offset_C = 0.25   (+0.25 K)
   %     liquid_water_scale    = 1.05   (+5 %)
   %
   % The resulting candidate bundle is NOT a real model output. The +5 %
   % storage / SWE biases visible in the suite metrics are the synthetic
   % perturbation, not a model error. This bridge proves the execution and
   % comparison path now, without forcing placeholder snow logic into the
   % core solver. Retirement of the synthetic-snow hook is tracked under
   % icemodel-tk6.7.
   %
   % The Colbeck verification path (icemodel.verification.colbeck.runCase /
   % compareSolutions) does NOT use this synthetic-candidate path; it
   % produces real numerical and analytical IceModel candidates from the
   % first-class infiltration kernel.

   arguments
      case_manifest (1, 1) struct
      kwargs.snow_depth_offset_m (1, 1) double = 0.02
      kwargs.swe_scale (1, 1) double = 1.05
      kwargs.surface_temp_offset_C (1, 1) double = 0.25
      kwargs.liquid_water_scale (1, 1) double = 1.05
   end

   simyears = caseSimulationYears(case_manifest);
   opts = icemodel.setopts('icemodel', 'verification', simyears, 'kanm', ...
      [], [], case_manifest.case_id, false, false);

   % These fields are deliberately not part of setopts. They are a narrow
   % verification-only request consumed by icemodel.verification.* and ignored
   % by normal model configuration.
   opts.verification_synthetic_snow = true;
   opts.verification_case_manifest = case_manifest;
   opts.verification_snow_depth_offset_m = kwargs.snow_depth_offset_m;
   opts.verification_swe_scale = kwargs.swe_scale;
   opts.verification_surface_temp_offset_C = kwargs.surface_temp_offset_C;
   opts.verification_liquid_water_scale = kwargs.liquid_water_scale;

   [ice1, ice2, opts] = icemodel(opts);
   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, case_manifest);
end

function simyears = caseSimulationYears(case_manifest)
   %CASESIMULATIONYEARS Infer icemodel simulation years from staged targets.

   targets = icemodel.verification.helpers.loadArtifact( ...
      case_manifest.evaluation_path, "targets");

   % Some cases stage multiple target sources keyed under one
   % evaluation.mat (e.g. colbeck1976 carries numerical_summa and
   % analytical_clark2017). Pick the default numerical_summa source
   % so the format dispatch below sees the inner experiment_bundle.
   if ~isfield(targets, 'format') && isfield(targets, 'numerical_summa')
      targets = targets.numerical_summa;
   end

   switch targets.format
      case "timeseries"
         time = targets.data.Properties.RowTimes;
      case "experiment_bundle"
         experiment_names = fieldnames(targets.experiments);
         first_experiment = targets.experiments.(experiment_names{1});
         time = first_experiment.Properties.RowTimes;
      otherwise
         error('unsupported verification target format: %s', targets.format)
   end

   simyears = unique(year(time));
end
