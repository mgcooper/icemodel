function candidate = runIcemodelSnowCandidate(case_manifest, kwargs)
   %RUNICEMODELSNOWCANDIDATE Run icemodel and return a verification candidate.
   %
   %  candidate = icemodel.verification.runIcemodelSnowCandidate(case_manifest)
   %
   % This is the bridge agents should use while developing the snow model. It
   % calls the normal icemodel entry point, receives ICE1/ICE2, then converts
   % those outputs into the same candidate bundle consumed by comparecase.
   %
   % Until production snow physics exists, icemodel is run with an explicit
   % verification synthetic-snow hook that perturbs the staged targets. That
   % proves the execution and comparison path now, without forcing placeholder
   % snow logic into the core solver.

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
