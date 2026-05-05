function candidate = runIcemodelSnowCandidate(case_manifest, kwargs)
   %RUNICEMODELSNOWCANDIDATE Run icemodel and return a verification candidate.
   %
   %  candidate = icemodel.verification.runIcemodelSnowCandidate(case_manifest)
   %  candidate = icemodel.verification.runIcemodelSnowCandidate(case_manifest, ...
   %     startdate=..., enddate=...)
   %
   % This is the bridge agents should use while developing the snow model. It
   % calls the normal icemodel entry point, receives ICE1/ICE2, then converts
   % those outputs into the same candidate bundle consumed by comparecase.
   %
   % The function constructs opts via icemodel.verification.helpers.caseSetopts
   % using the case_id as both sitename and forcings, so the standard chain
   % (configureRun + createMetFileNames + loadmet) resolves the staged met
   % file at <ICEMODEL_INPUT_PATH>/met/met_<case>_<case>_<year>_1hr.mat
   % without verification-only branches. This is the same scaffolding the
   % regression and perf suites use; once the snow model is operational the
   % synthetic-snow hook below retires and this function becomes a thin
   % wrapper around icemodel(opts).
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
      kwargs.startdate = NaT('TimeZone', 'UTC')
      kwargs.enddate   = NaT('TimeZone', 'UTC')
   end

   % Build standard model options. caseSetopts uses case_id as sitename
   % and forcings, derives simyears from the comparison window, and
   % populates opts.startdate / opts.enddate so the standard loadmet
   % path resolves the staged met file and subsets to the comparison
   % window.
   opts = icemodel.verification.helpers.caseSetopts(case_manifest, ...
      startdate=kwargs.startdate, enddate=kwargs.enddate);

   % Activate the synthetic-snow hook. These fields are deliberately not
   % part of setopts: they are a narrow verification-only request consumed
   % by icemodel.verification.syntheticSnowModelRun and ignored by normal
   % model configuration. They retire together with the hook itself once
   % production snow physics lands.
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
