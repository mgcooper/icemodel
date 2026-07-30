function paths = esmRuntimeMetFiles(case_manifest, input_root)
   %ESMRUNTIMEMETFILES Resolve an atomic ESM case's standard runtime met paths.
   %
   %  paths = icemodel.verification.helpers.esmRuntimeMetFiles( ...
   %     case_manifest, input_root)
   %
   % ESM-SnowMIP cases intentionally omit forcing_sources, colocation, and
   % met_files because forcing and observations are staged atomically. Reuse the
   % verification runner's normal option-resolution chain so audit and plotting
   % select the same nested-or-flat artifact as an actual model run. Missing
   % files remain in PATHS; callers decide whether to report or skip them.

   % Carry the caller-selected input root into the same helper used by
   % runIcemodelSnowCandidate. configureRun then delegates filename and path
   % selection to createMetFileNames and sourceSearchDirs.
   case_manifest.input_data_root = char(input_root);
   opts = icemodel.test.helpers.setModelOptsForCase(case_manifest);
   paths = reshape(string(opts.metfname), 1, []);
end
