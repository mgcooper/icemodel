function [met_outdir, userdata_outdir] = rcmArtifactOutputDirs( ...
      met_outdir, userdata_outdir)
   %RCMARTIFACTOUTPUTDIRS Resolve shared default RCM artifact output roots.
   %
   %  [met_outdir, userdata_outdir] = ...rcmArtifactOutputDirs( ...
   %     met_outdir, userdata_outdir)
   %
   % Explicit roots pass through unchanged. Blank roots resolve exactly as the
   % shared met/userdata writers do, keeping staging discovery and prior-leg
   % preservation on the same filesystem contract.

   arguments
      met_outdir (1, 1) string = ""
      userdata_outdir (1, 1) string = ""
   end

   if met_outdir == ""
      met_outdir = string(fullfile(icemodel.getpath('input'), 'met'));
   end
   if userdata_outdir == ""
      userdata_outdir = string(icemodel.getpath('userdata'));
   end
end
