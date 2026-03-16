function filepath = restartfile(opts, simyear)
%RESTARTFILE Return the canonical restart-state file path for a run/year.
%
%  filepath = icemodel.restartfile(opts, simyear)
%
% Restart files live under opts.pathrestart, which by default is:
%  opts.pathoutput/restart

   arguments
      opts struct
      simyear (1, 1) double {mustBeInteger, mustBePositive}
   end

   filepath = fullfile(opts.pathrestart, ...
      sprintf('restart_%s_%d.mat', opts.casename, simyear));
end
