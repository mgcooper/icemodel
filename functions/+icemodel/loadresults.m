function [ice1, ice2, met] = loadresults(opts, varargin)
   
   % If the simulation year is not specified, load the last one
   if nargin < 2
      thisyear = int2str(opts.simyears(end));
   else
      thisyear = varargin{1};
      if isnumeric(thisyear)
         thisyear = int2str(thisyear);
      end
   end
   
   % load the met data
   met = icemodel.loadmet(opts, 1);
   
   met = icemodel.processmet(met);
   
   % load the simulation results
   try
      load(fullfile(opts.pathoutput, thisyear, ...
         ['ice1_' opts.casename]), 'ice1');
   
      load(fullfile(opts.pathoutput, thisyear, ...
         ['ice2_' opts.casename]), 'ice2');
   catch e
      if strcmp(e.identifier, 'MATLAB:load:couldNotReadFile')
         rethrow(e)
      end
   end
end




