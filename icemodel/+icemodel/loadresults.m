function [ice1, ice2, met] = loadresults(opts, varargin)
   %LOADRESULTS Load saved yearly output files and concatenate them if needed.
   %
   %  [ice1, ice2, met] = icemodel.loadresults(opts)
   %  [ice1, ice2, met] = icemodel.loadresults(opts, simyears)
   %
   % Without SIMYEARS, this function loads every post-spinup output year. For
   % more than one year, it concatenates the yearly postprocessed output files
   % in time, so the returned data match the no-save multi-year run.point path.

   % If the simulation year is not specified, load all saved output years.
   if nargin < 2
      if isfield(opts, 'output_years') && ~isempty(opts.output_years)
         simyears = opts.output_years;
      else
         % A saved OPTS struct can lack the OUTPUT_YEARS field.
         simyears = icemodel.outputYears(opts);
      end
   else
      simyears = varargin{1};
      if ischar(simyears) || isstring(simyears)
         simyears = str2double(simyears);
      end
      simyears = simyears(:);
   end

   % Load the met data for the requested years and post-process it.
   met = icemodel.loadmet(opts);
   met = met(ismember(year(met.Time), simyears), :);
   met = icemodel.processmet(met, newTimeStep="hourly");

   % Load the simulation results, concatenating yearly files if needed.
   ice1 = [];
   ice2 = [];
   for n = 1:numel(simyears)
      thisyear = int2str(simyears(n));
      try
         S = load(fullfile(opts.pathoutput, thisyear, ...
            ['ice1_' opts.casename]), 'ice1');
         T = load(fullfile(opts.pathoutput, thisyear, ...
            ['ice2_' opts.casename]), 'ice2');
      catch e
         if strcmp(e.identifier, 'MATLAB:load:couldNotReadFile')
            rethrow(e)
         end
      end

      [ice1, ice2] = icemodel.concatoutput(ice1, ice2, S.ice1, T.ice2);
   end
end
