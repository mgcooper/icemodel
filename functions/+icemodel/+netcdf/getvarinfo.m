function [varnames, varsize] = getvarinfo(filepath, datafile, gridcell)
   %GETVARINFO Load one data file to get the shape and variable names
   %
   %
   % See also:

   arguments
      filepath (1, :) char {mustBeFolder}
      datafile (1, :) char {mustBeMember(datafile, {'ice1', 'ice2', 'met'})}
      gridcell (1, 1) double {mustBeNumeric}
   end

   tmp = load(fullfile(filepath, [datafile '_' num2str(gridcell) '.mat'])).(datafile);

   switch datafile
      case 'ice1'

         data = tmp.Tsfc.'; % transpose to numlayers x numtimes
         varnames = tmp.Properties.VariableNames;

      case 'ice2'
         data = tmp.Tice;
         varnames = fieldnames(tmp);
   end

   % Get the dimensions of the data
   varsize = size(data); % numlayers x numtimes

end
