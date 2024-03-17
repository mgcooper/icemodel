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

   % Nominally varsize is designed to catch the case where ice2 files have
   % different Z within the same year. However, since datasize is set here by
   % loading one test file, this won't actually catch the case where the size
   % changes from file to file.
   %
   % The value is used in getchunksize, but needs to also be used in defdimid,
   % or a method to
   %
   % If each individual ice2 file is written to an nc file, this would
   % need to be updated in writeice2. If multiple ice2 files are written
   % to one nc file, then the supplied opts.dz/Z should be used and
   % possibly NOFILL removed to account for different sized arrays.

   tmp = load(fullfile( ...
      filepath, [datafile '_' num2str(gridcell) '.mat'])).(datafile);

   switch datafile
      case 'ice1'

         data = tmp.Tsfc.'; % transpose to numlayers x numtimes
         varnames = tmp.Properties.VariableNames';

      case 'ice2'
         data = tmp.Tice;
         varnames = fieldnames(tmp)';
   end

   % Get the dimensions of the data
   varsize = size(data); % numlayers x numtimes
end
