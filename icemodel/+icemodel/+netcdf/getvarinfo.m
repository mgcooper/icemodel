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

   % varsize is meant to catch the case where ice2 files have different Z
   % within the same year. This function sets the size from one loaded test
   % file, so it does not catch a size change from file to file.
   %
   % getchunksize uses the value. defdimid must use it too, or get it
   % another way.
   %
   % If the code writes each ice2 file to its own nc file, writeice2 must
   % update this. If the code writes multiple ice2 files to one nc file,
   % use the supplied opts.dz/Z, and possibly remove NOFILL, to allow
   % different array sizes.

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
