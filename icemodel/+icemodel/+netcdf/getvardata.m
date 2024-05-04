function data = getvardata(filepath, varnames, dimdata, xtype, ...
      startcell, countcell, smbmodel)
   %GETVARDATA Read icemodel data into memory and fill the arrays.
   %
   % Allocate arrays in column-major format for efficiency.
   % Permute them to match the defined dimensions when writing to netcdf.
   %
   % For now, only allow chunking over cells. Note this is not netcdf chunking,
   % this preallocates and reads the ice2 arrays into memory in chunks.
   %
   % See also:

   arguments
      filepath (1, :) char {mustBeFolder}
      varnames (1, :) cell {mustBeText}
      dimdata (1, :) struct {mustBeStruct}
      xtype (1, :) char {mustBeText}
      startcell (1, 1) double {mustBeNumeric} = 1
      countcell (1, 1) double {mustBeNumeric} = NaN
      smbmodel (1, :) char {mustBeText} = ""
   end

   % Get the dimension sizes
   dimsize = icemodel.netcdf.getdimsize(dimdata);

   % Check if this is ice1 or ice2
   if dimsize.depth == 1
      sdata = 'ice1';
   else
      sdata = 'ice2';
   end

   % Short circuit if this is an ice2 file. If an entire year of ice2
   % data can fit in memory, this function will allocate it.
   if strcmp(sdata, 'ice2') && isnan(countcell)
      data = [];
      warning('For ICE2, set COUNTCELL. Returning empty data.')
      return
   end

   % Define the dimension sizes for storing the data in memory.
   numtime = dimsize.time;
   numlyrs = dimsize.depth;

   % Set numcells, the number of grid cells processed per write.
   if isnan(countcell)
      gridcells = dimdata.gridcell; % all cells
   else
      gridcells = startcell:(startcell + countcell - 1);
   end
   numcell = numel(gridcells);

   % Convert netcdf type to matlab type.
   mtype = icemodel.netcdf.nctype2mat(xtype);
   if ismember(mtype, {'char', 'string', 'cell'})
      error( ...
         'Preallocation for %s is not supported in this function.', xtype);
   end

   % Preallocate the data arrays. Squeeze to remove singleton dimensions.
   for v = 1:numel(varnames)
      thisvar = varnames{v};
      data.(thisvar) = squeeze(nan(numlyrs, numtime, numcell, mtype));
   end

   % Fill the data arrays from gridcell(1) to gridcell(numcell)
   for n = 1:numcell

      thiscell = num2str(gridcells(n));
      tmp = load(fullfile(filepath, [sdata '_' thiscell '.mat'])).(sdata);

      % Set freeze zero for skinmodel (freeze is Qf, not refreeze)
      if strcmp(sdata, 'ice1') && strcmp('skinmodel', smbmodel) && isvariable('freeze', tmp)
         tmp.freeze = 0 * tmp.freeze;
      end

      for v = 1:numel(varnames)
         thisvar = varnames{v};
         shape = size(tmp.(thisvar)); % nhrs x 1 (ice1), nlyrs x nhrs (ice2)

         if strcmp(sdata, 'ice1')
            data.(thisvar)(1:shape(1), n) = tmp.(thisvar); % nhrs x ncells

         elseif strcmp(sdata, 'ice2')
            data.(thisvar)(1:shape(1), :, n) = tmp.(thisvar); % nlyrs x nhrs x ncells
         end
      end
   end
end
