function data = getvardata(filepath, vars, dims, xtype, simmodel)
   %GETVARDATA Read icemodel data and fill the arrays.
   %
   % Allocate arrays in column-major format for efficiency.
   % Permute them to match the defined dimensions when writing to netcdf.
   %
   % See also:

   arguments
      filepath (1, :) char {mustBeFolder}
      vars
      dims
      xtype
      simmodel
   end

   dimsize = icemodel.netcdf.getdimsize(dims);
   numcell = dimsize.gridcell;
   numtime = dimsize.time;
   numlyrs = dimsize.depth;

   % Check if this is ice1 or ice2
   if dimsize.depth == 1
      sdata = 'ice1';
   else
      sdata = 'ice2';
   end

   % Short circuit if this is an ice2 file. If an entire year of ice2
   % data can fit in memory, this function will allocate it.
   if strcmp(sdata, 'ice2')
      data = [];
      return
   end

   % Convert netcdf type to matlab type
   mtype = nctype2mat(xtype);
   if ismember(mtype, {'char', 'string', 'cell'})
      error( ...
         'Preallocation for %s is not supported in this function.', xtype);
   end

   % Preallocate the data arrays using NaNs. Squeeze to remove singleton dimensions.
   for v = 1:numel(vars)
      thisvar = vars{v};
      data.(thisvar) = squeeze(nan(numlyrs, numtime, numcell, mtype));
   end

   % Fill the data arrays
   for n = 1:numcell
      tmp = load(fullfile(filepath, [sdata '_' num2str(n) '.mat'])).(sdata);

      % Set freeze zero for skinmodel (freeze is Qf, not refreeze)
      if strcmp('skinmodel', simmodel) && isvariable('freeze', tmp)
         tmp.freeze = 0 * tmp.freeze;
      end

      for v = 1:numel(vars)
         thisvar = vars{v};
         shape = size(tmp.(thisvar)); % nhrs x 1 (ice1), nlyrs x nhrs (ice2)

         if strcmp(sdata, 'ice1')
            data.(thisvar)(1:shape(1), n) = tmp.(thisvar); % nhrs x ncells

         elseif strcmp(sdata, 'ice2')
            data.(thisvar)(1:shape(1), :, n) = tmp.(thisvar); % nlyrs x nhrs x ncells
         end
      end
   end
end
