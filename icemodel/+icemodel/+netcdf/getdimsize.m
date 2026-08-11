function dimsize = getdimsize(dimdata, dimnames)
   %GETDIMSIZE Return the size of each named dimension.
   %
   %  DIMSIZE = GETDIMSIZE(DIMDATA)
   %  DIMSIZE = GETDIMSIZE(DIMDATA, DIMNAMES)
   %
   % For example:
   %
   % dimsize.gridcell = numel(dimdata.gridcell);  % number of grid cells
   % dimsize.depth = numel(dimdata.depth);        % number of vertical layers
   % dimsize.time = numel(dimdata.time);          % number of timesteps per year
   %
   % See also:

   arguments
      dimdata (1, :) struct {mustBeStruct}
      dimnames (1, :) string {mustBeText} = string(fieldnames(dimdata))
   end

   % The max(1, ...) prevents dimsize.depth = 0 in the ice1 case, where
   % dimdata.depth is an empty array because Z = 0 and dz = 0. To allow
   % dimsize.depth = 0, for example for data that truly has 1 depth, remove
   % the max condition. Then update defdimid and every place that uses
   % dimsize.depth == 1 to identify ice1 data.

   for f = dimnames(:)'
      dimsize.(f) = max(1, numel(dimdata.(f)));
   end

end
