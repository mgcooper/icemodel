function dimsize = getdimsize(dimdata, dimnames)
   %GETDIMSIZE
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
      dimdata
      dimnames (1, :) string {mustBeText} = string(fieldnames(dimdata))
   end

   % The max(1, ...) prevents dimsize.depth = 0, for the ice1 case where
   % dimdata.depth is an empty array because Z = 0, dz = 0. If it is useful to
   % allow dimsize.depth = 0, e.g., for data which truly has 1 depth, remove
   % the max condition and update defdimid and anywhere dimsize.depth == 1
   % is used to identify ice1 data.

   for f = dimnames(:)'
      dimsize.(f) = max(1, numel(dimdata.(f)));
   end

end
