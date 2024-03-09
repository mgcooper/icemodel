function out = longnames(whichdata)

   switch whichdata

      case 'ice1'
         
         

      case 'ice2'

         out = {
            'fraction of frozen water in control volume', ...
            'fraction of unfrozen water in control volume', ...
            'change in fraction of unfrozen water in control volume',...
            'thermodynamic temperature of control volume'};

      case 'met'

      otherwise

         error('unrecognized icemodel data file name')
   end
end
