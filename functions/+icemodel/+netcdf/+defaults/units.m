function out = units(whichdata)

   switch whichdata

      case 'ice1'

      case 'ice2'

         out = {'1','1','1','K'};

      case 'met'

      otherwise
         error('unrecognized icemodel data file name')
   end
end
