function out = varnames(whichdata)

   switch whichdata

      case 'ice1'

      case 'ice2'

         out = {
            'f_ice','f_liq','df_liq','Tice',...
            'latitude','longitude','x_easting','y_northing','elevation','time'};

      case 'met'
         
      case 'grid'
         

      otherwise
         error('unrecognized icemodel data file name')
   end
   
   % For all file types, append the grid and time dimensions
   'latitude','longitude',   ...
         'x_easting','y_northing','elevation','time'};
end
