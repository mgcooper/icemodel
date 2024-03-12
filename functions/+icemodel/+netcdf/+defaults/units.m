function out = units(whichdata)

   % Define the grid and time dimension units
   dimunits = {
      'degrees_north', ... lat
      'degrees_east', ...  lon
      'm', ...             x
      'm', ...             y
      'm', ...             elev
      'm', ...             depth
      'time', ...          placeholder, replaced with 'seconds since YYYY-01-01 00:00:00' in code
      };

   % Define the variable units
   switch whichdata

      case 'ice1'

         out = {
            'degC', ...    Tsfc (changed from K due to conversion in POSTPROC)
            'm w.e.', ...  melt
            'm w.e.', ...  freeze
            'm w.e.', ...  subl
            'm w.e.', ...  cond
            'm w.e.', ...  runoff
            'm w.e.', ...  depth_melt
            'm w.e.', ...  depth_freeze
            'm w.e.', ...  surf_runoff
            'm w.e.', ...  column_runoff
            };

      case 'ice2'

         out = {
            'degC', ... % changed from K due to conversion in POSTPROC
            '1', ...
            '1', ...
            '1', ...
            };

      case 'met'

      case {'dimensions', 'dims'}

         out = dimunits;

      otherwise
         error('unrecognized icemodel data file name')
   end

   % Return as a column
   out = out(:);
end
