function out = axes(whichdata)

   % Define the grid and time dimension units. Use empty char '' for variables
   % which are not technically "dimensions".
   dimaxes = {
      '', ...     gridcell
      '', ...     latitude  - I removed Y based on Example 5.2
      '', ...     longitude - I removed X based on Example 5.2
      'X', ...    x_easting
      'Y', ...    y_northing
      '', ...     elevation (not a dimension)
      'Z', ...    depth
      'T', ...    time
      };
   % Regarding lat lon, specify them as 'coordinate' attributes
   % Also see Example 5.17

   % Only the dimensions have the 'axis' attribute, the blanks below are dummy
   % values so icemodel.netcdf.getdefaults receives an output from this function
   % when called with 'ice1', 'ice2', etc.
   switch whichdata

      case 'ice1'

         out = {''};

      case 'ice2'

         out = {''};

      case 'met'

         out = {''};

      case {'dimensions', 'dims'}

         out = dimaxes;

      otherwise
         error('unrecognized icemodel data file name')
   end

   % Return as a column
   out = out(:);
end
