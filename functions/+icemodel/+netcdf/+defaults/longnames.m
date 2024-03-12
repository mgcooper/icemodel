function out = longnames(whichdata)

   % Define the grid and time dimension names
   dimnames = {
      'latitude of grid cell center', ...
      'longitude of grid cell center', ...
      'x-coordinate of grid cell center in EPSG:3413 WGS 84 / NSIDC Sea Ice Polar Stereographic North', ...
      'y-coordinate of grid cell center in EPSG:3413 WGS 84 / NSIDC Sea Ice Polar Stereographic North', ...
      'elevation above mean sea level', ...
      'depth below the surface', ...
      'seconds since 1 January, 00:00:00 (UTC)'
      };

   % Define the variable names
   switch whichdata

      case 'ice1'

         out = {
            'thermodynamic temperature of ice surface at end of timestep', ...            Tsfc
            'cumulative melt amount in meters liquid water equivalent', ...               melt
            'cumulative refreeze amount in meters liquid water equivalent', ...           freeze
            'cumulative sublimation amount in meters liquid water equivalent', ...        subl
            'cumulative condensation amount in meters liquid water equivalent', ...       cond
            'cumulative runoff amount in meters liquid water equivalent', ...             runoff
            'cumulative melt amount in meters liquid water equivalent', ...               depth_melt
            'cumulative refreeze amount in meters liquid water equivalent', ...           depth_freeze
            'cumulative surface runoff amount in meters liquid water equivalent', ...     surf_runoff
            'cumulative subsurface runoff amount in meters liquid water equivalent', ...  column_runoff
            };

      case 'ice2'

         out = {
            'thermodynamic temperature of control volume (layer) at end of timestep', ...
            'total change in volume fraction of unfrozen water in control volume (layer) during timestep inclusive of drained water',...
            'volume fraction of frozen water in control volume (layer) at end of timestep', ...
            'volume fraction of unfrozen water in control volume (layer) at end of timestep', ...
            };

      case 'met'

      case {'dimensions', 'dims'}

         out = dimnames;

      otherwise
         error('unrecognized icemodel data file name')
   end

   % Return as a column
   out = out(:);
end
