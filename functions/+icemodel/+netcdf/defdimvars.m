function defdimvars(ncid, dimid, vars, standardnames, longnames, units, axes, ncprops)
   %DEFDIMVARS Define icemodel netcdf grid and time dims and attributes
   %
   %  DEFDIMVARS(NCID, DIMID, VARS, STANDARDNAMES, LONGNAMES, UNITS, AXES, NCPROPS)
   %
   % See also:

   arguments
      ncid (1, 1)
      dimid
      vars
      standardnames
      longnames
      units
      axes
      ncprops.shuffle (1, 1) logical = true
      ncprops.deflate (1, 1) logical = true
      ncprops.deflateLevel (1, 1) double = 1
   end

   for v = 1:numel(vars)
      thisvar = vars{v};

      % This check is only made on depth b/c it's the only dimension that
      % differs between ice1/ice2. If dimid had 1:1 mapping with vars, like the
      % other attrs, then this check could apply to all vars.
      if strcmp(thisvar, 'depth') && ~isfield(dimid, 'depth')
         continue
      end

      % Explicitly assign the fundamental dimensions. The "if isfield(...)"
      switch thisvar

         case "gridcell"

            thisvid = netcdf.defVar(ncid, thisvar, 'NC_INT', dimid.gridcell);

         case "depth"

            thisvid = netcdf.defVar(ncid, thisvar, 'NC_FLOAT', dimid.depth);

         case "time"
            thisvid = netcdf.defVar(ncid, thisvar, 'NC_DOUBLE', dimid.time);

         otherwise % lat, lon, x, y, elev
            thisvid = netcdf.defVar(ncid, thisvar, 'NC_FLOAT', dimid.gridcell);
      end

      netcdf.defVarDeflate(ncid, thisvid, ...
         ncprops.shuffle, ncprops.deflate, ncprops.deflateLevel);

      % standard names
      if ~isempty(standardnames{v})
         netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames{v});
      end

      % long names
      if ~isempty(longnames{v})
         netcdf.putAtt(ncid, thisvid, 'long_name', longnames{v});
      end

      % Units
      if ~isempty(units{v})
         netcdf.putAtt(ncid, thisvid, 'units', units{v});
      end

      % netcdf.putAtt(ncid, thisvid, '_FillValue', -9999);

      % Axis
      if ~isempty(axes{v})
         netcdf.putAtt(ncid, thisvid, 'axis', axes{v});
      end

      % Additional attributes e.g. coordinates (linked to auxiliary coordinate
      % variables), positive direction for vertical coordinates, calendar, etc.
      switch thisvar

         case "gridcell"
            netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');
            netcdf.putAtt(ncid, thisvid, 'comment', 'A unique identifier for each grid cell');

         case ["x_easting", "y_northing"]

            netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');

         case "elevation"
            netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');
            netcdf.putAtt(ncid, thisvid, 'positive', 'up');

         case "depth"
            netcdf.putAtt(ncid, thisvid, 'positive', 'down');

         case "time"
            netcdf.putAtt(ncid, thisvid, 'calendar', 'noleap');

         otherwise
            % lat, lon - do not add the coordinates
      end
   end
end
