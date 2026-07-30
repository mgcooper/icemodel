function units = variableUnits(names)
   %VARIABLEUNITS Canonical unit string for each forcing-builder channel.
   %
   %  units = icemodel.forcing.helpers.variableUnits(names)
   %
   % Thin wrapper over the single canonical variable-metadata source
   % icemodel.netcdf.defaults.variable: returns the unit field for each
   % channel in NAMES. There is NO duplicated unit list here anymore; units,
   % CF standard_names, and long_names all come from one map
   % (icemodel.netcdf.defaults.variables) so every forcing source (MAR,
   % MERRA-2, RACMO, PROMICE) and the netcdf writers agree.
   %
   % NAMES is a string array of channel names; UNITS is a cellstr of the
   % matching unit strings, ready to assign to a timetable's
   % Properties.VariableUnits. The indexed ice-temperature string (ticeN, K)
   % and thermistor-depth string (dticeN, m) channels resolve by pattern.
   %
   % An emitted channel missing from the canonical map is an error: the
   % builders must never ship an unlabeled column. Add the channel to
   % icemodel.netcdf.defaults.variables when a builder gains a new output.
   %
   % See also: icemodel.netcdf.defaults.variable,
   %  icemodel.netcdf.defaults.variables,
   %  icemodel.forcing.helpers.metvariables, icemodel.forcing.data2met

   arguments
      names (1, :) string
   end

   units = strings(1, numel(names));
   for k = 1:numel(names)
      try
         info = icemodel.netcdf.defaults.variable(names(k));
      catch err
         if strcmp(err.identifier, 'icemodel:netcdf:variable:unknownChannel')
            error('icemodel:forcing:variableUnits:unmappedChannel', ...
               ['no canonical unit for channel "%s"; add it to ' ...
               'icemodel.netcdf.defaults.variables'], names(k));
         end
         rethrow(err)
      end
      units(k) = string(info.unit);
   end
   units = cellstr(units);
end
