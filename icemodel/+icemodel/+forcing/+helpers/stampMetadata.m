function tt = stampMetadata(tt)
   %STAMPMETADATA Embed canonical CF-ish metadata in a timetable's properties.
   %
   %  tt = icemodel.forcing.helpers.stampMetadata(tt)
   %
   % Stamps each variable of the timetable TT with its canonical metadata
   % from the single source icemodel.netcdf.defaults.variable, so met/Data
   % files are self-describing:
   %
   %    Properties.VariableUnits        <- unit
   %    Properties.VariableDescriptions <- long_name
   %    Properties.CustomProperties.StandardNames <- CF standard_name
   %
   % Timetables have native slots for units and descriptions but none for CF
   % standard names, so the standard_name strings are carried in a
   % table-level CustomProperty StandardNames (a string array aligned to the
   % variable order). A channel with no CF name carries "" in that slot.
   %
   % An unmapped channel errors via the canonical map: every shipped column
   % must be labelled.
   %
   % See also: icemodel.netcdf.defaults.variable,
   %  icemodel.forcing.helpers.variableUnits, icemodel.forcing.data2met

   arguments
      tt timetable
   end

   names = string(tt.Properties.VariableNames);
   info = icemodel.netcdf.defaults.variable(names);

   tt.Properties.VariableUnits = {info.unit};
   tt.Properties.VariableDescriptions = {info.long_name};

   if ~isprop(tt.Properties.CustomProperties, "StandardNames")
      tt = addprop(tt, "StandardNames", "table");
   end
   tt.Properties.CustomProperties.StandardNames = string({info.standard_name});
end
