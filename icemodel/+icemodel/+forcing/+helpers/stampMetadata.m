function tt = stampMetadata(tt, kwargs)
   %STAMPMETADATA Embed canonical CF-ish metadata in a timetable's properties.
   %
   %  tt = icemodel.forcing.helpers.stampMetadata(tt)
   %  tt = icemodel.forcing.helpers.stampMetadata(tt, strict=false)
   %
   % Stamps each variable of the table or timetable TT with its canonical
   % metadata from icemodel.netcdf.defaults.variable, so met, Data, and
   % observation files describe themselves:
   %
   %    Properties.VariableUnits        <- unit
   %    Properties.VariableDescriptions <- long_name
   %    Properties.CustomProperties.StandardNames <- CF standard_name
   %
   % Timetables have native slots for units and descriptions, but none for CF
   % standard names. This function therefore stores the standard_name strings
   % in the table-level CustomProperty StandardNames. That property is a string
   % array in the same order as the variables. A channel with no CF name
   % carries "" in its slot.
   %
   % By default, an unmapped channel raises the error from the canonical map,
   % because every shipped forcing column must carry a label. A verification
   % table can pass strict=false. Then a non-science string key keeps blank
   % units and descriptions, and a known science variable still gets its label.
   %
   % See also: icemodel.netcdf.defaults.variable,
   %  icemodel.forcing.helpers.variableUnits, icemodel.forcing.data2met

   arguments
      tt
      kwargs.strict (1, 1) logical = true
   end

   if ~istable(tt) && ~istimetable(tt)
      error('icemodel:forcing:stampMetadata:badInput', ...
         'stampMetadata requires a table or timetable')
   end

   % All icemodel timetables expose their row times as TT.Time. Stamping is a
   % common write-path step, so normalize the row-dimension name here.
   if istimetable(tt)
      tt.Properties.DimensionNames{1} = 'Time';
   end

   names = string(tt.Properties.VariableNames);
   info = metadataFor(names, kwargs.strict);

   tt.Properties.VariableUnits = {info.unit};
   tt.Properties.VariableDescriptions = {info.long_name};

   if ~isprop(tt.Properties.CustomProperties, "StandardNames")
      tt = addprop(tt, "StandardNames", "table");
   end
   tt.Properties.CustomProperties.StandardNames = string({info.standard_name});
end

function info = metadataFor(names, strict)
   %METADATAFOR Return canonical metadata, optionally blanking unknown columns.

   info = repmat(emptyInfo(), 1, numel(names));
   for k = 1:numel(names)
      try
         info(k) = icemodel.netcdf.defaults.variable(names(k));
      catch err
         if strict || ~strcmp(err.identifier, ...
               'icemodel:netcdf:variable:unknownChannel')
            rethrow(err)
         end
      end
   end
end

function info = emptyInfo()
   %EMPTYINFO Prototype for unknown non-science verification columns.
   info = struct('standard_name', '', 'long_name', '', 'unit', '', ...
      'is_cf', false);
end
