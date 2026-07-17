function tt = stampMetadata(tt, kwargs)
   %STAMPMETADATA Embed canonical CF-ish metadata in a timetable's properties.
   %
   %  tt = icemodel.forcing.helpers.stampMetadata(tt)
   %  tt = icemodel.forcing.helpers.stampMetadata(tt, strict=false)
   %
   % Stamps each variable of the table or timetable TT with its canonical
   % metadata from the single source icemodel.netcdf.defaults.variable, so
   % met/Data/observation files are self-describing:
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
   % An unmapped channel errors via the canonical map by default: every shipped
   % forcing column must be labelled. Verification tables can pass strict=false
   % so non-science string keys keep blank units/descriptions while known
   % science variables are still labelled.
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
