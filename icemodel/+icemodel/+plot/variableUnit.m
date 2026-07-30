function unit = variableUnit(T, varname)
   %VARIABLEUNIT Return source table units with canonical metadata fallback.
   unit = "";
   names = string(T.Properties.VariableNames);
   idx = find(names == varname, 1);
   if ~isempty(idx) && numel(T.Properties.VariableUnits) >= idx
      unit = string(T.Properties.VariableUnits{idx});
   end
   if unit ~= ""
      return
   end

   % Unknown variables can still be plotted; they simply receive no unit suffix.
   try
      info = icemodel.netcdf.defaults.variable(varname);
      unit = string(info.unit);
   catch
      unit = "";
   end
end
