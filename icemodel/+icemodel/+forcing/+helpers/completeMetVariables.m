function met = completeMetVariables(met, kwargs)
   %COMPLETEMETVARIABLES Add missing met-contract variables as NaN placeholders.
   %
   %  met = icemodel.forcing.helpers.completeMetVariables(met)
   %  met = ... completeMetVariables(met, include_split_precip=true)
   %
   % Role
   %  Source builders use this at the met-building boundary when a source is a
   %  useful native forcing record but lacks one or more required channels. The
   %  NaNs mean "missing by source design and available for runtime substitution",
   %  not zero flux.

   arguments
      met timetable
      kwargs.include_split_precip (1, 1) logical = false
   end

   required = icemodel.forcing.helpers.metvariables();
   if kwargs.include_split_precip
      required = [required, "rainf", "snowf"];
   end

   % Add absent channels after existing columns; data2met reorders required
   % variables afterwards, while direct source builders may keep source order.
   varnames = string(met.Properties.VariableNames);
   missing = setdiff(required, varnames, 'stable');
   for varname = missing
      met.(varname) = nan(height(met), 1);
   end
end
