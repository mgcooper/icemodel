function met = completeMetVariables(met, kwargs)
   %COMPLETEMETVARIABLES Add missing met-contract variables as NaN placeholders.
   %
   %  met = icemodel.forcing.helpers.completeMetVariables(met)
   %  met = ... completeMetVariables(met, include_split_precip=true)
   %
   % Role
   %  Source builders call this at the met-building boundary. A source can be a
   %  useful native forcing record and still lack one or more required channels.
   %  A NaN means the source does not supply that channel and the runtime can
   %  substitute a value. A NaN does not mean zero flux.

   arguments
      met timetable
      kwargs.include_split_precip (1, 1) logical = false
   end

   required = icemodel.forcing.helpers.metvariables();
   if kwargs.include_split_precip
      required = [required, "rainf", "snowf"];
   end

   % Add the absent channels after the existing columns. data2met then reorders
   % the required variables. Direct source builders can keep the source order.
   varnames = string(met.Properties.VariableNames);
   missing = setdiff(required, varnames, 'stable');
   for varname = missing
      met.(varname) = nan(height(met), 1);
   end
end
