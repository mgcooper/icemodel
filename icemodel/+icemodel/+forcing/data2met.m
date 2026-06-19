function met = data2met(Data, kwargs)
   %DATA2MET Convert a Data timetable to an icemodel met timetable.
   %
   %  met = icemodel.forcing.data2met(Data)
   %  met = ... data2met(Data, validate=false)
   %
   % Generalizes the legacy marData2Met: selects the met-contract
   % variables from a Data timetable (any source), derives the total
   % precipitation channel when the source carries a rain/snow split,
   % and orders required variables first. Variables outside the met
   % contract pass through after the contract set; bookkeeping columns
   % (date) are dropped.
   %
   %    ppt = snow + rain (or rainf + snowf), when not already present
   %
   % Inputs
   %  Data - timetable from a build<Source>Data builder (or a legacy
   %         userdata Data file)
   %
   % Name-value
   %  validate : assert the met contract on the result (default true)
   %
   % Outputs
   %  met - timetable ready for icemodel.forcing.helpers.writemet
   %
   % Legacy: reimplements runoff/functions/marData2Met.m (retained,
   % unchanged, as the legacy reference). Unlike the legacy version it derives
   % ppt from the rain/snow split, drops the bookkeeping `date` column (the
   % output is a timetable with a Time axis), and orders required variables
   % first.
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.forcing.helpers.writemet,
   %  icemodel.forcing.helpers.metvariables

   arguments
      Data timetable
      kwargs.validate (1, 1) logical = true
   end

   met = Data;
   varnames = string(met.Properties.VariableNames);

   % Total precipitation from the source's rain/snow split.
   if ~ismember("ppt", varnames)
      if all(ismember(["rain", "snow"], varnames))
         met.ppt = met.rain + met.snow;
      elseif all(ismember(["rainf", "snowf"], varnames))
         met.ppt = met.rainf + met.snowf;
      end
   end

   % Drop bookkeeping columns.
   met = removevars(met, intersect("date", ...
      string(met.Properties.VariableNames)));

   % Required contract variables first, everything else after.
   required = icemodel.forcing.helpers.metvariables();
   varnames = string(met.Properties.VariableNames);
   ordered = [required(ismember(required, varnames)), ...
      varnames(~ismember(varnames, required))];
   met = met(:, cellstr(ordered));

   if kwargs.validate
      icemodel.forcing.helpers.validatemet(met)
   end
end
