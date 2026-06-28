function met = data2met(Data, kwargs)
   %DATA2MET Convert a Data timetable to an icemodel met timetable.
   %
   %  met = icemodel.forcing.data2met(Data)
   %  met = ... data2met(Data, validate=false, fillwithmissing=true)
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
   % Precipitation rate: the gridded-source Data builders already emit their
   % precipitation channels in the canonical m s-1 water-equivalent rate (see
   % icemodel.forcing.helpers.metvariables), the unit
   % icemodel.surface.advective_heat_flux consumes and the one ESM-SnowMIP met
   % already uses. So the derived total ppt = rain + snow inherits m s-1, and
   % the result's Properties.VariableUnits comes from the shared canonical
   % unit map.
   %
   % Inputs
   %  Data - timetable from a build<Source>Data builder (or a legacy
   %         userdata Data file)
   %
   % Name-value
   %  validate : assert the met contract on the result (default true)
   %  fillwithmissing : add absent required met channels as NaN placeholders
   %                    before validation (default false)
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
      kwargs.fillwithmissing (1, 1) logical = false
   end

   met = Data;
   varnames = string(met.Properties.VariableNames);

   % Total precipitation from the source's rain/snow split (in m s-1, the unit
   % the split channels already carry).
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

   % Optional completion lets source adapters write native met products even
   % when a required channel must be supplied by runtime substitution later.
   if kwargs.fillwithmissing
      met = icemodel.forcing.helpers.completeMetVariables(met);
   end

   % Required contract variables first, everything else after.
   required = icemodel.forcing.helpers.metvariables();
   varnames = string(met.Properties.VariableNames);
   ordered = [required(ismember(required, varnames)), ...
      varnames(~ismember(varnames, required))];
   met = met(:, cellstr(ordered));

   % Self-describing metadata from the single canonical source: unit,
   % long_name (VariableDescriptions), and CF standard_name (StandardNames
   % custom property). Precipitation channels carry the canonical m s-1 rate
   % (including the derived ppt); every surviving channel is labelled.
   met = icemodel.forcing.helpers.stampMetadata(met);

   if kwargs.validate
      icemodel.forcing.helpers.validatemet(met)
   end
end
