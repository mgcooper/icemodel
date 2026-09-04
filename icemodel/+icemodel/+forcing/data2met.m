function met = data2met(Data, kwargs)
   %DATA2MET Convert a Data timetable to an icemodel met timetable.
   %
   %  met = icemodel.forcing.data2met(Data)
   %  met = ... data2met(Data, validate=false, fillwithmissing=true)
   %
   % Selects the required met variables from a Data timetable. It derives total
   % precipitation from a rain/snow split, removes `date`, and places required
   % variables first. Other variables follow the required variables.
   %
   %    ppt = snow + rain (or rainf + snowf), when not already present
   %
   % Precipitation is a water-equivalent rate in m s-1 (see
   % icemodel.forcing.helpers.metvariables). icemodel.surface.advective_heat_flux
   % takes that unit, and the ESM-SnowMIP met uses it. The derived total
   % ppt = rain + snow therefore has the unit m s-1. The result's
   % Properties.VariableUnits comes from the shared unit map.
   %
   % Inputs
   %  Data - timetable from a build<Source>Data builder (or a legacy
   %         userdata Data file)
   %
   % Name-value
   %  validate : run icemodel.forcing.helpers.validatemet on the result
   %             (default true)
   %  fillwithmissing : add absent required met channels as NaN placeholders
   %                    before validation (default true)
   %
   % Outputs
   %  met - timetable ready for icemodel.forcing.helpers.writemet
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.forcing.helpers.writemet,
   %  icemodel.forcing.helpers.metvariables

   arguments
      Data timetable
      kwargs.validate (1, 1) logical = true
      kwargs.fillwithmissing (1, 1) logical = true
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

   % Required met variables first, everything else after.
   required = icemodel.forcing.helpers.metvariables();
   varnames = string(met.Properties.VariableNames);
   ordered = [required(ismember(required, varnames)), ...
      varnames(~ismember(varnames, required))];
   met = met(:, cellstr(ordered));

   % Label each channel with units, a long name, and a CF standard name.
   % Precipitation, including derived ppt, uses m s-1.
   met = icemodel.forcing.helpers.stampMetadata(met);

   if kwargs.validate
      icemodel.forcing.helpers.validatemet(met)
   end
end
