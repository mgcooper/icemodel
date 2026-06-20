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
   % Precipitation-rate harmonization: the gridded-source Data builders carry
   % their precipitation channels in mWE/h (metres water equivalent per hour),
   % but the canonical met-file precipitation unit is m s-1 (a water-equivalent
   % RATE; see icemodel.forcing.helpers.metvariables), the unit
   % icemodel.surface.advective_heat_flux consumes and the one ESM-SnowMIP met
   % already uses. So the precipitation channels (ppt and any rain/snow split)
   % are converted from the source unit to m s-1 here, and the result's
   % Properties.VariableUnits records m s-1 for them. The conversion is driven
   % by PPTSOURCEUNIT so a caller can pass already-rate Data through unchanged.
   %
   % Inputs
   %  Data - timetable from a build<Source>Data builder (or a legacy
   %         userdata Data file)
   %
   % Name-value
   %  validate      : assert the met contract on the result (default true)
   %  pptsourceunit : unit of the incoming precipitation channels, one of
   %      "mWE/h" (default, the gridded-builder convention) or "m s-1"
   %      (already a rate; no conversion). Converted to the canonical
   %      m s-1 met unit.
   %
   % Outputs
   %  met - timetable ready for icemodel.forcing.helpers.writemet
   %
   % Legacy: reimplements runoff/functions/marData2Met.m (retained,
   % unchanged, as the legacy reference). Unlike the legacy version it derives
   % ppt from the rain/snow split, drops the bookkeeping `date` column (the
   % output is a timetable with a Time axis), orders required variables
   % first, and harmonizes the precipitation rate to the canonical m s-1.
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.forcing.helpers.writemet,
   %  icemodel.forcing.helpers.metvariables

   arguments
      Data timetable
      kwargs.validate (1, 1) logical = true
      kwargs.pptsourceunit (1, 1) string {mustBeMember( ...
         kwargs.pptsourceunit, ["mWE/h", "m s-1"])} = "mWE/h"
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

   % Harmonize the precipitation rate to the canonical m s-1. mWE/h is an
   % hourly water-equivalent depth rate, so dividing by 3600 s/h yields m s-1;
   % "m s-1" input is already canonical and passes through. Applied to every
   % precipitation channel present (ppt and any rain/snow split) so the
   % derived ppt = rain + snow identity is preserved (both sides scale
   % identically). Non-precipitation rates (melt/runoff/smb) keep mWE/h.
   if kwargs.pptsourceunit == "mWE/h"
      pptchannels = ["ppt", "rain", "snow", "rainf", "snowf"];
      for ch = pptchannels(ismember(pptchannels, ...
            string(met.Properties.VariableNames)))
         met.(ch) = met.(ch) / 3600;
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

   % Record the canonical precipitation-rate unit on the met output so the
   % m s-1 convention is self-describing on disk. Preserve any units the
   % source Data carried for the surviving columns; stamp m s-1 on the
   % precipitation channels and leave the rest unset where unknown.
   met = labelPrecipUnits(met, Data);

   if kwargs.validate
      icemodel.forcing.helpers.validatemet(met)
   end
end

%% Local functions
function met = labelPrecipUnits(met, Data)
   %LABELPRECIPUNITS Set met VariableUnits, stamping m s-1 on precip channels.
   %
   % Carries over the source Data's VariableUnits for columns that survive
   % into the met output (matched by name, so reordering is handled), then
   % overwrites the precipitation channels with the canonical m s-1 rate. A
   % column with no known source unit is left as an empty string.
   names = string(met.Properties.VariableNames);
   units = repmat("", 1, numel(names));

   % Inherit source units by name where the Data carried them.
   srcnames = string(Data.Properties.VariableNames);
   srcunits = string(Data.Properties.VariableUnits);
   if ~isempty(srcunits)
      [tf, loc] = ismember(names, srcnames);
      units(tf) = srcunits(loc(tf));
   end

   % Canonical precipitation rate.
   [~, ~, pptunit] = icemodel.forcing.helpers.metvariables();
   pptchannels = ["ppt", "rain", "snow", "rainf", "snowf"];
   units(ismember(names, pptchannels)) = pptunit;

   met.Properties.VariableUnits = cellstr(units);
end
