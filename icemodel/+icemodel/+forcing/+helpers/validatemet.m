function validatemet(met)
   %VALIDATEMET Assert that MET satisfies the icemodel met-file contract.
   %
   %  icemodel.forcing.helpers.validatemet(met)
   %
   % Errors unless MET is a timetable with a regular time axis containing every
   % required met variable (see icemodel.forcing.helpers.metvariables). All-NaN
   % required variables are allowed as explicit placeholders for runtime
   % substitution, but absent variables remain invalid.
   %
   % See also: icemodel.forcing.helpers.metvariables,
   %  icemodel.forcing.helpers.writemet, icemodel.loadmet

   arguments
      met timetable
   end

   required = icemodel.forcing.helpers.metvariables();
   varnames = string(met.Properties.VariableNames);

   missing = setdiff(required, varnames);
   if ~isempty(missing)
      error('icemodel:forcing:validatemet:missingVariables', ...
         'met file is missing required variable(s): %s', ...
         strjoin(missing, ', '));
   end

   if height(met) < 2
      error('icemodel:forcing:validatemet:tooFewSamples', ...
         'met file must contain at least two samples');
   end

   steps = diff(met.Time);
   if any(steps ~= steps(1))
      error('icemodel:forcing:validatemet:irregularTimeAxis', ...
         'met file time axis must have a uniform timestep');
   end

   % Precipitation-rate unit. When the met timetable records VariableUnits,
   % the ppt channel must carry the canonical water-equivalent rate (m s-1;
   % see icemodel.forcing.helpers.metvariables) so every source agrees. A met
   % file with no VariableUnits is accepted (legacy artifacts predate the
   % metadata), but a ppt unit that is set and wrong is rejected.
   units = string(met.Properties.VariableUnits);
   if ~isempty(units)
      [~, ~, pptunit] = icemodel.forcing.helpers.metvariables();
      pptidx = varnames == "ppt";
      if any(pptidx) && strlength(units(pptidx)) > 0 ...
            && units(pptidx) ~= pptunit
         error('icemodel:forcing:validatemet:pptUnit', ...
            'ppt unit must be the canonical "%s", got "%s"', ...
            pptunit, units(pptidx));
      end
   end
end
