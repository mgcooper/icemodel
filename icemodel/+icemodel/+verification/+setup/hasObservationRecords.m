function tf = hasObservationRecords(observations)
   %HASOBSERVATIONRECORDS True when any observation sub-bundle carries rows.
   %
   %  tf = icemodel.verification.setup.hasObservationRecords(observations)
   %
   % SUMup-derived observation bundles are nested structs of tables/timetables.
   % Importers use this helper to convert all-empty selections into skippable
   % missing-observation cases instead of writing empty observations.mat files.
   tf = false;
   if ~isstruct(observations)
      tf = hasRows(observations);
      return
   end

   fields = setdiff(string(fieldnames(observations)), "format", 'stable');
   for field = reshape(fields, 1, [])
      if hasRows(observations.(char(field)))
         tf = true;
         return
      end
   end
end

function tf = hasRows(value)
   %HASROWS Return true for non-empty observation tables or nested payloads.
   if istable(value) || istimetable(value)
      tf = height(value) > 0;
   elseif isstruct(value)
      tf = icemodel.verification.setup.hasObservationRecords(value);
   else
      tf = ~isempty(value);
   end
end
