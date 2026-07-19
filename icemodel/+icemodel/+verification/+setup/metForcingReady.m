function [tf, reason, complete_windows] = metForcingReady(met)
   %METFORCINGREADY Test unfilled readiness and inventory complete windows.
   %
   % Verification importers may stage met files with explicit all-NaN
   % placeholders so missing channels can be filled later. Such artifacts are
   % useful cache/data products, but they should not be advertised as runnable
   % forcing sources until every required met channel has finite data on one
   % regular time axis. COMPLETE_WINDOWS reports inclusive contiguous runs where
   % every required channel is finite; callers choose the scientifically useful
   % run length rather than inheriting a hidden minimum-duration threshold.

   arguments
      met timetable
   end

   required = icemodel.forcing.helpers.metvariables();
   names = string(met.Properties.VariableNames);
   missing = strings(0, 1);
   row_complete = true(height(met), 1);

   % Intersect per-row finite support while retaining channel names for the
   % actionable readiness reason.
   for v = reshape(required, 1, [])
      if ~ismember(v, names)
         missing(end + 1, 1) = v; %#ok<AGROW>
         row_complete(:) = false;
         continue
      end

      % A matrix-valued timetable variable is complete only when every value
      % associated with that posting is finite.
      channel_complete = all(isfinite(met.(char(v))), 2);
      row_complete = row_complete & channel_complete;
      if ~all(channel_complete)
         missing(end + 1, 1) = v; %#ok<AGROW>
      end
   end

   % Split finite runs at omitted, duplicate, reversed, or off-cadence row
   % times. A gap in the coordinate is as non-runnable as an explicit NaN row.
   adjacent = true(max(height(met) - 1, 0), 1);
   if height(met) > 1
      spacing = seconds(diff(met.Time));
      positive_spacing = spacing(isfinite(spacing) & spacing > 0);
      if isempty(positive_spacing)
         adjacent(:) = false;
      else
         cadence = median(positive_spacing);
         adjacent = spacing == cadence;
      end
   end

   if height(met) == 0
      start_index = zeros(0, 1);
      stop_index = zeros(0, 1);
   else
      starts = row_complete & [true; ~row_complete(1:end - 1) | ~adjacent];
      stops = row_complete & [~row_complete(2:end) | ~adjacent; true];
      start_index = find(starts);
      stop_index = find(stops);
   end
   complete_windows = table(met.Time(start_index), met.Time(stop_index), ...
      stop_index - start_index + 1, 'VariableNames', ...
      {'start_time', 'end_time', 'sample_count'});

   % Full-artifact readiness remains deliberately strict: 100% finite required
   % channels and no inferred coordinate gap over the declared artifact.
   tf = height(met) > 0 && all(row_complete) && all(adjacent);
   reasons = strings(0, 1);
   if height(met) == 0
      reasons(end + 1, 1) = "empty met timetable";
   end
   if ~isempty(missing)
      reasons(end + 1, 1) = "required met placeholder/gap channel(s): " ...
         + strjoin(unique(missing, 'stable'), ", ");
   end
   if any(~adjacent)
      reasons(end + 1, 1) = "met time coordinate is not contiguous at its native cadence";
   end
   reason = strjoin(reasons, "; ");
end
