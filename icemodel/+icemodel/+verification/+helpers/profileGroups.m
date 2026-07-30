function groups = profileGroups(value, kwargs)
   %PROFILEGROUPS Split profile rows by stable source identity and UTC date.
   %
   %  groups = icemodel.verification.helpers.profileGroups(value)
   %  groups = icemodel.verification.helpers.profileGroups( ...
   %     value, time_resolution="timestamp")
   %
   % MAR tables use their explicit profile_id. SUMup tables use the source
   % name_key that its builder already treats as physical profile identity,
   % falling back to name or coordinates for older tables. Row depth and
   % measurement/reference/method identifiers are deliberately excluded so one
   % physical profile remains one group. Tables without a datetime/profile
   % identity remain one ungrouped profile for backward-compatible plots.

   arguments
      value
      kwargs.time_resolution (1, 1) string ...
         {mustBeMember(kwargs.time_resolution, ["date", "timestamp"])} ...
         = "date"
   end

   if ~(istable(value) || istimetable(value))
      error('icemodel:verification:profileGroups:badInput', ...
         'profile data must be a table or timetable')
   end

   % Empty inputs have no physical profile groups.
   prototype = struct('id', "", ...
      'datetime', NaT(1, 1, 'TimeZone', 'UTC'), 'data', table());
   groups = repmat(prototype, 0, 1);
   if isempty(value)
      return
   end

   names = string(value.Properties.VariableNames);
   datetime_name = ["datetime", "date", "time", "Time"];
   datetime_name = datetime_name(find(ismember(datetime_name, names), 1));
   has_datetime = ~isempty(datetime_name) ...
      && isdatetime(value.(char(datetime_name)));
   has_row_times = istimetable(value) ...
      && isdatetime(value.Properties.RowTimes);

   % Explicit MAR identity takes precedence. SUMup itself counts physical
   % profiles by name_key, so use that stable site/core key plus date rather
   % than row-level reference or method provenance. Those keys can legitimately
   % vary down one composite core and would split its depth support into
   % one-row groups. Older tables fall back to the resolved name, then location.
   if ismember("profile_id", names)
      identity = string(value.profile_id);
   elseif ismember("name_key", names) ...
         && any(~ismissing(value.name_key), 'all')
      identity = "name_key|" + string(value.name_key);
   elseif ismember("name", names) ...
         && any(~ismissing(string(value.name)) ...
         & strlength(string(value.name)) > 0, 'all')
      identity = "name|" + string(value.name);
   else
      identity = repmat("profile", height(value), 1);
      if all(ismember(["latitude", "longitude"], names))
         identity = identity + "|" + string(value.latitude) ...
            + "|" + string(value.longitude);
      end
   end

   % A dated profile is keyed by its UTC calendar date. Time-of-day remains in
   % the rows and output metadata, but cannot split one daily MAR snapshot.
   if has_datetime
      timestamps = icemodel.verification.setup.ensureUtc( ...
         value.(char(datetime_name)));
   elseif has_row_times
      timestamps = icemodel.verification.setup.ensureUtc( ...
         value.Properties.RowTimes);
   else
      timestamps = NaT(height(value), 1, 'TimeZone', 'UTC');
   end
   dates = dateshift(timestamps, 'start', 'day');
   if has_datetime || has_row_times
      if kwargs.time_resolution == "timestamp"
         time_key = string(timestamps, ...
            'yyyy-MM-dd''T''HH:mm:ss.SSS');
      else
         time_key = string(dates, 'yyyy-MM-dd');
      end
      keys = identity + "|" + time_key;
   else
      keys = identity;
   end

   % Stable grouping preserves source order for deterministic plotting and
   % additive merge behavior.
   [unique_keys, first_rows] = unique(keys, 'stable');
   [~, group_index] = ismember(keys, unique_keys);
   groups = repmat(prototype, numel(unique_keys), 1);
   for n = 1:numel(unique_keys)
      keep = group_index == n;
      groups(n).id = identity(first_rows(n));
      if has_datetime || has_row_times
         groups(n).datetime = timestamps(first_rows(n));
      end
      groups(n).data = value(keep, :);
   end
end
