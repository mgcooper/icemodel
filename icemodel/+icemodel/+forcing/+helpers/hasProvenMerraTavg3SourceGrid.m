function tf = hasProvenMerraTavg3SourceGrid(T, metadata)
   %HASPROVENMERRATAVG3SOURCEGRID True for an exact native tavg3 inventory.
   %
   %  tf = icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(T, metadata)
   %
   % Validates the durable native-grid proof used to recover legacy MERRA glc
   % channels. The proof records every absent UTC 00/03/... source stamp, so a
   % saved value at an omitted source time cannot be mistaken for a native row.
   % Artifacts without runoff/albedo/snowd/swe do not consume tavg3 and require
   % no glacier-collection grid proof.

   tf = istimetable(T) && isstruct(metadata) && isscalar(metadata);
   if ~tf
      return
   end
   channels = intersect(["runoff", "albedo", "snowd", "swe"], ...
      string(T.Properties.VariableNames), 'stable');
   if isempty(channels)
      return
   end

   % Source stamps must be identifiable on a sorted UTC application axis.
   row_times = T.Properties.RowTimes;
   if ~isdatetime(row_times) || isempty(row_times) ...
         || string(row_times.TimeZone) ~= "UTC" || ~issorted(row_times)
      tf = false;
      return
   end
   expected = row_times(minute(row_times) == 0 & second(row_times) == 0 ...
      & mod(hour(row_times), 3) == 0);
   if isempty(expected)
      tf = false;
      return
   end

   % Counts and the explicit missing-time list must describe this exact grid.
   required = ["merra_tavg3_source_grid_policy", ...
      "merra_tavg3_expected_source_row_count", ...
      "merra_tavg3_source_row_count", ...
      "merra_tavg3_source_time_gap_count", ...
      "merra_tavg3_missing_source_times"];
   if ~all(isfield(metadata, required)) ...
         || ~isequal(string(metadata.merra_tavg3_source_grid_policy), ...
         "native_glc_timestamp_inventory")
      tf = false;
      return
   end
   expected_count = double(metadata.merra_tavg3_expected_source_row_count);
   source_count = double(metadata.merra_tavg3_source_row_count);
   gap_count = double(metadata.merra_tavg3_source_time_gap_count);
   if ~isscalar(expected_count) || ~isscalar(source_count) ...
         || ~isscalar(gap_count)
      tf = false;
      return
   end
   counts = [expected_count, source_count, gap_count];
   if any(~isfinite(counts)) || any(counts < 0) ...
         || any(counts ~= fix(counts)) ...
         || expected_count ~= numel(expected) ...
         || source_count + gap_count ~= expected_count
      tf = false;
      return
   end

   % Missing stamps are UTC, unique, sorted, and drawn from the expected grid.
   missing = metadata.merra_tavg3_missing_source_times;
   if ~isdatetime(missing) || string(missing.TimeZone) ~= "UTC"
      tf = false;
      return
   end
   missing = missing(:);
   if numel(missing) ~= gap_count || ~issorted(missing) ...
         || numel(unique(missing)) ~= numel(missing) ...
         || ~all(ismember(missing, expected))
      tf = false;
   end
end
