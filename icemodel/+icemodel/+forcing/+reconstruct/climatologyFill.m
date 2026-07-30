function [estimate, n_support] = climatologyFill(times, x, query_times, kwargs)
   %CLIMATOLOGYFILL Day-of-year climatology estimate of one channel.
   %
   %  [estimate, n_support] = ...
   %     icemodel.forcing.reconstruct.climatologyFill(times, x, query_times)
   %
   % Role
   %  Station day-of-year median climatology (provenance code 8): the
   %  policy's admission baseline for gaps longer than the short-gap cap,
   %  and an admissible last-resort method where its own gates pass (e.g.
   %  long albedo runs). A smoothing window in day-of-year pools nearby
   %  days so single-year records still yield stable medians; the target's
   %  time-of-day bins are preserved for diurnal channels.
   %  Leap years use a 365-day no-leap calendar (February 29 pools with
   %  February 28) so March-December observations retain the same calendar
   %  bins in every year. Implemented as one precomputed calendar-day x
   %  hour lookup table so a
   %  full-record query axis costs O(n), not O(n^2).
   %
   % Name-value
   %  window_days : half-width of the day-of-year pooling window
   %     (default 7, i.e. a 15-day pool).
   %  min_support : minimum pooled samples for a finite estimate
   %     (default 5); queries with less support return NaN.
   %  diurnal : pool within the query's hour-of-day (default true; set
   %     false for slowly varying channels like albedo).
   %
   % Returns
   %  estimate : double column, one climatology value per query time (NaN
   %     where support is insufficient).
   %  n_support : pooled sample count per query.
   %
   % See also: icemodel.forcing.reconstruct.admissionGate

   arguments
      times datetime
      x (:, 1) double
      query_times datetime
      kwargs.window_days (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().climatology_window_days
      kwargs.min_support (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().climatology_min_support
      kwargs.diurnal (1, 1) logical = true
   end

   % Bin the observed pool once by normalized calendar day and exact posting
   % time. Deriving bins from the data preserves hourly, half-hourly, and
   % quarter-hourly family semantics without a cadence-specific constant.
   finite = isfinite(x);
   pool_doy = normalizedDayOfYear(times(finite));
   pool_clock = timeofday(times(finite));
   query_clock = timeofday(query_times);
   pool_values = x(finite);
   clock_bins = unique([pool_clock(:); query_clock(:)]);
   [~, pool_bin] = ismember(pool_clock, clock_bins);
   [~, query_bin] = ismember(query_clock, clock_bins);
   n_bins = numel(clock_bins);
   if ~kwargs.diurnal
      pool_bin = ones(size(pool_bin));
      query_bin = ones(size(query_bin));
      n_bins = 1;
   end

   % Collect the pooled values per (doy, posting-time) cell, then aggregate each
   % query cell over its circular day window.
   cell_values = accumarray([pool_doy, pool_bin], pool_values, ...
      [365, n_bins], @(v) {v}, {zeros(0, 1)});

   % Precompute the windowed median and support per calendar/time cell once;
   % every query then indexes the table.
   window = -kwargs.window_days:kwargs.window_days;
   table_median = nan(365, n_bins);
   table_support = zeros(365, n_bins);
   for d = 1:365
      rows = mod(d - 1 + window, 365) + 1;
      for b = 1:n_bins
         pooled = vertcat(cell_values{rows, b});
         table_support(d, b) = numel(pooled);
         if numel(pooled) >= kwargs.min_support
            table_median(d, b) = median(pooled);
         end
      end
   end

   query_doy = normalizedDayOfYear(query_times);
   linear = sub2ind([365, n_bins], query_doy(:), query_bin(:));
   estimate = table_median(linear);
   n_support = table_support(linear);
end

function doy = normalizedDayOfYear(times)
   %NORMALIZEDDAYOFYEAR Map dates to a stable 365-day climatology calendar.
   month_offsets = [0, 31, 59, 90, 120, 151, ...
      181, 212, 243, 273, 304, 334].';
   month_number = month(times(:));
   day_number = day(times(:));
   doy = month_offsets(month_number) + day_number;
   doy(month_number == 2 & day_number == 29) = 59;
end
