function [window_start, window_end, enabled] = pairedWindow(startdate, enddate)
   %PAIREDWINDOW Normalize one optional start/end pair to UTC datetimes.
   %
   %  [window_start, window_end, enabled] = ...
   %     icemodel.internal.pairedWindow(startdate, enddate)
   %
   % Blank bounds disable the window and return UTC-zoned NaT values. A
   % supplied window must contain both finite scalar bounds in chronological
   % order. Keeping this policy in the neutral internal namespace lets forcing
   % builders and verification importers share one public-input boundary.

   % Determine endpoint presence before parsing so half-window errors are
   % independent of whether callers use strings, chars, or datetimes.
   has_start = hasBound(startdate);
   has_end = hasBound(enddate);
   if has_start ~= has_end
      invalidWindow('startdate and enddate must be provided together')
   end

   % Return UTC NaT sentinels for the all-available/default path.
   enabled = has_start;
   window_start = NaT('TimeZone', 'UTC');
   window_end = NaT('TimeZone', 'UTC');
   if ~enabled
      return
   end

   % Parse and normalize once at the public boundary. Zoned datetime inputs
   % preserve their instant while being represented in UTC.
   try
      window_start = normalizeUtc(startdate);
      window_end = normalizeUtc(enddate);
   catch
      invalidWindow( ...
         'startdate and enddate must be parseable datetime bounds')
   end
   if ~isscalar(window_start) || ~isscalar(window_end) ...
         || isnat(window_start) || isnat(window_end) ...
         || window_end < window_start
      invalidWindow( ...
         'startdate and enddate must be finite scalar bounds in chronological order')
   end
end

function tf = hasBound(value)
   %HASBOUND True for a finite datetime or a nonblank datetime string.
   if isempty(value)
      tf = false;
   elseif isdatetime(value)
      tf = any(~isnat(value(:)));
   else
      text = strtrim(string(value));
      tf = any(strlength(text(:)) > 0 & text(:) ~= "NaT");
   end
end

function value = normalizeUtc(value)
   %NORMALIZEUTC Parse text or convert a datetime instant to UTC.
   if isdatetime(value)
      % Assigning TimeZone interprets an unzoned wall time as UTC and converts a
      % zoned datetime's existing instant into its UTC representation.
      value.TimeZone = 'UTC';
   else
      value = datetime(value, 'TimeZone', 'UTC');
   end
end

function invalidWindow(message)
   %INVALIDWINDOW Raise the single neutral optional-window contract error.
   error('icemodel:internal:pairedWindow:invalidWindow', '%s', message)
end
