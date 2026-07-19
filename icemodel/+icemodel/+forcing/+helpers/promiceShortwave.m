function [swd, swu, metadata] = promiceShortwave(aws, kwargs)
   %PROMICESHORTWAVE Select physical public PROMICE shortwave channels.
   %
   %  [swd, swu, metadata] = icemodel.forcing.helpers.promiceShortwave(aws)
   %  [swd, swu, metadata] = ... promiceShortwave(aws, fill_darkness=true, ...
   %     latitude=..., longitude=...)
   %
   % `readPromiceAws` intentionally exposes both the source-faithful raw
   % pyranometer channels (`swd`/`swu`, from `dsr`/`usr`) and pypromice's
   % tilt/bias-corrected products (`swd_cor`/`swu_cor`). The PROMICE variable
   % dictionary allows raw shortwave down to -10 W m-2, so small negative raw
   % values are valid source records but not physical public forcing/evaluation
   % fluxes. Builders therefore prefer the corrected product sample-by-sample,
   % fall back to the raw measurement where the corrected product is missing,
   % and clamp any remaining finite negative selected value to zero.
   %
   % Missing samples remain missing by default. With fill_darkness=true, missing
   % samples are set to physical zero only when the complete hourly interval is
   % below -6 degrees NOAA solar elevation. Twilight and daylight gaps remain
   % missing, and an absent or whole-source-file all-missing channel remains a
   % placeholder. Whole-file support can be supplied when `aws` is a surgical
   % slice so selection is independent of the requested window.
   % The source file remains the authoritative raw-value record, while returned
   % metadata preserves exact source, fallback, negative-input, darkness-fill,
   % and remaining-missing counts for staged artifacts.
   %
   % See also: icemodel.forcing.readPromiceAws,
   %  icemodel.forcing.buildPromiceMet, icemodel.forcing.buildPromiceData

   arguments
      aws timetable
      kwargs.fill_darkness (1, 1) logical = false
      kwargs.latitude (1, 1) double = NaN
      kwargs.longitude (1, 1) double = NaN
      kwargs.swd_source_file_observations_present (1, 1) double = NaN
      kwargs.swu_source_file_observations_present (1, 1) double = NaN
   end

   % Deep civil night is the conservative empirical boundary at KAN_L: finite
   % source SWD is identically zero below it, whereas twilight SWD can be
   % positive. Evaluate the complete source-hour interval in canonical UTC.
   darkness_threshold_degrees = -6;
   darkness_method = ...
      "NOAA solar elevation at UTC interval start, quarter-hours, and end; " + ...
      "zero fill when whole-hour maximum <= -6 degrees";
   deep_dark = false(height(aws), 1);
   if kwargs.fill_darkness
      if ~isfinite(kwargs.latitude) || ~isfinite(kwargs.longitude)
         error('icemodel:forcing:helpers:promiceShortwave:missingLocation', ...
            'fill_darkness=true requires finite latitude and longitude')
      end
      deep_dark = wholeHourDeepCivilNight(aws.Properties.RowTimes, ...
         kwargs.latitude, kwargs.longitude, darkness_threshold_degrees);
   end

   % Normalize downwelling and upwelling consistently so derived net
   % shortwave in observational Data uses one coherent corrected-first pair.
   [swd, swd_metadata] = selectChannel(aws, "swd", "swd_cor", "dsr", ...
      "dsr_cor", deep_dark, kwargs.fill_darkness, ...
      darkness_threshold_degrees, darkness_method, ...
      kwargs.swd_source_file_observations_present);
   [swu, swu_metadata] = selectChannel(aws, "swu", "swu_cor", "usr", ...
      "usr_cor", deep_dark, kwargs.fill_darkness, ...
      darkness_threshold_degrees, darkness_method, ...
      kwargs.swu_source_file_observations_present);

   % Keep metadata flat because writer and audit contracts store scalar
   % artifact metadata rather than nested per-channel structs.
   metadata = struct();
   metadata = copyFields(metadata, swd_metadata);
   metadata = copyFields(metadata, swu_metadata);
end

function [data, metadata] = selectChannel(aws, raw_name, corrected_name, ...
      raw_source, corrected_source, deep_dark, fill_darkness, ...
      darkness_threshold_degrees, darkness_method, ...
      source_file_observations_status)
   %SELECTCHANNEL Prefer corrected finite values, then physical raw fallback.
   names = string(aws.Properties.VariableNames);
   has_raw = ismember(raw_name, names);
   has_corrected = ismember(corrected_name, names);

   % Absent source channels are represented only inside this selection helper;
   % observational builders decide whether an absent public channel is omitted.
   raw = nan(height(aws), 1);
   corrected = nan(height(aws), 1);
   if has_raw
      raw = aws.(raw_name);
   end
   if has_corrected
      corrected = aws.(corrected_name);
   end

   % The pypromice correction is the preferred public product. Raw fallback
   % retains source coverage without allowing a negative radiative flux.
   data = corrected;
   raw_fallback = ~isfinite(data) & isfinite(raw);
   data(raw_fallback) = raw(raw_fallback);
   data(~isfinite(data)) = NaN;
   negative = isfinite(data) & data < 0;
   data(negative) = 0;

   % A source-backed channel may safely replace missing deep-night radiation
   % with physical zero. Whole-file support makes this decision independent of
   % window selection, while absent/all-missing source files remain placeholders.
   window_observations_present = ...
      any(isfinite(raw) | isfinite(corrected), 'all');
   source_file_observations_present = window_observations_present;
   if ~isnan(source_file_observations_status)
      source_file_observations_present = window_observations_present ...
         || logical(source_file_observations_status);
   end
   darkness_fill = false(height(aws), 1);
   if fill_darkness && (has_raw || has_corrected) ...
         && source_file_observations_present
      darkness_fill = ~isfinite(data) & deep_dark;
      data(darkness_fill) = 0;
   end

   % Preserve quantitative provenance for both the raw product and the exact
   % public selection. This is sufficient to recover rejected samples from the
   % named source file without staging a physically invalid public channel.
   raw_finite = isfinite(raw);
   raw_minimum = NaN;
   if any(raw_finite, 'all')
      raw_minimum = min(raw(raw_finite), [], 'all');
   end
   prefix = char(raw_name);
   metadata = struct();
   metadata.([prefix '_source_present']) = has_raw;
   metadata.([prefix '_corrected_source_present']) = has_corrected;
   metadata.([prefix '_observations_present']) = window_observations_present;
   metadata.([prefix '_source_file_observations_present']) = ...
      source_file_observations_present;
   metadata.([prefix '_corrected_used_count']) = nnz(isfinite(corrected));
   metadata.([prefix '_raw_fallback_count']) = nnz(raw_fallback);
   metadata.([prefix '_raw_negative_count']) = nnz(raw_finite & raw < 0);
   metadata.([prefix '_corrected_negative_count']) = ...
      nnz(isfinite(corrected) & corrected < 0);
   metadata.([prefix '_raw_fallback_negative_count']) = ...
      nnz(raw_fallback & raw < 0);
   metadata.([prefix '_negative_clamped_count']) = nnz(negative);
   metadata.([prefix '_darkness_fill_enabled']) = fill_darkness;
   metadata.([prefix '_darkness_fill_method']) = darkness_method;
   metadata.([prefix '_darkness_threshold_degrees']) = ...
      darkness_threshold_degrees;
   metadata.([prefix '_darkness_zero_filled_count']) = nnz(darkness_fill);
   metadata.([prefix '_remaining_missing_count']) = nnz(~isfinite(data));
   metadata.([prefix '_raw_minimum']) = raw_minimum;
   metadata.([prefix '_raw_source']) = ...
      "PROMICE L3 " + raw_source + " [W m-2]";
   metadata.([prefix '_corrected_source']) = ...
      "PROMICE L3 " + corrected_source + " [W m-2]";

   % A channel-specific policy string also makes a whole-file all-missing or
   % absent required met channel an intentional placeholder to artifact QA.
   if ~(has_raw || has_corrected) || ~source_file_observations_present
      metadata.([prefix '_policy']) = raw_name + ...
         " = NaN placeholder (PROMICE L3 " + raw_source + " and " + ...
         corrected_source + ...
         " are absent or all missing in the source file); no observations invented";
   elseif ~window_observations_present
      metadata.([prefix '_policy']) = raw_name + ...
         " has no finite PROMICE L3 source observations in the requested window";
   elseif any(isfinite(corrected), 'all')
      metadata.([prefix '_policy']) = raw_name + ...
         " uses PROMICE L3 " + corrected_source + ...
         " where finite, falls back to " + raw_source + ...
         ", and clamps any remaining finite negative selected value to zero";
   else
      metadata.([prefix '_policy']) = raw_name + ...
         " uses PROMICE L3 " + raw_source + ...
         " and clamps finite negative sensor-offset values to zero";
   end
   if fill_darkness && source_file_observations_present
      metadata.([prefix '_policy']) = metadata.([prefix '_policy']) + ...
         "; missing whole-hour deep-civil-night samples are set to zero; " + ...
         "twilight and daylight missing samples remain missing";
   end
end

function deep_dark = wholeHourDeepCivilNight( ...
      Time, latitude, longitude, threshold_degrees)
   %WHOLEHOURDEEPCIVILNIGHT Identify bins wholly below civil twilight.

   % PROMICE times label hourly interval starts. Canonicalize zoned inputs by
   % instant, while interpreting the reader's unzoned UTC timestamps as UTC.
   Time.TimeZone = 'UTC';
   solar_times = Time + minutes(0:15:60);
   solar_elevation = icemodel.forcing.helpers.solarElevation( ...
      solar_times, latitude, longitude);
   deep_dark = max(solar_elevation, [], 2) <= threshold_degrees;
end

function target = copyFields(target, source)
   %COPYFIELDS Copy one scalar metadata struct into another.
   fields = fieldnames(source);
   for k = 1:numel(fields)
      target.(fields{k}) = source.(fields{k});
   end
end
