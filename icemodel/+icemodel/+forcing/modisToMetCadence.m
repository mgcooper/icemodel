function [albedo, support] = modisToMetCadence(albedo_daily, time_daily, ...
      time_met, kwargs)
   %MODISTOMETCADENCE Attach staged daily MODIS albedo onto a met time axis.
   %
   %  [albedo, support] = icemodel.forcing.modisToMetCadence( ...
   %     albedo_daily, time_daily, time_met)
   %  [albedo, support] = ... modisToMetCadence(_, max_gap=days(5))
   %
   % Single source of the daily-MODIS -> met-cadence conversion mandated by
   % the gap-fill reconstruction policy (reconstruct/POLICY.md B12, D-15):
   % attachment happens at staging/reconstruction time through THIS helper,
   % never at model runtime (icemodel.loadmet rejects daily swap userdata),
   % so every consumer resolves the staged daily artifact identically. Met
   % samples are linear interpolations between finite daily samples; daily
   % satellite albedo is adequate sub-daily forcing because albedo varies
   % with snowfall and synoptic events, not the hour of day (B12).
   %
   % Gap awareness: the staged daily series carries NaN holes (masked 999
   % sentinels, polar darkness, bounds-masked samples). Interpolation never
   % bridges a pair of finite daily samples separated by more than MAX_GAP
   % (default days(5)); met samples inside a wider hole, and met samples
   % outside the finite daily support (no extrapolation), return NaN with
   % SUPPORT false. Met samples that land exactly on a finite daily sample
   % are always supported.
   %
   % The reconstruction albedo bounds (physicalBounds("albedo")) are
   % enforced fail-closed on the finite daily input: staged artifacts are
   % already masked to those bounds, so an out-of-bounds sample here means
   % the caller bypassed staging and the conversion refuses to run.
   %
   % Inputs
   %  albedo_daily - daily MODIS albedo [-], NaN where missing
   %  time_daily   - strictly increasing datetimes of the daily samples
   %  time_met     - target met time axis (any cadence); when its TimeZone
   %                 differs from TIME_DAILY the met axis wins, matching the
   %                 loadmet swap-data convention
   %
   % Name-value
   %  max_gap - widest finite-sample spacing the linear bridge may cross
   %            (duration, default days(5))
   %
   % Outputs
   %  albedo  - albedo on TIME_MET [-], NaN where unsupported
   %  support - logical column, true where ALBEDO carries a MODIS-derived
   %            value
   %
   % See also: icemodel.forcing.stageModisAlbedo,
   %  icemodel.forcing.reconstruct.physicalBounds,
   %  icemodel.forcing.helpers.dailyToHourly, icemodel.loadmet

   arguments
      albedo_daily (:, 1) double
      time_daily (:, 1) datetime
      time_met (:, 1) datetime
      kwargs.max_gap (1, 1) duration = days(5)
   end

   % An unordered or duplicated daily axis makes interval bookkeeping (and
   % interp1) meaningless, so fail loudly instead of silently sorting.
   if numel(time_daily) ~= numel(albedo_daily)
      error('icemodel:forcing:modisToMetCadence:sizeMismatch', ...
         'albedo_daily (%d) and time_daily (%d) must align', ...
         numel(albedo_daily), numel(time_daily))
   end
   if any(diff(time_daily) <= 0)
      error('icemodel:forcing:modisToMetCadence:unsortedDailyAxis', ...
         'time_daily must be strictly increasing')
   end
   % A non-positive gap allowance would refuse even consecutive samples.
   if kwargs.max_gap <= days(0)
      error('icemodel:forcing:modisToMetCadence:invalidMaxGap', ...
         'max_gap must be a positive duration')
   end

   % Follow the loadmet swap convention (swap data adopts mettime.TimeZone):
   % the met axis owns the time zone so mixed zoned/unzoned callers convert
   % here, at the single attachment boundary, not ad hoc downstream.
   if ~strcmp(time_daily.TimeZone, time_met.TimeZone)
      time_daily.TimeZone = time_met.TimeZone;
   end

   albedo = nan(numel(time_met), 1);

   % Interpolate between finite daily samples only: NaN holes are handled by
   % the explicit gap policy below rather than by interp1 NaN propagation,
   % which would refuse every hole regardless of width.
   finite = isfinite(albedo_daily);
   time_finite = time_daily(finite);
   value_finite = albedo_daily(finite);
   if isempty(time_finite)
      % No finite support anywhere: the whole met axis is unsupported.
      support = false(numel(time_met), 1);
      return
   end
   if isscalar(time_finite)
      % One finite sample defines no interval, so only met samples landing
      % exactly on that day carry a value.
      on_sample = time_met == time_finite;
      albedo(on_sample) = value_finite;
      support = on_sample;
      return
   end

   % Reuse the canonical daily->met interpolation without endpoint
   % extrapolation so met samples outside the finite support stay NaN, and
   % enforce the reconstruction albedo bounds (SSOT) fail-closed on the
   % staged input rather than restating the limits here.
   albedo = icemodel.forcing.helpers.dailyToHourly( ...
      value_finite, time_finite, time_met, extrapolate=false, ...
      bounds=icemodel.forcing.reconstruct.physicalBounds("albedo"));

   % Refuse bridges across holes wider than max_gap: locate each met sample's
   % bracketing finite-sample interval and blank the ones whose interval is
   % too wide. Samples exactly on a finite daily sample stay supported even
   % when their neighbouring intervals are wide.
   spacing = diff(time_finite);
   interval = discretize(time_met, time_finite);
   on_sample = ismember(time_met, time_finite);
   inside = ~isnan(interval) & ~on_sample;
   wide = false(numel(time_met), 1);
   wide(inside) = spacing(interval(inside)) > kwargs.max_gap;
   albedo(wide) = NaN;

   % The support mask is exactly the finite output: in-support linear values
   % of finite operands are finite, and every refusal above wrote NaN.
   support = isfinite(albedo);
end
