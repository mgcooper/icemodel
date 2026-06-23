function [corrected, record, flags] = destepSurface(t, surf, kwargs)
   %DESTEPSURFACE Detect (and optionally correct) step-shifts in a surface series.
   %
   %  [corrected, record, flags] = icemodel.forcing.destepSurface(t, surf)
   %  [...] = icemodel.forcing.destepSurface(t, surf, mode="unambiguous", ...
   %     transition_times=..., max_rate=..., tol_days=..., season="ablation")
   %
   % An OPT-IN, evaluation-time transform for PROMICE/GC-Net cumulative
   % surface-height series (ablation [+down] or surface_height [+up]). It
   % DETECTS candidate single-timestep step-shifts using multiple independent
   % lines of evidence, CLASSIFIES each as UNAMBIGUOUS or AMBIGUOUS, and (per
   % the requested mode) levels UNAMBIGUOUS steps by subtracting the offset from
   % the post-step segment. The staged .mat is never altered by this function;
   % buildPromiceData stages the raw series plus per-sample flags, and a
   % consumer calls this transform at analysis time.
   %
   % DETECTION (a candidate step is a single finite-to-finite jump d = surf(k+1)
   % - surf(k) between adjacent FINITE samples; the gap-bridged interior is not a
   % step). A candidate must first clear the per-step RATE bound (magnitude
   % line), then carries up to four evidence lines:
   %   1. magnitude   |d| exceeds the physically-plausible per-step bound implied
   %                  by max_rate over the inter-sample interval (a jump no real
   %                  ablation/accumulation rate can produce in that time). This
   %                  is the GATE: only rate-implausible jumps are candidates.
   %   2. implausible |d| exceeds an absolute single-step ceiling (max_step): a
   %                  jump so large no physical process produces it across a
   %                  single sampling interval (a sensor reset / re-survey /
   %                  installation jump), independent of the rate test.
   %   3. transition  the jump time falls within tol_days of a known station
   %                  handover (transition_times), the documented cause of an
   %                  expected discrete offset when a site merges several AWS.
   %   4. season      a melt-signed jump in the accumulation/winter season is
   %                  suspect: ice ablation occurs in the melt season, so a
   %                  melt-signed step in winter violates season consistency. To
   %                  avoid promoting hourly noise, the season line requires the
   %                  jump to also clear max_step (a winter MICRO-jump is noise).
   %
   % CLASSIFICATION
   %   UNAMBIGUOUS  the rate-magnitude gate fires AND at least one independent
   %                corroborating line (implausible, transition, or season)
   %                agrees -> a single physically-impossible jump with an
   %                independent cause. Corrected when mode allows.
   %   AMBIGUOUS    the rate-magnitude gate alone fires (a jump merely faster
   %                than the smooth bound but neither grossly implausible, nor
   %                near a transition, nor season-inconsistent) -> FLAGGED, never
   %                auto-corrected. Hourly noise lands here.
   %
   % CORRECTION levels the post-step segment to the pre-step level by subtracting
   % the cumulative offset of the applied steps from all later samples, so the
   % de-stepped series preserves the within-segment shape and only removes the
   % artificial discontinuity. Ambiguous steps are left in place (optionally
   % censored to NaN via mode="censor_ambiguous").
   %
   % Inputs
   %  t    - datetime column, the series time axis
   %  surf - double column, the cumulative surface-height series (may hold NaN)
   %
   % Name-value
   %  mode             : "detect"             detect only, no correction
   %                     "unambiguous"        (DEFAULT) correct unambiguous steps
   %                     "all"                correct unambiguous AND ambiguous
   %                     "censor_ambiguous"   correct unambiguous, NaN ambiguous
   %  transition_times : datetime array of known station-transition times ([])
   %  max_rate         : max plausible |surface rate| [m/day] (default 0.6;
   %                     covers extreme bare-ice melt ~0.1-0.2 m/day with margin)
   %  max_step         : absolute single-step ceiling [m] for the implausible
   %                     line; |d| above it is physically impossible across one
   %                     interval regardless of duration (default 1.0)
   %  tol_days         : transition-coincidence tolerance [days] (default 14)
   %  season           : "ablation" (+down series) or "accumulation" (+up); sets
   %                     the melt-signed direction for the season evidence line
   %
   % Outputs
   %  corrected - double column, the (possibly) de-stepped series; equals surf
   %              when mode="detect" or no step is corrected
   %  record    - struct array, one entry per detected candidate step:
   %              .time .magnitude .evidence (string array of fired lines)
   %              .classification ("unambiguous"/"ambiguous") .applied (logical)
   %  flags     - struct of per-sample logical/double masks aligned to surf:
   %              .step_detected   1 at the post-step sample of any candidate
   %              .step_correctable 1 at the post-step sample of an unambiguous
   %                                (correctable) candidate
   %              .step_magnitude  signed jump magnitude at each detected sample,
   %                                0 elsewhere
   %              .step_reason     string per sample: "" or the evidence summary
   %
   % See also: icemodel.forcing.buildPromiceData,
   %  icemodel.forcing.helpers.surfaceFlags

   arguments
      t (:, 1) datetime
      surf (:, 1) double
      kwargs.mode (1, 1) string {mustBeMember(kwargs.mode, ...
         ["detect", "unambiguous", "all", "censor_ambiguous"])} = "unambiguous"
      kwargs.gap_flag (:, 1) double = []
      kwargs.transition_times (:, 1) datetime = datetime.empty(0, 1)
      kwargs.max_rate (1, 1) double {mustBePositive} = 0.6
      kwargs.max_step (1, 1) double {mustBePositive} = 1.0
      kwargs.tol_days (1, 1) double {mustBeNonnegative} = 14
      kwargs.season (1, 1) string ...
         {mustBeMember(kwargs.season, ["ablation", "accumulation"])} = "ablation"
   end

   n = numel(surf);
   % Per-sample output masks, all zero/empty until a candidate is recorded.
   step_detected = false(n, 1);
   step_correctable = false(n, 1);
   step_magnitude = zeros(n, 1);
   step_reason = strings(n, 1);
   record = emptyEntry();
   record = record([]);   % 1x0 record with the canonical field order

   % A candidate is a jump between ADJACENT DIRECT-OBSERVATION samples. A sample
   % is usable only when it is finite AND not gap-bridged (gap_flag == 0): the
   % slope-bridged / interpolated surface over a data gap is not a direct
   % observation (readme), so a "jump" landing on a bridged sample is an
   % interpolation artifact, not a discrete step. When gap_flag is omitted all
   % finite samples are usable (back-compatible).
   usable = isfinite(surf);
   if ~isempty(kwargs.gap_flag)
      usable = usable & (kwargs.gap_flag(:) == 0);
   end
   finite_idx = find(usable);
   if numel(finite_idx) < 2
      corrected = surf;
      flags = packFlags(step_detected, step_correctable, ...
         step_magnitude, step_reason);
      return
   end

   % For ablation [+down] the melt-signed step is POSITIVE (surface lowering ->
   % ablation increases); for accumulation [+up] melt is NEGATIVE (surface drop).
   melt_sign = 1;
   if kwargs.season == "accumulation"
      melt_sign = -1;
   end

   % Pass 1 (vectorized over adjacent finite pairs): score every jump against
   % the three evidence lines, so the candidate count is known before allocating
   % the record (no growing-array pattern). k0/k1 index the pre/post samples.
   k0 = finite_idx(1:end - 1);
   k1 = finite_idx(2:end);
   d = surf(k1) - surf(k0);

   % Plausible per-step bound from max_rate over the ACTUAL interval (a longer
   % gap legitimately accumulates more change, so the bound scales with time).
   dt_days = max(days(t(k1) - t(k0)), eps);
   bound = kwargs.max_rate * dt_days;

   % Evidence 1 (GATE): magnitude beyond the physically-plausible per-step rate
   % bound. Only rate-implausible jumps are candidates (a within-rate jump is
   % ordinary surface change, not a step-shift).
   mag_fired = abs(d) > bound;

   % Evidence 2: gross single-step implausibility. A jump above max_step cannot
   % be produced by any physical process across one interval (sensor reset /
   % re-survey / installation jump), independent of the rate test.
   implausible_fired = abs(d) > kwargs.max_step;

   % Evidence 3: coincidence (within tol_days) with a known station transition.
   trans_fired = false(size(d));
   if ~isempty(kwargs.transition_times)
      for c = 1:numel(d)
         near = min(abs(days(t(k1(c)) - kwargs.transition_times)));
         trans_fired(c) = near <= kwargs.tol_days;
      end
   end

   % Evidence 4: a melt-signed jump in the non-melt season (Nov..Apr) is
   % season-inconsistent (ice ablation requires melt energy, absent in winter),
   % so a winter surface-lowering step is physically implausible. Now that
   % gap-bridged samples are excluded, a winter melt-signed jump beyond the rate
   % bound is strong evidence; it only needs a small absolute floor (winter_floor)
   % to keep sub-decimetre sensor noise from promoting itself - it no longer
   % needs to clear the gross max_step ceiling.
   winter_floor = 0.3;
   mon = month(t(k1));
   in_winter = mon >= 11 | mon <= 4;
   season_fired = mag_fired & in_winter & (sign(d) == melt_sign) ...
      & (abs(d) > winter_floor);

   % UNAMBIGUOUS = the rate gate fires AND an independent corroborating line
   % (implausible, transition, or season); the rate gate alone -> AMBIGUOUS
   % (flagged, never auto-corrected). Hourly noise lands here.
   cand = find(mag_fired);
   n_cand = numel(cand);
   corroborated = implausible_fired(cand) | trans_fired(cand) | season_fired(cand);

   % Pass 2: allocate the record at the known candidate count and fill it, and
   % set the per-sample masks at each candidate's post-step sample.
   record = repmat(emptyEntry(), 1, n_cand);
   for c = 1:n_cand
      i = cand(c);
      post = k1(i);
      evidence = ["magnitude"; ...
         repmat("implausible", implausible_fired(i), 1); ...
         repmat("transition", trans_fired(i), 1); ...
         repmat("season", season_fired(i), 1)];
      if corroborated(c)
         klass = "unambiguous";
      else
         klass = "ambiguous";
      end
      step_detected(post) = true;
      step_magnitude(post) = d(i);
      step_reason(post) = strjoin(evidence, "+");
      step_correctable(post) = corroborated(c);
      record(c) = struct('time', t(post), 'magnitude', d(i), ...
         'evidence', {evidence}, 'classification', klass, 'applied', false);
   end

   % Apply corrections per the mode. Leveling subtracts the applied step offset
   % from every later sample, removing the discontinuity while preserving the
   % within-segment shape. Detect-only returns surf unchanged.
   corrected = surf;
   switch kwargs.mode
      case "detect"
         apply_mask = false(1, n_cand);
      case "unambiguous"
         apply_mask = corroborated(:)';
      case {"all", "censor_ambiguous"}
         apply_mask = true(1, n_cand);
   end

   % Walk candidates in time order, leveling each applied step's later segment to
   % the pre-step baseline.
   for c = 1:n_cand
      if ~apply_mask(c)
         continue
      end
      record(c).applied = true;
      later = t >= record(c).time;
      corrected(later) = corrected(later) - record(c).magnitude;
   end

   % Censor mode: blank out ambiguous steps' post-step samples up to the next
   % candidate (the ambiguous discontinuity is removed from scoring, not faked).
   if kwargs.mode == "censor_ambiguous"
      for r = 1:numel(record)
         if record(r).classification == "ambiguous"
            corrected(t == record(r).time) = NaN;
         end
      end
   end

   flags = packFlags(step_detected, step_correctable, ...
      step_magnitude, step_reason);
end

%% Local functions
function entry = emptyEntry()
   %EMPTYENTRY One step-record entry template (canonical field order).
   %
   % Scalar so repmat can preallocate the record at the known candidate count;
   % record([]) yields the 1x0 form for the no-candidate case.
   entry = struct('time', NaT, 'magnitude', 0, 'evidence', strings(0, 1), ...
      'classification', "", 'applied', false);
end

function flags = packFlags(detected, correctable, magnitude, reason)
   %PACKFLAGS Assemble the per-sample flag struct in canonical field order.
   flags = struct( ...
      'step_detected', double(detected), ...
      'step_correctable', double(correctable), ...
      'step_magnitude', magnitude, ...
      'step_reason', reason);
end
