function [x, provenance, audit, quality] = smoothShortwaveSeams( ...
      times, x, provenance, latitude, longitude, kwargs)
   %SMOOTHSHORTWAVESEAMS Repair empirical outlier boundaries in filled SWD.
   %
   %  [x, provenance, audit, quality] = ...
   %     icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
   %     times, x, provenance, latitude, longitude)
   %
   % Role
   %  Post-final SWD quality repair under POLICY D-32. This function
   %  compares method boundaries with the observed hourly-step distribution
   %  of the station in the same season and solar-elevation band. A
   %  boundary above the configured percentile and above the adjacent
   %  same-direction local slope gets a repair. The repair masks only the
   %  reconstructed posting beside the boundary and reconnects the
   %  remaining anchors with a straight line. Native observations stay
   %  unchanged. The local-slope floor and the observed-anchor bridge floor
   %  exempt smooth solar ramps and changes that no interpolation can make
   %  smaller. The repair never exceeds the short-gap cap, and never
   %  crosses deep darkness or a season.
   %
   %  One-posting repairs run only up to max_passes. KANL evidence showed
   %  that one posting reduced 75 outliers to 3. A two-posting window
   %  changed more data with no further benefit. The second one-posting
   %  pass reduced the count to the expected tail rate.
   %
   % Name-value
   %  percentile, min_reference_steps, max_passes, cap_hours : defaults
   %     from the central reconstruction options.
   %
   % Returns
   %  x, provenance : repaired source-posting series and codes.
   %  audit : bounded_interp audit rows for changed segments.
   %  quality : one-row table with before/after boundary diagnostics.

   arguments
      times (:, 1) datetime
      x (:, 1) double
      provenance (:, 1) uint8
      latitude (1, 1) double {mustBeFinite}
      longitude (1, 1) double {mustBeFinite}
      kwargs.percentile (1, 1) double ...
         {mustBeInRange(kwargs.percentile, 50, 100)} = ...
         icemodel.forcing.reconstruct.setopts().seam_qa_percentile
      kwargs.min_reference_steps (1, 1) double ...
         {mustBeInteger, mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().seam_qa_min_reference_steps
      kwargs.max_passes (1, 1) double ...
         {mustBeInteger, mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().seam_qa_max_passes
      kwargs.cap_hours (1, 1) double ...
         {mustBePositive, ...
         icemodel.forcing.reconstruct.mustBeCapHours(kwargs.cap_hours)} = ...
         icemodel.forcing.reconstruct.setopts().cap_hours
   end

   if numel(times) ~= numel(x) || numel(x) ~= numel(provenance)
      error('icemodel:reconstruct:smoothShortwaveSeams:sizeMismatch', ...
         'times, x, and provenance must share one sample axis');
   end
   if numel(times) < 2
      error('icemodel:reconstruct:smoothShortwaveSeams:shortAxis', ...
         'at least two source postings are required for seam QA');
   end

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   original = x;
   original_provenance = provenance;
   native = nativeShortwave(original_provenance, codes);
   protected = native ...
      | original_provenance == codes.twilight_climatology ...
      | original_provenance == codes.darkness;
   % One audit slot per QA pass, filled with that pass's repair rows and
   % concatenated once at the end. Each slot starts as the typed empty column,
   % so a pass that repairs nothing adds nothing.
   pass_audits = repmat({cell(0, 1)}, kwargs.max_passes, 1);
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   dt_hours = hours(median(diff(times)));
   elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, latitude, longitude, hours(dt_hours));
   regime = 1 + (elevation > bands.civil_twilight_deg) ...
      + (elevation >= bands.horizon_deg);
   [initial_outliers, initial_ratio, n_boundaries, ...
      n_reference_steps, initial_unsupported] = boundaryOutliers( ...
      times, x, provenance, elevation, regime, ...
      kwargs.percentile, kwargs.min_reference_steps, codes);
   initial_max_ratio = max([0; initial_ratio(initial_outliers)]);
   changed = false(size(x));
   n_passes = 0;

   % Each pass repairs only the reconstructed posting immediately beside
   % a flagged transition. A recomputed screen finds a rare cascaded
   % boundary, and no repair becomes wider than one posting.
   for pass = 1:kwargs.max_passes
      [outliers, ~] = boundaryOutliers(times, x, provenance, ...
         elevation, regime, kwargs.percentile, ...
         kwargs.min_reference_steps, codes);
      if ~any(outliers)
         break
      end
      boundary = find(outliers);
      mask = false(size(x));
      current_protected = nativeShortwave(provenance, codes) ...
         | provenance == codes.twilight_climatology ...
         | provenance == codes.darkness;
      right_idx = boundary + 1;
      left_editable = ~current_protected(boundary);
      right_editable = ~current_protected(right_idx);
      selected = boundary;
      selected(~left_editable & right_editable) = ...
         right_idx(~left_editable & right_editable);
      both_editable = left_editable & right_editable;
      selected(both_editable) = right_idx(both_editable);
      selected(~left_editable & ~right_editable) = [];
      mask(selected) = true;
      if ~any(mask)
         break
      end

      % Reconnect each resulting mask run on its own. A run that fails
      % stays unchanged, so QA reports it and the repair does not become
      % wider.
      proposal = x;
      proposal(mask) = NaN;
      repaired_this_pass = false(size(x));
      edges = diff([false; mask; false]);
      starts = find(edges == 1);
      stops = find(edges == -1) - 1;
      % One slot per mask run, for the same reason as the per-pass slots.
      run_audits = repmat({cell(0, 1)}, numel(starts), 1);
      for g = 1:numel(starts)
         target = (starts(g):stops(g)).';
         before = target(1) - 1;
         after = target(end) + 1;
         if numel(target) * dt_hours > kwargs.cap_hours ...
               || before < 1 || after > numel(x) ...
               || ~isfinite(proposal(before)) ...
               || ~isfinite(proposal(after))
            proposal(target) = x(target);
            continue
         end
         support = [before; target; after];
         seasons = icemodel.forcing.reconstruct.seasonOf(times(support));
         if any(seasons ~= seasons(1)) ...
               || any(regime(support) ~= regime(target(1))) ...
               || any(regime(target) == 1)
            proposal(target) = x(target);
            continue
         end
         candidate = interp1([before; after], ...
            [proposal(before); proposal(after)], target);
         if ~all(icemodel.forcing.reconstruct.scalarValidity( ...
               "swd", candidate))
            proposal(target) = x(target);
            continue
         end

         proposal(target) = candidate;
         repaired_this_pass(target) = true;
         segment = false(size(x));
         segment(target) = true;
         run_audits{g} = icemodel.forcing.reconstruct.auditSegments( ...
            times, segment, "swd", "bounded_interp", sprintf( ...
            ['post-final seam QA pass %d; one-posting synthetic-side ' ...
            'flux-linear repair; observed-step percentile %.3g'], ...
            pass, kwargs.percentile));
      end
      pass_audits{pass} = vertcat(cell(0, 1), run_audits{:});
      if ~any(repaired_this_pass)
         break
      end
      x = proposal;
      provenance(repaired_this_pass) = codes.bounded_interp;
      changed = changed | repaired_this_pass;
      n_passes = pass;
   end
   audit = vertcat(cell(0, 1), pass_audits{:});

   [final_outliers, final_ratio, final_boundaries, ...
      ~, final_unsupported] = boundaryOutliers( ...
      times, x, provenance, elevation, regime, ...
      kwargs.percentile, kwargs.min_reference_steps, codes);
   final_max_ratio = max([0; final_ratio(final_outliers)]);
   tail_fraction = max(0, (100 - kwargs.percentile) / 100);
   expected_outliers = ceil(tail_fraction * final_boundaries);
   status = "pass";
   if final_unsupported > 0 ...
         || nnz(final_outliers) > expected_outliers
      status = "review";
   end
   delta = x(changed) - original(changed);
   rms_change = 0;
   max_change = 0;
   if ~isempty(delta)
      rms_change = sqrt(mean(delta .^ 2));
      max_change = max(abs(delta));
   end
   quality = table(status, kwargs.percentile, n_reference_steps, ...
      n_boundaries, initial_unsupported, nnz(initial_outliers), ...
      initial_max_ratio, final_boundaries, final_unsupported, ...
      nnz(final_outliers), expected_outliers, final_max_ratio, ...
      nnz(changed), rms_change, max_change, n_passes, ...
      'VariableNames', {'status', 'percentile', ...
      'observed_reference_steps', 'initial_boundaries', ...
      'initial_unsupported_boundaries', 'initial_outliers', ...
      'initial_max_ratio', 'final_boundaries', ...
      'final_unsupported_boundaries', 'final_outliers', ...
      'expected_tail_outliers', 'final_max_ratio', ...
      'changed_samples', 'rms_change_wm2', 'max_change_wm2', 'passes'});

   % Source-backed SWD and the separately validated D-44 twilight tier are
   % immutable. Seam QA may diagnose their boundaries but must never relabel
   % or replace them with a generic interpolation.
   if ~isequaln(x(protected), original(protected)) ...
         || ~isequal(provenance(protected), original_provenance(protected))
      error('icemodel:reconstruct:smoothShortwaveSeams:protectedMutation', ...
         'post-final seam repair modified protected SWD');
   end
end

function [outliers, ratio, n_boundaries, n_reference_steps, ...
      n_unsupported] = ...
      boundaryOutliers( ...
      times, x, provenance, elevation, regime, percentile, ...
      min_reference_steps, codes)
   %BOUNDARYOUTLIERS Compare method transitions with observed local steps.
   signed_steps = diff(x);
   jumps = abs(signed_steps);
   seasons = icemodel.forcing.reconstruct.seasonOf(times(1:end - 1));
   same_season = seasons == ...
      icemodel.forcing.reconstruct.seasonOf(times(2:end));
   native = nativeShortwave(provenance, codes);
   same_regime = regime(2:end) == regime(1:end - 1);
   near_regime_edge = ~same_regime;
   near_regime_edge(2:end) = near_regime_edge(2:end) ...
      | ~same_regime(1:end - 1);
   near_regime_edge(1:end - 1) = near_regime_edge(1:end - 1) ...
      | ~same_regime(2:end);
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   elevation_bin = discretize(elevation, ...
      bands.calibration_bin_edges_deg);
   same_band = elevation_bin(2:end) == elevation_bin(1:end - 1);
   observed_step = native(2:end) & native(1:end - 1) ...
      & isfinite(x(2:end)) & isfinite(x(1:end - 1)) ...
      & same_regime & same_band & same_season;
   coarse_observed_step = native(2:end) & native(1:end - 1) ...
      & isfinite(x(2:end)) & isfinite(x(1:end - 1)) ...
      & same_regime & same_season;
   n_reference_steps = nnz(observed_step);
   limit = empiricalStepLimit(jumps, seasons, observed_step, ...
      coarse_observed_step, elevation(2:end), regime(2:end), ...
      percentile, min_reference_steps);

   % A method transition on a smooth monotonic sunrise or sunset ramp is
   % not a seam. Keep the immediately adjacent same-direction step as a
   % local continuity floor, but never borrow it across a solar regime or
   % season boundary.
   eligible_step = isfinite(x(2:end)) & isfinite(x(1:end - 1)) ...
      & same_regime & same_season;
   previous_matches = eligible_step(2:end) & eligible_step(1:end - 1) ...
      & signed_steps(2:end) .* signed_steps(1:end - 1) > 0;
   previous_floor = zeros(size(jumps));
   previous_floor(2:end) = ...
      abs(signed_steps(1:end - 1)) .* previous_matches;
   next_matches = eligible_step(1:end - 1) & eligible_step(2:end) ...
      & signed_steps(1:end - 1) .* signed_steps(2:end) > 0;
   next_floor = zeros(size(jumps));
   next_floor(1:end - 1) = ...
      abs(signed_steps(2:end)) .* next_matches;
   local_floor = max(previous_floor, next_floor);
   independent_previous = previous_matches ...
      & native(1:end - 2) & native(2:end - 1);
   independent_previous_floor = zeros(size(jumps));
   independent_previous_floor(2:end) = ...
      abs(signed_steps(1:end - 1)) .* independent_previous;
   independent_next = next_matches ...
      & native(2:end - 1) & native(3:end);
   independent_next_floor = zeros(size(jumps));
   independent_next_floor(1:end - 1) = ...
      abs(signed_steps(2:end)) .* independent_next;
   independent_local_floor = max( ...
      independent_previous_floor, independent_next_floor);

   % The native-anchor change across one contiguous synthetic run imposes
   % a minimum possible linear step. Do not label that unavoidable bridge
   % as a seam defect.
   synthetic = ~native & provenance ~= codes.missing;
   bridge_floor = zeros(size(jumps));
   edges = diff([false; synthetic; false]);
   starts = find(edges == 1);
   stops = find(edges == -1) - 1;
   for g = 1:numel(starts)
      before = starts(g) - 1;
      after = stops(g) + 1;
      if before < 1 || after > numel(x) ...
            || ~isfinite(x(before)) || ~isfinite(x(after))
         continue
      end
      support = (before:after).';
      support_seasons = ...
         icemodel.forcing.reconstruct.seasonOf(times(support));
      if all(regime(support) == regime(support(1))) ...
            && all(support_seasons == support_seasons(1))
         bridge = abs(x(after) - x(before)) / (after - before);
         bridge_floor(before:after - 1) = ...
            max(bridge_floor(before:after - 1), bridge);
         limit(before:after - 1) = max(limit(before:after - 1), bridge);
      end
   end

   % A bounded-interpolation label inside one otherwise continuous source
   % marks a provenance transition, not a flux seam. This holds when the
   % reconstructed postings equal the exact line between finite anchors of
   % the same provenance.
   interpolation_floor = zeros(size(jumps));
   interpolated = provenance == codes.bounded_interp;
   edges = diff([false; interpolated; false]);
   starts = find(edges == 1);
   stops = find(edges == -1) - 1;
   for g = 1:numel(starts)
      before = starts(g) - 1;
      after = stops(g) + 1;
      if before < 1 || after > numel(x) ...
            || ~isfinite(x(before)) || ~isfinite(x(after))
         continue
      end
      support = (before:after).';
      expected = interp1([before; after], ...
         [x(before); x(after)], support);
      if all(abs(x(support) - expected) <= 1e-12)
         bridge = abs(x(after) - x(before)) / (after - before);
         interpolation_floor(before:after - 1) = ...
            max(interpolation_floor(before:after - 1), bridge);
      end
   end

   % A change across dark/twilight/daylight regimes is driven by the Sun,
   % not by a fill-method seam. Diagnose method boundaries only when both
   % postings occupy the same whole-interval solar regime.
   boundary = provenance(2:end) ~= provenance(1:end - 1) ...
      & isfinite(x(2:end)) & isfinite(x(1:end - 1)) ...
      & (~native(2:end) | ~native(1:end - 1)) ...
      & same_regime & same_band & same_season & ~near_regime_edge;
   % Exact continuity, an adjacent native-only slope, and the unavoidable
   % slope between two native anchors are independent support even when
   % the narrow season/elevation stratum misses the reference minimum.
   support_floor = max( ...
      max(independent_local_floor, bridge_floor), interpolation_floor);
   n_unsupported = nnz(boundary ...
      & jumps > support_floor + 1e-12 & ~isfinite(limit));
   ratio = jumps ./ max(limit, eps);
   outliers = boundary & jumps > limit + 1e-12 ...
      & jumps > local_floor + 1e-12;
   n_boundaries = nnz(boundary);
end

function native = nativeShortwave(provenance, codes)
   %NATIVESHORTWAVE Identify every source-backed PROMICE SWD provenance.
   native = ismember(provenance, [ ...
      codes.observed, codes.raw_shortwave, codes.clamped_shortwave]);
end

function limit = empiricalStepLimit(jumps, seasons, observed_step, ...
      coarse_observed_step, elevation, regimes, percentile, ...
      min_reference_steps)
   %EMPIRICALSTEPLIMIT Observed step quantile by season and solar stratum.
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   bin = discretize(elevation, bands.calibration_bin_edges_deg);
   limit = zeros(size(jumps));
   for season = ["DJF", "MAM", "JJA", "SON"]
      in_season = seasons == season;
      limit(in_season) = Inf;
      for b = unique(bin(isfinite(bin))).'
         for r = unique(regimes(in_season & bin == b)).'
            target = in_season & bin == b & regimes == r;
            reference = jumps(observed_step & target);
            if numel(reference) >= min_reference_steps
               limit(target) = prctile(reference, percentile);
               continue
            end
            % Sparse fine bins fall back only to native steps in the same
            % season and coarse dark/twilight/daylight regime.
            reference = jumps(coarse_observed_step ...
               & in_season & regimes == r);
            if numel(reference) >= min_reference_steps
               limit(target) = prctile(reference, percentile);
            end
         end
      end
   end
end
