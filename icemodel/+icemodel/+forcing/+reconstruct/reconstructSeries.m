function result = reconstructSeries(series, channel_methods, kwargs)
   %RECONSTRUCTSERIES Compose admitted fill methods into one target series.
   %
   %  result = icemodel.forcing.reconstruct.reconstructSeries( ...
   %     series, channel_methods, latitude=67.0, longitude=-48.8)
   %
   % Role
   %  The engine orchestrator (DesignSpec Resolution 6; POLICY B1 tier
   %  order). For swd a darkness pre-pass zero-fills missing
   %  below-civil-twilight samples first (a known value per POLICY B2
   %  refined by D-28, stamped with its own provenance code) so long
   %  outages decompose into daylight fragments; twilight-band samples
   %  carry real diffuse light and are left to the fill tiers. Tier 1 (bounded interior
   %  interpolation, CSI-space for swd) applies next to the
   %  interpolation-eligible channels; every
   %  remaining missing run then walks the caller's ordered, ADMITTED
   %  method estimates — donor transfers, calibrated proxies, climatology,
   %  constants — taking the first method admitted for the run's stratum
   %  (season x duration bucket) whose estimate passes the physical bounds
   %  check. Anchored boundary mismatches beyond the jump limit are
   %  BLENDED per POLICY B6: the excess offset tapers to zero across
   %  the blend window, landing the seam at half the limit, while seams
   %  already inside the limit stay untouched so skilled estimates keep
   %  their exact values. Refusal is reserved for bounds violations and
   %  missing estimates. Native finite samples are never modified; runs
   %  no admitted method fills stay missing with a per-method refusal
   %  reason in the audit. Fitting and admission happen UPSTREAM (harness
   %  experiments + fit/apply helpers); this function only composes and
   %  stamps, which is what keeps it family-generic under the role
   %  contract.
   %
   % Inputs
   %  series : target timetable (native, gaps preserved).
   %  channel_methods : struct array, one element per channel to fill —
   %     channel : canonical channel name present in the series.
   %     methods : ordered struct array of admitted estimates —
   %        name : method label for audit/registry (e.g. "donor:aws10").
   %        code : uint8 provenance code from provenanceCodes.
   %        estimate : Nx1 estimate on the series axis (NaN = declines).
   %        seasons : admitted seasons (string vector; "all" for every).
   %        buckets : admitted duration buckets (double vector; [] = all).
   %        max_validated_hours : longest held-out gap tested for each
   %           admitted bucket (optional only for direct/manual callers).
   %
   % Name-value
   %  latitude, longitude : target site point (swd CSI tier and census).
   %  interp_channels : channels eligible for tier 1 (swu is derived and
   %     never interpolated independently).
   %  cap_hours : tier-1 interior cap in wall-clock hours.
   %  cap_hours_by_channel : per-channel tier-1 cap overrides (POLICY
   %     B3/D-21); fields name interp channels and replace cap_hours for
   %     that channel only, within interpolationCapHours.
   %  jump_factor : boundary-jump multiplier on the seasonal step scale,
   %     used by the tier-1 interpolation check and the post-blend
   %     safety verify.
   %  blend_hours : seam-blend taper window; an anchored seam mismatch
   %     beyond the jump limit is corrected by an offset decaying
   %     linearly to zero across this many hours.
   %  toa_dark_wm2 : irradiance threshold below which a sample counts as
   %     dark for the tier-1 CSI mask (fillShortGaps). The darkness
   %     zero-fill itself keys on the civil-twilight solar elevation from
   %     solarElevationBands (D-28), not on this parameter, so twilight
   %     samples with real diffuse light are never forced to zero.
   %  max_validation_duration_factor : maximum real-gap duration divided
   %     by the longest held-out gap that admitted the method.
   %  native_provenance : optional struct of source-backed uint8 code
   %     vectors keyed by channel; unspecified native values are observed.
   %  Defaults for all of the above come from the central
   %  icemodel.forcing.reconstruct.setopts contract.
   %
   % Returns
   %  result : struct with fields
   %     series : the filled timetable (same shape as the input).
   %     provenance : timetable of uint8 codes, one channel per filled
   %        channel, on the same time axis.
   %     audit : segment-audit table (channel, start_time, end_time,
   %        duration_hours, method, detail, context_id). A production
   %        context_id joins to the persisted station plan's exact fitted
   %        parameters and held-out metrics.
   %     registry : provenanceCodes() struct recorded with the product.
   %
   % See also: icemodel.forcing.reconstruct.fillShortGaps,
   %  icemodel.forcing.reconstruct.applyDonorTransfer,
   %  icemodel.forcing.reconstruct.applyProxyCalibration

   arguments
      series timetable
      channel_methods (1, :) struct
      kwargs.latitude (1, 1) double = NaN
      kwargs.longitude (1, 1) double = NaN
      kwargs.interp_channels (1, :) string = ...
         icemodel.forcing.reconstruct.setopts().interp_channels
      kwargs.cap_hours (1, 1) double {mustBePositive, ...
         icemodel.forcing.reconstruct.mustBeCapHours(kwargs.cap_hours)} = ...
         icemodel.forcing.reconstruct.setopts().cap_hours
      % Per-channel cap overrides (POLICY B3/D-21); fields replace
      % cap_hours for that channel only.
      kwargs.cap_hours_by_channel (1, 1) struct = ...
         icemodel.forcing.reconstruct.setopts().cap_hours_by_channel
      kwargs.jump_factor (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().jump_factor
      kwargs.blend_hours (1, 1) double {mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().blend_hours
      kwargs.toa_dark_wm2 (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().toa_dark_wm2
      kwargs.max_validation_duration_factor (1, 1) double ...
         {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().max_validation_duration_factor
      kwargs.native_provenance (1, 1) struct = struct()
   end

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   times = series.Properties.RowTimes;
   dt_hours = hours(median(diff(times)));
   % The swu relational bound depends on the current swd channel. Process
   % swu last so caller ordering cannot make it validate against stale NaNs.
   method_channels = string({channel_methods.channel});
   channel_methods = channel_methods( ...
      [find(method_channels ~= "swu"), find(method_channels == "swu")]);
   provenance = timetable(times);
   % One audit-cell slot per channel, sized inside the loop; assembled once
   % at the end so nothing grows sample-by-sample.
   channel_audits = cell(numel(channel_methods), 1);

   for c = 1:numel(channel_methods)
      channel = string(channel_methods(c).channel);
      if ~ismember(channel, string(series.Properties.VariableNames))
         error('icemodel:reconstruct:reconstructSeries:unknownChannel', ...
            'channel %s is not in the series', channel);
      end
      x = series.(channel);
      native_finite = isfinite(x);
      code = repmat(codes.missing, numel(x), 1);
      code(native_finite) = codes.observed;
      if isfield(kwargs.native_provenance, channel)
         native_code = kwargs.native_provenance.(channel);
         if ~isa(native_code, 'uint8') || ~iscolumn(native_code) ...
               || numel(native_code) ~= numel(x)
            error(['icemodel:reconstruct:reconstructSeries:' ...
               'invalidNativeProvenance'], ...
               'native provenance for %s must be an axis-length uint8 column', ...
               channel);
         end
         code(native_finite) = native_code(native_finite);
      end
      % Freeze the observed step scale before darkness, interpolation, or
      % any later method can contribute synthetic values to seam limits.
      scale = icemodel.forcing.reconstruct.stepScale(times, x);

      % Policy swd darkness rule (B2, refined by D-28): a missing sample
      % with the sun below civil twilight is a KNOWN zero, not an
      % estimate. Zero-filling deep darkness FIRST decomposes a multi-day
      % outage into per-day daylight fragments the admitted duration
      % buckets cover — validation strata are daylight-bounded, so
      % buckets above ~one day can never admit a shortwave method.
      % Samples in the twilight band (civil twilight up to 0 deg) carry
      % real diffuse light stations measure at 15-38 W m^-2 (D-28), so
      % they stay missing here and reach the normal fill tiers instead of
      % being forced to zero.
      % Hoisted per-channel solar geometry: the run x method walk below
      % validates every swd candidate against the TOA ceiling and the
      % whole-posting twilight allowance. One interval-maximum evaluation
      % on the complete axis is reused by darkness and sliced per run.
      axis_elevation = zeros(0, 1);
      axis_toa = zeros(0, 1);
      if channel == "swd" && isfinite(kwargs.latitude) ...
            && isfinite(kwargs.longitude)
         axis_elevation = ...
            icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
            times, kwargs.latitude, kwargs.longitude, hours(dt_hours));
         axis_toa = icemodel.forcing.reconstruct.toaIrradiance(times, ...
            kwargs.latitude, kwargs.longitude);
      end
      audit_dark = cell(0, 1);
      if channel == "swd" ...
            && isfinite(kwargs.latitude) && isfinite(kwargs.longitude)
         % Darkness uses the maximum signed elevation over the complete
         % interval support, matching PROMICE staging: a posting is known
         % dark only when its start, quarter points, and end all remain at
         % or below civil twilight.
         bands = icemodel.forcing.reconstruct.solarElevationBands();
         dark = ~isfinite(x) ...
            & axis_elevation <= bands.civil_twilight_deg;
         if any(dark)
            x(dark) = 0;
            code(dark) = codes.darkness;
            audit_dark = icemodel.forcing.reconstruct.auditSegments( ...
               times, dark, channel, "darkness_zero", ...
               sprintf('interval maximum solar elevation <= %g deg', ...
               bands.civil_twilight_deg));
         end
      end

      % Tier 1: bounded interior interpolation for eligible channels,
      % honoring any per-channel cap override (POLICY B3/D-21).
      channel_cap = kwargs.cap_hours;
      if isfield(kwargs.cap_hours_by_channel, channel)
         channel_cap = kwargs.cap_hours_by_channel.(channel);
         % Direct engine callers bypass setopts, so the override honors
         % the same B3 ceiling contract here as well.
         icemodel.forcing.reconstruct.mustBeCapHours(channel_cap);
      end
      audit1 = cell(0, 1);
      if ismember(channel, kwargs.interp_channels)
         [x, filled1, audit1] = ...
            icemodel.forcing.reconstruct.fillShortGaps(times, x, ...
            channel, cap_hours=channel_cap, ...
             latitude=kwargs.latitude, longitude=kwargs.longitude, ...
             jump_factor=kwargs.jump_factor, ...
             blend_hours=kwargs.blend_hours, ...
             toa_dark_wm2=kwargs.toa_dark_wm2, step_scale=scale);
         code(filled1) = codes.bounded_interp;
      end

      % Remaining missing runs are split at season boundaries so every
      % production segment matches one held-out admission stratum.
      miss = ~isfinite(x);
      seasons = icemodel.forcing.reconstruct.seasonOf(times);
      starts = find(miss & [true; ~miss(1:end - 1) ...
         | seasons(2:end) ~= seasons(1:end - 1)]);
      stops = find(miss & [~miss(2:end) ...
         | seasons(1:end - 1) ~= seasons(2:end); true]);
      methods = channel_methods(c).methods;
       % Method attempts can fill disjoint fragments, so audit row groups
       % are collected and flattened after the run loop.
       run_audits = cell(0, 1);
      blend_len = max(1, round(kwargs.blend_hours / dt_hours));
      for g = 1:numel(starts)
         idx = (starts(g):stops(g)).';
         run_len = numel(idx);
         duration = run_len * dt_hours;
         bucket = icemodel.forcing.reconstruct.gapDurationBucket(duration);
         season = seasons(starts(g));
         % The run's boundary anchors are fixed OUTSIDE the run before any
         % method fills it; deriving them from a shrunken leftover index
         % would test seams against another method's fresh fill.
         run_before = starts(g) - 1;
         run_after = stops(g) + 1;
         if run_before >= 1 && miss(run_before)
            run_before = 0;
         end
         if run_after <= numel(x) && miss(run_after)
            run_after = numel(x) + 1;
         end
         % Per-run seam limit: floored at the linear bridge's implied
         % per-step change between the two anchors — when the native
         % record itself steps across the gap by more than the seasonal
         % limit allows (a dawn shortwave transition, a frontal jump),
         % no fill value can satisfy both seams and refusal would be
         % unphysical. A fill at least as smooth as the straight bridge
         % is never refused; the tiny slack keeps float-exact midpoints
         % on the passing side of the strict > test.
         run_limit = kwargs.jump_factor * scale.(char(season));
         if run_before >= 1 && run_after <= numel(x) ...
               && isfinite(x(run_before)) && isfinite(x(run_after))
            bridge_step = abs(x(run_after) - x(run_before)) ...
               / (run_len + 1);
            run_limit = max(run_limit, bridge_step * (1 + 1e-6));
         end
         % Refusal bookkeeping: every attempt that declines the run (or
         % part of it) records why, so residual missing runs are
         % explained in the audit instead of silently omitted.
         reasons = strings(numel(methods), 1);
         n_reason = 0;
         for m = 1:numel(methods)
            method = methods(m);
            [admitted, admission_reason] = stratumAdmitted(method, ...
               season, bucket, duration, ...
               kwargs.max_validation_duration_factor);
            if ~admitted
               n_reason = n_reason + 1;
               reasons(n_reason) = method.name + ": " + admission_reason;
               continue
            end
            candidate = method.estimate(idx);
            have = isfinite(candidate);
            if ~any(have)
               n_reason = n_reason + 1;
               reasons(n_reason) = method.name + ": no finite estimate";
               continue
            end
             % The shared validator enforces both scalar and relational
             % limits; invalid samples cascade to the next admitted method.
             if channel == "swu"
                physically_valid = ...
                   icemodel.forcing.reconstruct.physicalValidity( ...
                   channel, candidate, times(idx), swd=series.swd(idx));
             elseif channel == "swd" && ~isempty(axis_toa)
                % Slices of the hoisted axis geometry keep the validator
                % bit-identical while skipping its per-call NOAA passes.
                physically_valid = ...
                   icemodel.forcing.reconstruct.physicalValidity( ...
                   channel, candidate, times(idx), ...
                   latitude=kwargs.latitude, longitude=kwargs.longitude, ...
                   toa=axis_toa(idx), elevation=axis_elevation(idx));
             else
                physically_valid = ...
                   icemodel.forcing.reconstruct.physicalValidity( ...
                   channel, candidate, times(idx), ...
                   latitude=kwargs.latitude, longitude=kwargs.longitude);
             end
             usable = have & physically_valid;
            if ~any(usable)
               n_reason = n_reason + 1;
               reasons(n_reason) = method.name + ": all out of bounds";
               continue
            end
            % POLICY B6: anchored seam offsets that exceed the jump
            % limit blend away across the taper window rather than
            % refusing the run — admitted-method point RMSE routinely
            % exceeds the seasonal step scale, so an exact-seam demand
            % would reject skilled fills wholesale. Seams already inside
            % the limit are left untouched, so skilled estimates keep
            % their exact values.
            % A prior method can leave disjoint residual fragments. Blend
            % each contiguous wall-clock segment independently so ordinal
            % position never carries a taper across an already filled span.
            segment_starts = [1; find(diff(idx) > 1) + 1];
            segment_stops = [segment_starts(2:end) - 1; numel(idx)];
            seam_note = "";
            for s = 1:numel(segment_starts)
               part = (segment_starts(s):segment_stops(s)).';
               [candidate(part), part_note] = ...
                  icemodel.forcing.reconstruct.blendSeams( ...
                  x, candidate(part), idx(part), usable(part), ...
                  starts(g), stops(g), run_before, run_after, ...
                  blend_len, run_limit);
               seam_note = seam_note + string(part_note);
            end
            seam_note = char(seam_note);
             % The taper can violate a scalar or relational bound; drop
             % those samples rather than clip so limits stay hard.
             if channel == "swu"
                physically_valid = ...
                   icemodel.forcing.reconstruct.physicalValidity( ...
                   channel, candidate, times(idx), swd=series.swd(idx));
             elseif channel == "swd" && ~isempty(axis_toa)
                % Same hoisted-geometry slices as the pre-blend check.
                physically_valid = ...
                   icemodel.forcing.reconstruct.physicalValidity( ...
                   channel, candidate, times(idx), ...
                   latitude=kwargs.latitude, longitude=kwargs.longitude, ...
                   toa=axis_toa(idx), elevation=axis_elevation(idx));
             else
                physically_valid = ...
                   icemodel.forcing.reconstruct.physicalValidity( ...
                   channel, candidate, times(idx), ...
                   latitude=kwargs.latitude, longitude=kwargs.longitude);
             end
             usable = usable & physically_valid;
            if ~any(usable)
               n_reason = n_reason + 1;
               reasons(n_reason) = method.name + ": out of bounds after blend";
               continue
            end
            % Safety verify: anchored seams match by construction after
            % blending, so a violation here flags an engine regression,
            % not a data decision.
            if jumpViolation(x, candidate, idx, usable, ...
                  starts(g), stops(g), run_before, run_after, run_limit)
               n_reason = n_reason + 1;
               reasons(n_reason) = method.name + ": seam jump after blend";
               continue
            end
            x(idx(usable)) = candidate(usable);
            code(idx(usable)) = method.code;
             filled_idx = idx(usable);
             filled_mask = false(numel(times), 1);
             filled_mask(filled_idx) = true;
             rows = icemodel.forcing.reconstruct.auditSegments( ...
                times, filled_mask, channel, method.name, sprintf( ...
                'season %s bucket %d, filled %d/%d run samples%s', ...
                season, bucket, numel(filled_idx), run_len, seam_note), ...
                context_id=methodContextId(method));
             run_audits = [run_audits; rows]; %#ok<AGROW>
            if all(usable)
               break
            end
            % Leftover samples advance to later methods; the run anchors
            % captured above keep the seam checks honest.
            idx = idx(~usable);
            if isempty(idx)
               break
            end
         end
         % Residual-missing runs get one audit row joining every refusal
         % reason, honoring the documented contract.
          if any(~isfinite(x(starts(g):stops(g))))
             residual = false(numel(times), 1);
             residual(starts(g):stops(g)) = ...
                ~isfinite(x(starts(g):stops(g)));
             rows = icemodel.forcing.reconstruct.auditSegments( ...
                times, residual, channel, "unfilled", ...
                char("residual missing; " + ...
                strjoin(reasons(1:n_reason), "; ")));
             run_audits = [run_audits; rows]; %#ok<AGROW>
         end
      end

      series.(channel) = x;
      provenance.(channel) = code;
       channel_audits{c} = [audit_dark; audit1; run_audits];
   end

   audit_rows = vertcat(channel_audits{:});
    if isempty(audit_rows)
       audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
          'datetime', 'datetime', 'double', 'cellstr', 'cellstr', ...
          'cellstr'}, ...
          'VariableNames', {'channel', 'start_time', 'end_time', ...
          'duration_hours', 'method', 'detail', 'context_id'});
       % Empty datetime columns must retain the target axis timezone so a
       % later last-resort or constant row can concatenate safely.
       audit.start_time.TimeZone = times.TimeZone;
       audit.end_time.TimeZone = times.TimeZone;
   else
      audit = cell2table(vertcat(audit_rows{:}), 'VariableNames', ...
         {'channel', 'start_time', 'end_time', 'duration_hours', ...
         'method', 'detail', 'context_id'});
   end
   result = struct('series', series, 'provenance', provenance, ...
      'audit', audit, 'registry', codes);
end

function context_id = methodContextId(method)
   %METHODCONTEXTID Link a segment to its persisted fitted-method record.
   context_id = string(method.name);
   if isfield(method, 'audit_context_id') ...
         && strlength(string(method.audit_context_id)) > 0
      context_id = string(method.audit_context_id);
   end
end

function [tf, reason] = stratumAdmitted(method, season, bucket, ...
      duration_hours, max_duration_factor)
   %STRATUMADMITTED True when the method is admitted for this stratum.
   seasons = string(method.seasons);
   tf = any(seasons == "all") || any(seasons == season);
   reason = sprintf("not admitted (%s b%d)", season, bucket);
   if tf && ~isempty(method.buckets)
      % Every duration bucket has its own measured admission. Bucket 1
      % reaches this walk only when tier 1 declines, and its candidates
      % must have beaten persistence on held-out <=6 h gaps.
      graded_bucket = bucket;
      tf = ismember(graded_bucket, method.buckets);
      if tf && isfield(method, 'max_validated_hours') ...
            && ~isempty(method.max_validated_hours)
         k = find(method.buckets == graded_bucket, 1);
         max_hours = max_duration_factor * method.max_validated_hours(k);
         if isfinite(max_hours) && duration_hours > max_hours
            tf = false;
            reason = sprintf( ...
               "gap %.3g h exceeds %.3g h held-out duration limit", ...
               duration_hours, max_hours);
         end
      end
   end
end

function tf = jumpViolation(x, candidate, idx, usable, run_start, ...
      run_stop, run_before, run_after, limit)
   %JUMPVIOLATION POLICY B6 check at the run's fixed boundary anchors.
   % Anchors sit outside the original run (native or tier-1 values that no
   % method in this loop modifies). Each seam is tested only when this
   % method actually fills the run-edge sample; interior seams between two
   % methods' partial fills carry no observed anchor to test against.
   tf = false;
   first_use = find(usable, 1, 'first');
   last_use = find(usable, 1, 'last');
   if idx(first_use) == run_start && run_before >= 1 ...
         && isfinite(x(run_before))
      tf = tf || abs(candidate(first_use) - x(run_before)) > limit;
   end
   if idx(last_use) == run_stop && run_after <= numel(x) ...
         && isfinite(x(run_after))
      tf = tf || abs(candidate(last_use) - x(run_after)) > limit;
   end
end
