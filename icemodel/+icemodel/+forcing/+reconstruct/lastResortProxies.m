function [filled, provenance, audit, denials] = lastResortProxies(filled, ...
      provenance, audit, proxies, codes, opts, kwargs)
   %LASTRESORTPROXIES Adopt aligned proxy values for residual gaps.
   %
   %  [filled, provenance, audit, denials] = ...
   %     icemodel.forcing.reconstruct.lastResortProxies(filled, ...
   %     provenance, audit, proxies, codes, opts)
   %
   % Role
   %  The policy's final tier: required-channel samples still missing
    %  after method composition take one proxy source per whole outage.
    %  The first catalog source that covers the whole outage wins; if none
    %  does, the first source with any usable values supplies a partial fill.
    %  A later source never fills leftovers inside the same outage, so one
   %  span keeps thermodynamically coupled channels consistent with each
   %  other, which per-channel method composition cannot guarantee.
   %  Native and method-filled samples are never overwritten; each
    %  adoption appends one audit row per contiguous channel/source segment.
   %
   % Inputs
    %  filled : composed timetable (post reconstructSeries).
   %  provenance : matching per-channel uint8 code timetable.
   %  audit : segment-audit table to append adoption rows to.
    %  proxies : proxy struct array (series, name, code_name) in
    %     adoption-preference order. A source/channel with a persisted
    %     overlap calibration adopts corrected values; without one the raw
    %     values adopt identity and are stamped low-confidence in the
    %     audit (POLICY A11/D-25).
   %  codes : provenanceCodes() registry.
    %  opts : reconstruction options (required_channels).
    %  latitude, longitude : station geometry for the shortwave bound.
    %  native : pre-reconstruction timetable used only for seam step scales.
    %  plan : stationMethodPlan output containing persisted proxy corrections.
   %
   % Returns
   %  filled, provenance, audit : the adopted product, codes, and rows.
   %  denials : struct of per-sample final-tier denial notes, one string
   %     array per handled channel ("" = no denial). The driver appends
   %     these to residual unfilled audit rows so the report states the
   %     ACTUAL last-tier cause (no usable source, source without valid
   %     values, post-blend rejection) instead of a stale provisional
   %     reason from an earlier tier.
   %
   % See also: icemodel.forcing.reconstruct.fillPromiceStation,
   %  icemodel.forcing.reconstruct.setopts

   arguments
      filled timetable
      provenance timetable
      audit table
      proxies (1, :) struct
      codes (1, 1) struct
      opts (1, 1) struct
      kwargs.latitude (1, 1) double = NaN
      kwargs.longitude (1, 1) double = NaN
      kwargs.native timetable = timetable()
      kwargs.plan (1, 1) struct = struct()
   end

   times = filled.Properties.RowTimes;

   % Precipitation channels have a dedicated adoption policy. The code
   % derives swu from the final albedo and swd, so this tier does not
   % adopt swu on its own.
   channels = setdiff(opts.required_channels, ...
      ["swu", icemodel.forcing.helpers.precipitationVariables()], 'stable');
   channels = intersect(channels, ...
      string(filled.Properties.VariableNames), 'stable');

   % Per-sample final-tier denial notes ("" = none) for every handled
   % channel. Audit reconciliation appends these notes to the residual
   % unfilled rows, so the reported reason is the last tier's refusal.
   denials = struct();
   for name = channels
      denials.(name) = strings(height(filled), 1);
   end

   % The code selects one source per contiguous span in which any
   % required channel is still missing. One source per span keeps the
   % coupled channels consistent.
   missing = false(height(filled), 1);
   for name = channels
      missing = missing | ~isfinite(filled.(name));
   end
   edge = diff([false; missing; false]);
   starts = find(edge == 1);
   stops = find(edge == -1) - 1;
   row_blocks = cell(numel(starts) * numel(channels), 1);
   n_blocks = 0;

   % Compute the swd solar geometry once for the whole axis, because
   % proxyPlan and the seam fix-up both validate per outage. The maximum
   % solar elevation over the whole posting keeps those decisions aligned
   % with staging and with the darkness rule.
   axis_toa = zeros(0, 1);
   axis_elevation = zeros(0, 1);
   if ismember("swd", channels) && isfinite(kwargs.latitude) ...
         && isfinite(kwargs.longitude)
      interval = median(diff(times));
      axis_elevation = ...
         icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
         times, kwargs.latitude, kwargs.longitude, interval);
      axis_toa = icemodel.forcing.reconstruct.toaIrradiance(times, ...
         kwargs.latitude, kwargs.longitude);
   end

   % Cache of the aligned, calibrated source series for each (proxy, channel)
   % pair. calibratedChannel depends only on the source and the channel, never
   % on the outage. The first probe computes the series. Every later outage
   % reuses the cached copy, and does not align and calibrate the full axis
   % again.
   source_cache = cell(numel(proxies), numel(channels));
   for g = 1:numel(starts)
      idx = (starts(g):stops(g)).';
      chosen = 0;
      chosen_values = cell(numel(channels), 1);
      chosen_masks = cell(numel(channels), 1);
      chosen_calibrated = true(numel(channels), 1);
      chosen_clamped = false(numel(channels), 1);
      fallback = 0;
      fallback_values = cell(numel(channels), 1);
      fallback_masks = cell(numel(channels), 1);
      fallback_calibrated = true(numel(channels), 1);
      fallback_clamped = false(numel(channels), 1);
      for p = 1:numel(proxies)
         [values, masks, n_usable, complete, calibrated, clamped, ...
            source_cache] = proxyPlan(proxies(p), channels, filled, ...
            times, idx, kwargs, p, source_cache, axis_toa, ...
            axis_elevation);
         if n_usable > 0 && fallback == 0
            fallback = p;
            fallback_values = values;
            fallback_masks = masks;
            fallback_calibrated = calibrated;
            fallback_clamped = clamped;
         end
         if complete
            chosen = p;
            chosen_values = values;
            chosen_masks = masks;
            chosen_calibrated = calibrated;
            chosen_clamped = clamped;
            break
         end
      end
      if chosen == 0
         chosen = fallback;
         chosen_values = fallback_values;
         chosen_masks = fallback_masks;
         chosen_calibrated = fallback_calibrated;
         chosen_clamped = fallback_clamped;
      end
      if chosen == 0
         % No catalog source has even one usable value for this span.
         % Record that final-tier cause for every channel sample that the
         % span leaves missing, so the audit reports the real refusal.
         for c = 1:numel(channels)
            name = channels(c);
            needs = idx(~isfinite(filled.(name)(idx)));
            denials.(name)(needs) = "no usable last-resort source";
         end
         continue
      end

      % Apply only the selected source and audit every contiguous segment.
      for c = 1:numel(channels)
         name = channels(c);
         adopt = chosen_masks{c};

         % Samples that the chosen source cannot supply stay missing under
         % the whole-outage policy, because a later source never fills the
         % leftovers. Record that final-tier cause for the audit.
         needs = ~isfinite(filled.(name)(idx));
         leftover = idx(needs & ~adopt);
         denials.(name)(leftover) = ...
            "last-resort source lacks valid values (whole-outage policy)";

         % Lone-island guard: one admitted posting with missing values on
         % BOTH sides becomes an isolated spike pair when the hourly
         % posting expands over its quarter-hour support (the KAN_M
         % 2011-03-22 class: one bright MAR hour between unfillable
         % samples). One posting between two holes carries no usable
         % structure, so the code refuses such runs and records them as
         % final-tier denials. A record edge counts as missing, because no
         % support exists beyond it.
         adopt_global = false(numel(times), 1);
         adopt_global(idx(adopt)) = true;
         run_edge = diff([false; adopt_global; false]);
         run_starts = find(run_edge == 1);
         run_stops = find(run_edge == -1) - 1;
         x_now = filled.(name);
         for r = 1:numel(run_starts)
            if run_stops(r) - run_starts(r) + 1 >= 2
               continue
            end
            sample = run_starts(r);

            % Runs are maximal, so this pass never adopts the neighbors.
            % Only their current (pre-adoption) state counts.
            left_missing = sample == 1 ...
               || ~isfinite(x_now(sample - 1));
            right_missing = sample == numel(times) ...
               || ~isfinite(x_now(sample + 1));
            if left_missing && right_missing
               adopt_global(sample) = false;
               denials.(name)(sample) = ...
                  "lone single-posting adoption between missing neighbors";
            end
         end
         adopt = adopt_global(idx);
         if ~any(adopt)
            continue
         end
         target = idx(adopt);
         x = filled.(name);
         candidate = x;
         candidate(target) = chosen_values{c}(adopt);
         unblended = candidate;
         adopted = false(numel(times), 1);
         adopted(target) = true;
         native = x;
         if ismember(name, string(kwargs.native.Properties.VariableNames))
            native = kwargs.native.(name);
         end
         [candidate, seam_note] = ...
            icemodel.forcing.reconstruct.blendFallbackSeams( ...
            times, native, x, candidate, adopted, ...
            jump_factor=opts.jump_factor, blend_hours=opts.blend_hours);

         % Blending can move a pre-validated estimate outside a scalar
         % or relational candidate bound. A15/B6 makes that post-blend
         % value inadmissible. Leave the sample missing. A clamp is not
         % approved here.
         if name == "swd" && ~isempty(axis_toa)

            % Slices of the precomputed axis geometry give bit-identical
            % validator results and skip the per-outage NOAA passes.
            valid = ...
               icemodel.forcing.reconstruct.physicalValidity( ...
               name, candidate(target), times(target), ...
               latitude=kwargs.latitude, longitude=kwargs.longitude, ...
               toa=axis_toa(target), elevation=axis_elevation(target));
         else
            valid = ...
               icemodel.forcing.reconstruct.physicalValidity( ...
               name, candidate(target), times(target), ...
               latitude=kwargs.latitude, longitude=kwargs.longitude);
         end

         % Wind speed must stay above its floor at runtime. Optional seam
         % optimization can push an already-valid source posting outside
         % that window. Retain the whole unblended source segment in that
         % case, so a calm clamp or a residual outage does not replace the
         % source variability. Other channels keep the general post-blend
         % refusal rule, because their bounds are physical. The channel
         % list is policy, so it comes from setopts.
         if any(name == opts.seam_retain_unblended_channels) && any(~valid)
            candidate(target) = unblended(target);
            valid = ...
               icemodel.forcing.reconstruct.physicalValidity( ...
               name, candidate(target), times(target), ...
               latitude=kwargs.latitude, longitude=kwargs.longitude);
            seam_note = string(seam_note) + ...
               ", seam blend crossed wind bounds; " + ...
               "prevalidated wind source retained";
         end
         denials.(name)(target(~valid)) = "post-blend seam rejection";
         target = target(valid);
         adopted(:) = false;
         adopted(target) = true;
         if isempty(target)
            continue
         end
         x(target) = candidate(target);
         filled.(name) = x;
         code = provenance.(name);
         code(target) = codes.(char(proxies(chosen).code_name));
         provenance.(name) = code;

         % Identity adoptions carry a low-confidence stamp. The policy accepts
         % them so the model can run where no overlap calibration exists.
         detail = "whole-outage source";
         if ~chosen_calibrated(c)
            detail = detail + ...
               "; identity adoption, no overlap calibration " + ...
               "(low confidence)";
         end

         % D-27: the code clamps corrected rh candidates back into their
         % bounds. The audit records that the clamp occurred.
         if chosen_clamped(c)
            detail = detail + "; calibration clamped to bounds";
         end
         segment_rows = icemodel.forcing.reconstruct.auditSegments( ...
            times, adopted, name, ...
            "proxy:" + proxies(chosen).name + ":last_resort", ...
            detail + string(seam_note));
         n_blocks = n_blocks + 1;
         row_blocks{n_blocks} = segment_rows;
      end
   end
   rows = vertcat(cell(0, 1), row_blocks{1:n_blocks});
   if ~isempty(rows)
      audit = [audit; cell2table(vertcat(rows{:}), ...
         'VariableNames', audit.Properties.VariableNames)];
   end
end

function [values, masks, n_usable, complete, calibrated, clamped, ...
      source_cache] = proxyPlan(proxy, channels, filled, times, idx, ...
      kwargs, proxy_index, source_cache, axis_toa, axis_elevation)
   %PROXYPLAN Test one proxy as the sole source for one whole outage.

   values = cell(numel(channels), 1);
   masks = cell(numel(channels), 1);

   % The function returns the per-channel calibration status, so identity
   % adoptions get their low-confidence audit stamp (POLICY A11/D-25). It
   % returns the clamp flag, so D-27 bound clamps get their audit note.
   calibrated = true(numel(channels), 1);
   clamped = false(numel(channels), 1);
   n_missing = 0;
   n_usable = 0;

   for c = 1:numel(channels)
      % swu never reaches this tier. The code derives swu after every
      % tier (POLICY B10), and the adoption channel set excludes it
      % (A11). Every channel here uses the geometric validity contract.
      name = channels(c);

      % Alignment and calibration depend only on (source, channel). The
      % first outage that probes this pair fills the cache, and every
      % later outage reuses the identical full-axis series.
      if isempty(source_cache{proxy_index, c})
         [source, was_calibrated, clamp_mask] = calibratedChannel( ...
            kwargs.plan, proxy, name, times, axis_elevation);
         source_cache{proxy_index, c} = struct('source', source, ...
            'calibrated', was_calibrated, 'clamp_mask', clamp_mask);
      end
      cached = source_cache{proxy_index, c};
      calibrated(c) = cached.calibrated;
      values{c} = cached.source(idx);
      needs = ~isfinite(filled.(name)(idx));
      n_missing = n_missing + nnz(needs);
      if name == "swd" && ~isempty(axis_toa)

         % Slices of the precomputed axis geometry give bit-identical
         % validator results and skip the per-outage NOAA passes.
         valid = icemodel.forcing.reconstruct.physicalValidity( ...
            name, values{c}, times(idx), ...
            latitude=kwargs.latitude, longitude=kwargs.longitude, ...
            toa=axis_toa(idx), elevation=axis_elevation(idx));
      else
         valid = icemodel.forcing.reconstruct.physicalValidity( ...
            name, values{c}, times(idx), ...
            latitude=kwargs.latitude, longitude=kwargs.longitude);
      end
      masks{c} = needs & valid;
      n_usable = n_usable + nnz(masks{c});

      % A clamp matters to the audit only when this outage adopts a
      % clamped sample.
      clamped(c) = any(cached.clamp_mask(idx) & masks{c});
   end
   complete = n_missing > 0 && n_usable == n_missing;
end

function [source, calibrated, clamp_mask] = calibratedChannel(plan, ...
      proxy, channel, times, target_elevation)
   %CALIBRATEDCHANNEL Apply the planner's stored overlap correction.
   % When no usable overlap calibration exists (the station never
   % observed the channel together with this source), the raw aligned
   % values adopt IDENTITY. If the code refused them, in-bounds gaps of
   % millions of samples would stay unfilled at stations such as DY2 and
   % SWC. The caller stamps such adoptions low-confidence in the audit
   % (POLICY A11/D-25).

   calibrated = false;
   source = alignedChannel(proxy, channel, times);

   % Per-sample record of D-27 clamping for the audit note.
   clamp_mask = false(numel(times), 1);
   if ~isfield(plan, 'channels')
      return
   end
   channel_index = find(string({plan.channels.channel}) == channel);
   if numel(channel_index) ~= 1 ...
         || ~isfield(plan.channels(channel_index), 'proxy_calibrations')
      return
   end
   calibrations = plan.channels(channel_index).proxy_calibrations;
   source_index = find(string({calibrations.source}) == string(proxy.name));
   if numel(source_index) ~= 1
      return
   end
   parameters = calibrations(source_index).parameters;
   if ~isfield(parameters, 'n_overlap') ...
         || ~isscalar(parameters.n_overlap) ...
         || ~isfinite(parameters.n_overlap) || parameters.n_overlap <= 0
      return
   end
   calibrated = true;

   % D-27: applyProxyCalibration applies the clamp, and the method tier
   % uses the same function. An overlap correction that pushes humidity
   % past its bounds is calibration ARITHMETIC, not physics. If the code
   % refused those candidates, 4.5% of rh at SWC would stay unfilled. The
   % second output reports where the clamp occurred, so the caller can
   % record it in the audit.
   if channel == "swd" && ~isempty(target_elevation)

      % D-28 elevation bins correct the shoulder bias at the last-resort
      % tier just as they do for admitted proxy methods.
      [source, clamp_mask] = ...
         icemodel.forcing.reconstruct.applyProxyCalibration( ...
         parameters, times, source, target_elevation=target_elevation);
   else
      [source, clamp_mask] = ...
         icemodel.forcing.reconstruct.applyProxyCalibration( ...
         parameters, times, source);
   end
end

function source = alignedChannel(proxy, channel, times)
   %ALIGNEDCHANNEL Align one proxy channel to the target axis.

   source = nan(numel(times), 1);
   if ~ismember(channel, string(proxy.series.Properties.VariableNames))
      return
   end
   [tf, loc] = ismember(times, proxy.series.Properties.RowTimes);
   source(tf) = proxy.series.(channel)(loc(tf));
end
