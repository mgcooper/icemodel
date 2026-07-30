function [x, filled, audit] = fillShortGaps(times, x, channel, kwargs)
   %FILLSHORTGAPS Tier-1 bounded interior interpolation of one channel.
   %
   %  [x, filled, audit] = icemodel.forcing.reconstruct.fillShortGaps( ...
   %     times, x, "tair")
   %  [x, filled, audit] = icemodel.forcing.reconstruct.fillShortGaps( ...
   %     times, x, "swd", latitude=67.0, longitude=-48.8)
   %
   % Role
   %  The approved short-gap tier (POLICY B3): linear
   %  interpolation bridges interior gaps up to the wall-clock cap for
   %  state variables and albedo, while shortwave is never raw-linearly interpolated
   %  across daylight — its short gaps interpolate in clear-sky-index
   %  space (swd / top-of-atmosphere irradiance) and reconvert, which
   %  preserves the diurnal cycle a raw line would flatten. Every candidate
   %  fill passes physical bounds (A15), and excess boundary mismatch
   %  is tapered under POLICY B6 rather than refused. Native finite
   %  samples are never modified.
   %
   % Name-value
   %  cap_hours : maximum interior gap bridged by this tier. Caps are
   %     per-channel (POLICY B3/D-21) and the caller resolves this
   %     channel's value. interpolationCapHours is the shared registry:
   %     six hours ordinarily, nine for SWD and RH, and 30 for albedo.
   %  latitude, longitude : site point, required for the swd CSI variant.
   %  jump_factor : boundary-jump multiplier (POLICY B6).
    %  blend_hours : seam-blend taper window (POLICY B6).
   %  toa_dark_wm2 : irradiance threshold below which a sample counts as
    %     dark in the swd CSI mask.
   %  allow_swd_flux_fallback : false by default. The post-final D-32
   %     pass sets true so residual SWD gaps for which CSI is undefined
   %     use a capped flux-linear bridge. A run of at least two postings
   %     may instead use one finite daylight CSI anchor and one physically
   %     known darkness zero; the target itself never enters deep darkness.
   %     If that CSI candidate is blocked, the same scalar-valid flux bridge
   %     closes the residual without clipping it to the diagnostic solar
   %     relation. Only SWD within its D-49 nine-hour cap and
   %     D-48 RH within its evidence-backed nine-hour cap may cross a
   %     calendar-season boundary.
    %  step_scale : optional native-only scale frozen by the orchestrator
    %     before any reconstruction values are inserted.
   %  Defaults come from the central
   %  icemodel.forcing.reconstruct.setopts contract.
   %
   % Returns
   %  x : the channel with admitted short gaps filled.
   %  filled : logical column, true where this tier filled a sample.
   %  audit : segment-audit cell rows {channel, start, end, duration_hours,
   %     method, detail} — one per filled gap — for the orchestrator's
   %     audit table.
   %
   % See also: icemodel.forcing.reconstruct.reconstructSeries,
   %  icemodel.forcing.reconstruct.stepScale,
   %  icemodel.forcing.reconstruct.physicalBounds

   arguments
      times datetime
      x (:, 1) double
      channel (1, 1) string
      kwargs.cap_hours (1, 1) double ...
         {mustBePositive, ...
         icemodel.forcing.reconstruct.mustBeCapHours(kwargs.cap_hours)} = ...
         icemodel.forcing.reconstruct.setopts().cap_hours
      kwargs.latitude (1, 1) double = NaN
      kwargs.longitude (1, 1) double = NaN
      kwargs.jump_factor (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().jump_factor
      kwargs.blend_hours (1, 1) double {mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().blend_hours
       kwargs.toa_dark_wm2 (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().toa_dark_wm2
       kwargs.allow_swd_flux_fallback (1, 1) logical = false
       kwargs.step_scale (1, 1) struct = struct()
   end

   % SWU has no independent tier-1 path: it follows albedo*swd downstream.
   if channel == "swu"
      error('icemodel:reconstruct:fillShortGaps:channelNotInterpolable', ...
         '%s is derived from albedo and swd; use the downstream rule', ...
         channel);
   end
   caps = icemodel.forcing.reconstruct.interpolationCapHours();
   channel_ceiling = caps.default;
   if isfield(caps, channel)
      channel_ceiling = caps.(channel);
   end
   if kwargs.cap_hours > channel_ceiling
      error('icemodel:reconstruct:fillShortGaps:capHours', ...
         'cap_hours for %s must not exceed %g', channel, channel_ceiling);
   end

   dt_hours = hours(median(diff(times)));
   blend_len = max(1, round(kwargs.blend_hours / dt_hours));
   scale = kwargs.step_scale;
   if isempty(fieldnames(scale))
      scale = icemodel.forcing.reconstruct.stepScale(times, x);
   end

   % Shortwave interpolates in clear-sky-index space; the policy forbids
   % raw-linear daylight interpolation for radiation-cycle channels.
   use_csi = channel == "swd";
   if use_csi
      if ~(isfinite(kwargs.latitude) && isfinite(kwargs.longitude))
         error('icemodel:reconstruct:fillShortGaps:missingSolarGeometry', ...
            'swd short-gap filling requires latitude and longitude');
      end
      location = struct('lat_wgs84', kwargs.latitude, ...
         'lon_wgs84', kwargs.longitude);
      [work, toa] = icemodel.forcing.reconstruct.clearSkyIndex( ...
         times, x, location, toa_dark_wm2=kwargs.toa_dark_wm2);
      % Hoisted whole-posting maximum elevation keeps candidate validity
      % aligned with staging and the reconstruction darkness rule.
      axis_elevation = ...
         icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
         times, kwargs.latitude, kwargs.longitude, hours(dt_hours));
      solar_bands = ...
         icemodel.forcing.reconstruct.solarElevationBands();
   else
      work = x;
   end

   % Interior missing runs of the working signal, capped by wall-clock
   % length; leading/trailing gaps have no interior anchors and stay for
   % later tiers.
   miss = ~isfinite(work);
   if use_csi
      % Darkness is a hard context boundary, not part of a missing daylight
      % run. A dawn/dusk fragment therefore needs daylight anchors on both
      % sides and can never borrow an anchor across the intervening night.
      miss = ~isfinite(x) & toa >= kwargs.toa_dark_wm2;
   end
   d = diff([0; miss; 0]);
   starts = find(d == 1);
   stops = find(d == -1) - 1;
   filled = false(size(x));
   audit = cell(0, 1);
   for g = 1:numel(starts)
      before = starts(g) - 1;
      after = stops(g) + 1;
      if before < 1 || after > numel(work)
         continue
      end
      if ~(isfinite(work(before)) && isfinite(work(after)))
         continue
      end

      % Linear bridge in working space, reconverted for shortwave.
      work_idx = (starts(g):stops(g)).';
      candidate = interp1([before; after], ...
         [work(before); work(after)], work_idx);
      if use_csi
         candidate = candidate .* toa(work_idx);
      end

      % Apply the duration cap and audit independently to each contiguous
      % native gap. For CSI, MISS already splits runs at every dark sample.
      native_missing = ~isfinite(x(work_idx));
      native_d = diff([false; native_missing; false]);
      native_starts = find(native_d == 1);
      native_stops = find(native_d == -1) - 1;
      for r = 1:numel(native_starts)
         target = work_idx(native_starts(r):native_stops(r));
         duration = numel(target) * dt_hours;
         [season_ok, season] = seasonBridgeAllowed( ...
            times, before, target, after, channel, duration);
         if ~season_ok
            % Held-out interpolation evidence is season-contained, so both
            % observed anchors and the missing samples normally share a
            % season; D-48/D-49 own the narrow calendar-edge exceptions.
            continue
         end
         if use_csi && any(x([before; after]) > ...
               solar_bands.toa_ceiling_multiplier ...
               * solar_bands.solar_constant_wm2)
            % A relationally invalid native sample stays immutable, but it
            % cannot authorize a synthetic bridge (D-49).
            continue
         end
         if duration > kwargs.cap_hours
            continue
         end
         candidate_fill = ...
            candidate(native_starts(r):native_stops(r));

         % Physical bounds remain a hard refusal. Boundary mismatches use the
         % same bridge floor and taper as later admitted methods. The CSI
         % path passes slices of the hoisted axis geometry so the shared
         % validator skips its per-fragment NOAA recomputation.
         if use_csi
            physically_valid = ...
               icemodel.forcing.reconstruct.physicalValidity( ...
               channel, candidate_fill, times(target), ...
               latitude=kwargs.latitude, longitude=kwargs.longitude, ...
               toa=toa(target), elevation=axis_elevation(target));
         else
            physically_valid = ...
               icemodel.forcing.reconstruct.physicalValidity( ...
               channel, candidate_fill, times(target), ...
               latitude=kwargs.latitude, longitude=kwargs.longitude);
         end
         if ~all(physically_valid)
            continue
         end
         limit = kwargs.jump_factor * scale.(char(season));
         if ~isfinite(limit)
            limit = 0;
         end
         target_before = target(1) - 1;
         target_after = target(end) + 1;
         if target_before >= 1 && target_after <= numel(x) ...
               && isfinite(x(target_before)) && isfinite(x(target_after))
            bridge_step = abs(x(target_after) - x(target_before)) ...
               / (numel(target) + 1);
            limit = max(limit, bridge_step * (1 + 1e-6));
         end
         [candidate_fill, seam_note] = ...
            icemodel.forcing.reconstruct.blendSeams( ...
            x, candidate_fill, target, true(size(target)), ...
            target(1), target(end), target_before, target_after, ...
            blend_len, limit);
         % Post-blend recheck under the same hoisted-geometry contract as
         % the pre-blend refusal above.
         if use_csi
            physically_valid = ...
               icemodel.forcing.reconstruct.physicalValidity( ...
               channel, candidate_fill, times(target), ...
               latitude=kwargs.latitude, longitude=kwargs.longitude, ...
               toa=toa(target), elevation=axis_elevation(target));
         else
            physically_valid = ...
               icemodel.forcing.reconstruct.physicalValidity( ...
               channel, candidate_fill, times(target), ...
               latitude=kwargs.latitude, longitude=kwargs.longitude);
         end
         if ~all(physically_valid)
            continue
         end
         jump_ok = true;
         if target_before >= 1 && isfinite(x(target_before))
            jump_ok = ...
               abs(candidate_fill(1) - x(target_before)) <= limit;
         end
         if target_after <= numel(x) && isfinite(x(target_after))
            jump_ok = jump_ok ...
               && abs(candidate_fill(end) - x(target_after)) <= limit;
         end
         if ~jump_ok
            continue
         end

         % Fill only the natively missing samples; never touch finite data.
       x(target) = candidate_fill;
       filled(target) = true;
       method = "bounded_interp";
         detail = sprintf('linear, cap %.3g h', kwargs.cap_hours);
         if use_csi
            detail = sprintf('csi-linear, cap %.3g h', kwargs.cap_hours);
         end
         detail = [detail seam_note]; %#ok<AGROW>
       segment = false(numel(times), 1);
       segment(target) = true;
       rows = icemodel.forcing.reconstruct.auditSegments( ...
          times, segment, channel, method, detail);
       audit = [audit; rows]; %#ok<AGROW>
      end
    end

   % CSI intentionally has no denominator near darkness. On the D-32
   % post-final pass only, retry the still-missing short SWD slivers in
   % flux space. Real KANL holdouts show this local bridge halves
   % shoulder-hour RMSE versus persistence while CSI has no coverage;
   % ordinary tier 1 remains CSI-only.
   if use_csi && kwargs.allow_swd_flux_fallback && any(~isfinite(x))
      [x, flux_filled, flux_audit] = fillShortwaveFluxResiduals( ...
         times, x, kwargs.cap_hours, kwargs.latitude, kwargs.longitude, ...
         kwargs.toa_dark_wm2);
      filled = filled | flux_filled;
      audit = [audit; flux_audit];
   end
end

function [x, filled, audit] = fillShortwaveFluxResiduals( ...
      times, x, cap_hours, latitude, longitude, toa_dark_wm2)
   %FILLSHORTWAVEFLUXRESIDUALS Bridge capped post-final SWD slivers.
   % Two finite solar anchors use flux-linear interpolation. A run of at
   % least two still-sunlit postings beside one known darkness zero uses
   % the opposite finite anchor's clear-sky index over the target TOA
   % curve. The target never enters deep darkness. The flux-linear bridge
   % retains its evidence-backed diagnostic exception to the pointwise
   % solar relation; the one-sided CSI path must pass the shared SWD
   % physical-validity ceiling.
   dt_hours = hours(median(diff(times)));
   maximum_elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, latitude, longitude, hours(dt_hours));
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   location = struct('lat_wgs84', latitude, 'lon_wgs84', longitude);
   [csi, toa] = icemodel.forcing.reconstruct.clearSkyIndex( ...
      times, x, location, toa_dark_wm2=toa_dark_wm2);
   eligible = ~isfinite(x) ...
      & maximum_elevation > bands.civil_twilight_deg;
   edges = diff([false; eligible; false]);
   starts = find(edges == 1);
   stops = find(edges == -1) - 1;
   filled = false(size(x));
   audit = cell(0, 1);
   for g = 1:numel(starts)
      target = (starts(g):stops(g)).';
      before = target(1) - 1;
      after = target(end) + 1;
      if numel(target) * dt_hours > cap_hours ...
            || before < 1 || after > numel(x) ...
            || ~isfinite(x(before)) || ~isfinite(x(after))
         continue
      end
      before_sunlit = ...
         maximum_elevation(before) > bands.civil_twilight_deg;
      after_sunlit = ...
         maximum_elevation(after) > bands.civil_twilight_deg;
      before_darkness = ...
         maximum_elevation(before) <= bands.civil_twilight_deg ...
         && x(before) == 0;
      after_darkness = ...
         maximum_elevation(after) <= bands.civil_twilight_deg ...
         && x(after) == 0;
      duration = numel(target) * dt_hours;
      [season_ok, ~] = seasonBridgeAllowed( ...
         times, before, target, after, "swd", duration);
      if ~season_ok
         continue
      end
      anchors = x([before; after]);
      anchors_plausible = all(anchors <= ...
         bands.toa_ceiling_multiplier * bands.solar_constant_wm2);

      if before_sunlit && after_sunlit && anchors_plausible
         candidate = interp1([before; after], ...
            [x(before); x(after)], target);
         detail = sprintf( ...
            'flux-linear fallback, cap %.3g h; ', cap_hours);
         candidate_valid = ...
            icemodel.forcing.reconstruct.scalarValidity("swd", candidate);
      elseif numel(target) >= 2 && before_sunlit && after_darkness ...
            && isfinite(csi(before))
         candidate = csi(before) .* toa(target);
         detail = sprintf( ...
            'one-sided CSI to known darkness, cap %.3g h; ', cap_hours);
         candidate_valid = ...
            icemodel.forcing.reconstruct.physicalValidity( ...
            "swd", candidate, times(target), latitude=latitude, ...
            longitude=longitude, toa=toa(target), ...
            elevation=maximum_elevation(target), ...
            interval=hours(dt_hours));
      elseif numel(target) >= 2 && before_darkness && after_sunlit ...
            && isfinite(csi(after))
         candidate = csi(after) .* toa(target);
         detail = sprintf( ...
            'one-sided CSI from known darkness, cap %.3g h; ', cap_hours);
         candidate_valid = ...
            icemodel.forcing.reconstruct.physicalValidity( ...
            "swd", candidate, times(target), latitude=latitude, ...
            longitude=longitude, toa=toa(target), ...
            elevation=maximum_elevation(target), ...
            interval=hours(dt_hours));
      else
         candidate_valid = false;
      end
      % D-32 makes the solar relation diagnostic for this final residual
      % bridge. If CSI is unavailable or fails its relational ceiling, use
      % the finite local anchors directly and retain only the hard scalar
      % bound; clipping would manufacture the very seam this pass closes.
      can_flux_bridge = numel(target) >= 2 ...
         && any(isfinite(csi([before; after]))) ...
         && anchors_plausible;
      if ~all(candidate_valid) && can_flux_bridge
         candidate = interp1([before; after], ...
            [x(before); x(after)], target);
         detail = sprintf( ...
            'flux-linear fallback, cap %.3g h; ', cap_hours);
         candidate_valid = ...
            icemodel.forcing.reconstruct.scalarValidity("swd", candidate);
      end
      if ~all(candidate_valid)
         continue
      end
      x(target) = candidate;
      filled(target) = true;
      segment = false(numel(times), 1);
      segment(target) = true;
      rows = icemodel.forcing.reconstruct.auditSegments( ...
         times, segment, "swd", "bounded_interp", ...
         detail + "solar relation diagnostic");
      audit = [audit; rows]; %#ok<AGROW>
   end
end

function [ok, season] = seasonBridgeAllowed( ...
      times, before, target, after, channel, duration_hours)
   %SEASONBRIDGEALLOWED Apply the shared interpolation season boundary.
   target_seasons = icemodel.forcing.reconstruct.seasonOf(times(target));
   bridge_seasons = icemodel.forcing.reconstruct.seasonOf( ...
      times([before; target; after]));
   crosses_boundary = any(bridge_seasons ~= bridge_seasons(1));
   caps = icemodel.forcing.reconstruct.interpolationCapHours();
   % A calendar-season label is not a physical discontinuity. D-49 applies
   % SWD's existing ceiling to the whole local bridge, including samples
   % that straddle the boundary; D-48 does the same within RH's
   % cap. Other state channels retain the season guard.
   boundary_exception = (channel == "swd" ...
      && duration_hours <= caps.swd_season_boundary) ...
      || (channel == "rh" && duration_hours <= caps.rh);
   ok = ~crosses_boundary || boundary_exception;
   season = target_seasons(1);
end
