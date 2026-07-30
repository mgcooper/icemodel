function [candidate, note] = blendFallbackSeams(times, native, current, ...
      candidate, fill_mask, kwargs)
   %BLENDFALLBACKSEAMS Apply the policy seam taper to fallback fills.
   %
   %  [candidate, note] = ...
   %     icemodel.forcing.reconstruct.blendFallbackSeams( ...
   %     times, native, current, candidate, fill_mask)
   %
   % `native` fixes the station/season step scale before reconstruction,
   % while `current` supplies the anchors bordering each fallback segment.
   % The returned note is suitable for appending to segment-audit detail.

   arguments
      times (:, 1) datetime
      native (:, 1) double
      current (:, 1) double
      candidate (:, 1) double
      fill_mask (:, 1) logical
      kwargs.jump_factor (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().jump_factor
      kwargs.blend_hours (1, 1) double {mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().blend_hours
   end

   n = numel(times);
   if numel(native) ~= n || numel(current) ~= n ...
         || numel(candidate) ~= n || numel(fill_mask) ~= n
      error('icemodel:reconstruct:blendFallbackSeams:sizeMismatch', ...
         'times, native, current, candidate, and fill_mask must align');
   end
   note = '';
   if n < 2 || ~any(fill_mask)
      return
   end

   % Reuse the admitted-method taper for each contiguous fallback segment.
   dt_hours = hours(median(diff(times)));
   if ~isfinite(dt_hours) || dt_hours <= 0
      error('icemodel:reconstruct:blendFallbackSeams:invalidCadence', ...
         'times must have a finite positive cadence');
   end
   scale = icemodel.forcing.reconstruct.stepScale(times, native);
   blend_len = max(1, round(kwargs.blend_hours / dt_hours));
   edge = diff([false; fill_mask; false]);
   starts = find(edge == 1);
   stops = find(edge == -1) - 1;
   notes = strings(numel(starts), 1);
   n_notes = 0;
   for g = 1:numel(starts)
      idx = (starts(g):stops(g)).';
      before = starts(g) - 1;
      after = stops(g) + 1;
      season = icemodel.forcing.reconstruct.seasonOf(times(starts(g)));
      limit = kwargs.jump_factor * scale.(char(season));
      if before >= 1 && after <= n ...
            && isfinite(current(before)) && isfinite(current(after))
         bridge_step = abs(current(after) - current(before)) ...
            / (numel(idx) + 1);
         limit = max(limit, bridge_step * (1 + 1e-6));
      end
      [candidate(idx), local_note] = ...
         icemodel.forcing.reconstruct.blendSeams( ...
         current, candidate(idx), idx, true(size(idx)), ...
         starts(g), stops(g), before, after, blend_len, limit);
      if ~isempty(local_note)
         n_notes = n_notes + 1;
         notes(n_notes) = string(local_note);
      end
   end
   % Era-scale adoptions (a multi-year outage union at stations whose
   % native record never carried a channel) produce dozens of seams; a
   % verbatim concatenation makes every audit row unreadable, so beyond a
   % handful the note collapses to a count. The per-sample product and
   % provenance remain the detailed record.
   if n_notes > 4
      note = char(", seam blend: " + n_notes + " seams tapered");
   else
      note = char(strjoin(notes(1:n_notes), ""));
   end
end
