function [candidate, note] = blendSeams(x, candidate, idx, usable, ...
      run_start, run_stop, run_before, run_after, blend_len, limit)
   %BLENDSEAMS Taper excess anchored boundary mismatch across one run.
   %
   % Each anchored edge whose sample this method fills is tested against
   % the jump limit. Larger mismatches are pulled to half the limit by a
   % linearly decaying offset; seams already inside the limit are unchanged.

   n = numel(idx);
   note = "";
   w = min(blend_len, n);
   first_use = find(usable, 1, 'first');
   last_use = find(usable, 1, 'last');

   % A candidate can be usable in another residual segment but have no
   % usable samples in this one; leave that segment for the next method.
   if isempty(first_use)
      note = char(note);
      return
   end

   % Zero offsets mean an unanchored, unfilled, or already-valid seam.
   off_start = 0;
   off_end = 0;
   if idx(first_use) == run_start && run_before >= 1 ...
         && isfinite(x(run_before))
      delta = x(run_before) - candidate(first_use);
      if abs(delta) > limit
         off_start = delta - sign(delta) * limit / 2;
         note = note + sprintf(" start %+.3g", off_start);
      end
   end
   if idx(last_use) == run_stop && run_after <= numel(x) ...
         && isfinite(x(run_after))
      delta = x(run_after) - candidate(last_use);
      if abs(delta) > limit
         off_end = delta - sign(delta) * limit / 2;
         note = note + sprintf(" end %+.3g", off_end);
      end
   end

   % Short runs use one correction connecting both offsets; longer runs
   % keep independent ramps and preserve their interior signal.
   if n < 2 * w
      if n == 1
         both_anchored = idx(first_use) == run_start ...
            && run_before >= 1 && isfinite(x(run_before)) ...
            && idx(last_use) == run_stop ...
            && run_after <= numel(x) && isfinite(x(run_after));
         if both_anchored && (off_start ~= 0 || off_end ~= 0)
            correction = (x(run_before) + x(run_after)) / 2 ...
               - candidate(first_use);
            note = note + " (bridge midpoint)";
         else
            correction = off_start + off_end;
         end
      else
         correction = off_start + (off_end - off_start) ...
            * ((0:n - 1).' / (n - 1));
      end
   else
      correction = zeros(n, 1);
      % Include zero at each interior taper endpoint; a one-sample window
      % can only correct its anchored boundary sample.
      if w == 1
         taper = 1;
      else
         taper = linspace(1, 0, w).';
      end
      correction(1:w) = off_start * taper;
      correction(n - w + 1:n) = off_end * flip(taper);
   end
   candidate = candidate + correction;
   if strlength(note) > 0
      note = ", seam blend:" + note;
   end
   note = char(note);
end
