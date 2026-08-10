function cadence_seconds = uniformCadenceSeconds(value)
   %UNIFORMCADENCESECONDS Derive one timetable's exact regular cadence.
   %
   %  cadence_seconds = icemodel.forcing.helpers.uniformCadenceSeconds(value)
   %
   % VALUE is a timetable. CADENCE_SECONDS is its sample spacing in seconds, or
   % NaN if the timetable has fewer than two rows or the spacing is not
   % constant.

   arguments
      value timetable
   end

   cadence_seconds = NaN;
   if height(value) < 2
      return
   end

   % Compare every step against the median rather than the first step so a
   % single corrupt leading sample cannot define the accepted cadence. The
   % 1e-6 s tolerance absorbs datetime round-off without admitting a real gap.
   steps = seconds(diff(value.Time));
   candidate = median(steps, 'omitnan');
   if isfinite(candidate) && candidate > 0 ...
         && all(isfinite(steps)) && all(abs(steps - candidate) < 1e-6)
      cadence_seconds = candidate;
   end
end
