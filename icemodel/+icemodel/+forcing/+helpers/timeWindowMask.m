function keep = timeWindowMask(Time, startdate, enddate)
   %TIMEWINDOWMASK Select an optional datetime window from a source time axis.

   % Enforce the same paired-window contract at every public builder that uses
   % this helper, then leave the all-available path independent of source bounds.
   [t0, t1, has_window] = icemodel.internal.pairedWindow(startdate, enddate);
   keep = true(size(Time));
   if has_window
      keep = Time >= t0 & Time <= t1;
   end
end
