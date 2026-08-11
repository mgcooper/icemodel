function keep = timeWindowMask(Time, startdate, enddate)
   %TIMEWINDOWMASK Select an optional datetime window from a source time axis.

   % Enforce the paired-window contract, then keep every sample when no window
   % is supplied, independent of the source bounds.
   [t0, t1, has_window] = icemodel.pairedWindow(startdate, enddate);
   keep = true(size(Time));
   if has_window
      keep = Time >= t0 & Time <= t1;
   end
end
