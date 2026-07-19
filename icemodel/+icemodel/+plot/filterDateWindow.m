function T = filterDateWindow(T, startdate, enddate)
   %FILTERDATEWINDOW Apply a shared inclusive date window to plot inputs.
   %
   %  T = icemodel.plot.filterDateWindow(T, startdate, enddate)
   %
   % Supports timetables, tables with a time/date column, and interval tables with
   % start_date/end_date columns. Empty tables return unchanged.

   if isequal(startdate, "") && isequal(enddate, "")
      return
   end
   if isempty(T)
      return
   end

   if istimetable(T)
      T = filterByVector(T, T.Time, startdate, enddate, "containment");
      return
   end

   names = string(T.Properties.VariableNames);
   time_name = names(find(ismember(lower(names), ["time", "datetime", ...
      "date"]), 1));
   if ~isempty(time_name)
      T = filterByVector(T, T.(char(time_name)), startdate, enddate, ...
         "containment");
      return
   end

   if all(ismember(["start_date", "end_date"], names))
      T = filterByInterval(T, startdate, enddate);
   end
end

function T = filterByVector(T, time, startdate, enddate, mode)
   %FILTERBYVECTOR Keep rows whose timestamp falls inside the requested window.
   if ~isdatetime(time)
      time = datetime(time);
   end
   t0 = icemodel.plot.parseDate(startdate, time(1));
   t1 = icemodel.plot.parseDate(enddate, time(end));
   if mode == "containment"
      T = T(time >= t0 & time <= t1, :);
   end
end

function T = filterByInterval(T, startdate, enddate)
   %FILTERBYINTERVAL Keep rows whose intervals overlap the requested window.
   start_time = datetime(T.start_date);
   end_time = datetime(T.end_date);
   t0 = icemodel.plot.parseDate(startdate, start_time(1));
   t1 = icemodel.plot.parseDate(enddate, end_time(end));
   T = T(end_time >= t0 & start_time <= t1, :);
end
