function x = computeCumulativeWindowDelta(ice1, varname, t1, t2, timevec)
%COMPUTECUMULATIVEWINDOWDELTA Return cumulative-variable change over a window.
%
%  x = test.helpers.computeCumulativeWindowDelta(ice1, "melt", t1, t2)
%  x = test.helpers.computeCumulativeWindowDelta(ice1, "runoff", t1, t2, timevec)
%
% Inputs:
%  ice1    - timetable or struct containing a cumulative variable
%  varname - variable name, e.g. "melt" or "runoff"
%  t1, t2  - inclusive evaluation window
%  timevec - optional datetime vector (defaults to ice1.Time if available)
%
% Output:
%  x - cumulative change over [t1, t2]

   arguments
      ice1
      varname (1, 1) string
      t1 (1, 1) datetime
      t2 (1, 1) datetime
      timevec = []
   end

   if nargin < 5 || isempty(timevec)
      if istimetable(ice1)
         timevec = ice1.Time;
      elseif isstruct(ice1) && isfield(ice1, 'Time')
         timevec = ice1.Time;
      else
         error('time vector is required when ice1 does not contain Time')
      end
   end

   x = nan;
   series = getSeries(ice1, varname);
   series = series(:);
   timevec = timevec(:);

   if numel(series) ~= numel(timevec)
      error('%s and time vector length mismatch', varname)
   end
   if ~isdatetime(timevec)
      error('timevec must be datetime')
   end

   inc = [0; diff(series, 1, 1)];
   keep = isbetween(timevec, t1, t2);
   if ~any(keep)
      return
   end

   inc_window = inc(keep);
   inc_window(~isfinite(inc_window)) = 0;
   if ~isempty(inc_window)
      x = sum(inc_window, 1);
   end
end

function x = getSeries(ice1, varname)
   name = char(varname);
   if istimetable(ice1)
      if ~ismember(name, ice1.Properties.VariableNames)
         error('ice1 timetable missing %s variable', varname)
      end
      x = ice1.(name);
   elseif isstruct(ice1) && isfield(ice1, name)
      x = ice1.(name);
   else
      error('ice1 must provide %s', varname)
   end

   if ~isnumeric(x)
      error('%s must be numeric', varname)
   end
   if size(x, 2) > 1
      x = x(:, 1);
   end
end
