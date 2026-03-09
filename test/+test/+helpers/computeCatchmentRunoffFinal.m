function final_m3 = computeCatchmentRunoffFinal(ice1, area_med_m2, t1, t2, timevec)
%COMPUTECATCHMENTRUNOFFFINAL Convert point runoff to catchment cumulative final.
%
%  final_m3 = test.helpers.computeCatchmentRunoffFinal(ice1, area_med_m2, t1, t2)
%  final_m3 = test.helpers.computeCatchmentRunoffFinal(ice1, area_med_m2, t1, t2, timevec)
%
% Inputs:
%  ice1        - timetable or struct containing cumulative point runoff [m]
%  area_med_m2 - catchment med area [m^2]
%  t1, t2      - evaluation window
%  timevec     - optional datetime vector (defaults to ice1.Time if available)
%
% Output:
%  final_m3 - final cumulative catchment runoff [m^3]

   if nargin < 5 || isempty(timevec)
      if istimetable(ice1)
         timevec = ice1.Time;
      elseif isstruct(ice1) && isfield(ice1, 'Time')
         timevec = ice1.Time;
      else
         error('time vector is required when ice1 does not contain Time')
      end
   end

   runoff = getRunoffSeries(ice1);
   runoff = runoff(:);
   timevec = timevec(:);

   if numel(runoff) ~= numel(timevec)
      error('runoff and time vector length mismatch')
   end

   if ~isdatetime(timevec)
      error('timevec must be datetime')
   end

   r_inc = [0; diff(runoff, 1, 1)];
   keep = isbetween(timevec, t1, t2);
   if ~any(keep)
      final_m3 = nan;
      return
   end

   r_inc_window = r_inc(keep);
   r_inc_window(~isfinite(r_inc_window)) = 0;

   catchment_cum_m3 = area_med_m2 * cumsum(r_inc_window, 1);
   if isempty(catchment_cum_m3)
      final_m3 = nan;
   else
      final_m3 = catchment_cum_m3(end);
   end
end

function runoff = getRunoffSeries(ice1)
   if istimetable(ice1)
      if ~ismember('runoff', ice1.Properties.VariableNames)
         error('ice1 timetable missing runoff variable')
      end
      runoff = ice1.runoff;
   elseif isstruct(ice1) && isfield(ice1, 'runoff')
      runoff = ice1.runoff;
   else
      error('ice1 must provide runoff')
   end

   if ~isnumeric(runoff)
      error('runoff must be numeric')
   end
   if size(runoff, 2) > 1
      runoff = runoff(:, 1);
   end
end
