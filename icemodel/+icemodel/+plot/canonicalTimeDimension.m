function T = canonicalTimeDimension(T)
   %CANONICALTIMEDIMENSION Enforce the icemodel timetable row-time name.
   T.Properties.DimensionNames{1} = 'Time';
end
