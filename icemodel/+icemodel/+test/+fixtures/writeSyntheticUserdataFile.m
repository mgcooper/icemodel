function filepath = writeSyntheticUserdataFile(userdatadir, simyear, kwargs)
   %WRITESYNTHETICUSERDATAFILE Write a synthetic yearly userdata timetable.
   %
   %  filepath = icemodel.test.fixtures.writeSyntheticUserdataFile(userdatadir, simyear)

   arguments
      userdatadir (1, :) char
      simyear (1, 1) double {mustBeInteger, mustBePositive}
      kwargs.sitename (1, :) char = 'kanm'
      kwargs.userdata (1, :) char = 'modis'
      kwargs.varname (1, :) char = 'modis'
      kwargs.nsteps (1, 1) double {mustBeInteger, mustBePositive} = 24
      kwargs.dt_seconds (1, 1) double {mustBePositive} = 3600
      kwargs.Time datetime = datetime.empty()
      kwargs.values double = []
   end

   % Build the synthetic time vector unless the caller supplied one.
   if isempty(kwargs.Time)
      dt = seconds(kwargs.dt_seconds);
      t0 = datetime(simyear, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
      Time = transpose(t0 + dt * (0:kwargs.nsteps-1));
   else
      Time = kwargs.Time;
   end

   % Fill one smooth userdata series unless explicit values were provided.
   if isempty(kwargs.values)
      tod_hours = hour(Time) + minute(Time) / 60 + second(Time) / 3600;
      values = min(max(0.55 + 0.03 * sin(2 * pi * tod_hours / 24), 0.05), 0.95);
   else
      values = kwargs.values;
   end
   values = reshape(values, [], 1);

   % Save the file in the same yearly naming scheme as real userdata inputs.
   Data = timetable(Time, values, 'VariableNames', {kwargs.varname});

   filepath = fullfile(userdatadir, ...
      [kwargs.sitename '_' kwargs.userdata '_' int2str(simyear) '.mat']);
   save(filepath, 'Data');
end
