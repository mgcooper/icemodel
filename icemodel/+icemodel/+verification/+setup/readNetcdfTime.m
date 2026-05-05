function time = readNetcdfTime(pathname, varname)
   %READNETCDFTIME Read a NetCDF time coordinate as UTC datetime.
   %
   %  time = icemodel.verification.setup.readNetcdfTime(pathname, varname)
   %
   %  Reads numeric time values from a NetCDF variable and converts them
   %  to UTC datetime using the variable's "units" attribute. Currently
   %  supports "hours since YYYY-MM-DD HH:MM:SS.S" reference strings,
   %  which is the convention used by the ESM-SnowMIP PANGAEA bundle.
   %
   %  Inputs
   %    pathname : string
   %        Absolute path to a NetCDF file.
   %    varname : string
   %        Name of the time-coordinate variable (typically "time").
   %
   %  Returns
   %    time : datetime column (TimeZone = 'UTC')
   %
   % See also: icemodel.verification.setup.readObsChannel,
   %  icemodel.verification.setup.buildEsmSnowmipForcing,
   %  icemodel.verification.setup.buildEsmSnowmipObservations

   raw = double(ncread(pathname, varname));
   units = string(ncreadatt(pathname, varname, 'units'));
   tref = datetime(extractAfter(units, 'hours since '), ...
      'InputFormat', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
   time = tref + hours(raw);
   time = time(:);
end
