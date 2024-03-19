function data = ncread(source, varnames, start, count, stride)
   %NCREAD Read icemodel nc file into memory
   %
   %  DATA = NCREAD(SOURCE, VARNAMES, START, COUNT, STRIDE)
   %
   % See also:

   % This is meant to be a lower-level function, see loadSectorRunoff for a
   % higher-level version that transposes the data to mimic the old data files

   arguments
      source (1, :) char {mustBeFile}
      varnames (1, :) string {mustBeText} = "runoff"
      start (1, :) double {mustBeNumeric} = [1 1]
      count (1, :) double {mustBeNumeric} = [2479 8760]
      stride (1, :) double {mustBeNumeric} = ones(size(start))
   end

   % info = ncinfo(source);

   % Read in and transpose for consistency with og data files
   for varname = varnames(:)'
      % ncread uses matlab 1-based indexing, I think getVar starts at 0.
      data.(varname) = transpose(ncread(source, varname, start, count, stride));
   end

   % Read the time
   Time = ncread(source, 'time', start(2), count(2), stride(2));
   unit = ncreadatt(source, 'time', 'units');
   Time = datetime(strrep(unit, 'seconds since ', '')) + seconds(Time);

   % Read the x, y
   data.x = ncread(source, 'x_easting');
   data.y = ncread(source, 'y_northing');
   data.Time = Time;
end
