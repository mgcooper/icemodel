function [albedo, Time] = readGeusModis(filename, lat, lon)
   %READGEUSMODIS Read the GEUS MODIS daily albedo at a point.
   %
   %  [albedo, Time] = icemodel.forcing.readGeusModis(filename, lat, lon)
   %
   % Reads the daily MODIS-derived surface albedo (fraction) from one
   % GEUS Greenland reflectivity file (Greenland_Reflectivity_<YYYY>_
   % 5km_C6.nc) at the grid cell nearest the requested point.
   %
   % Inputs
   %  filename - GEUS reflectivity NetCDF for one year
   %  lat, lon - point location [degrees]
   %
   % Outputs
   %  albedo - daily albedo series at the nearest cell [-]
   %  Time   - UTC daily datetime axis (Jan 1 of the file year onward)
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.forcing.helpers.dailyToHourly

   arguments
      filename (1, 1) string
      lat (1, 1) double
      lon (1, 1) double
   end

   LON = double(ncread(filename, 'lon'));
   LAT = double(ncread(filename, 'lat'));

   % Nearest cell by projected distance (the 5 km grid is regular in
   % polar stereographic space).
   proj = icemodel.forcing.helpers.psnProjection();
   [X, Y] = projfwd(proj, LAT, LON);
   [xq, yq] = projfwd(proj, lat, lon);
   [~, idx] = min(hypot(X(:) - xq, Y(:) - yq));
   [i, j] = ind2sub(size(X), idx);

   info = ncinfo(filename, 'albedo');
   ndays = info.Size(end);
   albedo = squeeze(double(ncread(filename, 'albedo', [i j 1], ...
      [1 1 ndays])));

   tok = regexp(filename, '_(\d{4})_', 'tokens', 'once');
   assert(~isempty(tok), ...
      'cannot parse the file year from %s', filename)
   t0 = datetime(str2double(tok{1}), 1, 1, 'TimeZone', 'UTC');
   Time = (t0:days(1):(t0 + days(ndays - 1)))';
end
