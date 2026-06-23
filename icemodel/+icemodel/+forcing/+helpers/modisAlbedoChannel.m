function modis = modisAlbedoChannel(modis_dir, years, location, method, ...
      remap, Time)
   %MODISALBEDOCHANNEL GEUS MODIS daily albedo on a time axis, per location.
   %
   %  modis = icemodel.forcing.helpers.modisAlbedoChannel(modis_dir, years, ...
   %     location, method, remap, Time)
   %
   % Reads the GEUS MODIS daily albedo for each requested year and interpolates
   % it onto the hourly TIME axis, returning one column aligned to TIME. This is
   % the single shared implementation behind the optional MODIS channel of every
   % gridded-source builder (buildMarData / buildMerraData / buildRacmoData), so
   % all three resolve MODIS identically.
   %
   % LOCATION is the ORIGINAL build request (a [lat lon] point or an EPSG:3413
   % polyshape); readGeusModis maps it onto the GEUS 5 km grid with the same
   % point (nearest/natural) or polygon (conservative/equal) selection as the
   % builder's other gridded channels, so a catchment build gets the
   % area-weighted ROI mean rather than the single nearest cell. One GEUS
   % reflectivity file per year is expected under MODIS_DIR (one match required).
   %
   % Inputs
   %  modis_dir - directory with GEUS Greenland_Reflectivity_<YYYY>_5km_C6.nc
   %  years     - calendar years to read (one MODIS file per year)
   %  location  - [lat lon] point or EPSG:3413 polyshape (the build request)
   %  method    - point sampling "nearest" | "natural"
   %  remap     - polygon aggregation "conservative" | "equal"
   %  Time      - target datetime axis the daily albedo is interpolated onto
   %
   % Outputs
   %  modis - daily MODIS albedo interpolated to TIME [-] (NaN where no year
   %          covers a sample)
   %
   % See also: icemodel.forcing.readGeusModis,
   %  icemodel.forcing.helpers.dailyToHourly, icemodel.forcing.buildMarData

   % NaN-initialized so any sample whose year is not requested stays missing.
   modis = nan(numel(Time), 1);
   for yyyy = years
      % One GEUS reflectivity file per year; a missing/ambiguous match is an
      % error rather than a silent gap.
      match = dir(fullfile(modis_dir, sprintf('*_%d_*.nc', yyyy)));
      if numel(match) ~= 1
         error('icemodel:forcing:modisAlbedoChannel:fileNotFound', ...
            'expected one MODIS file for %d in %s, found %d', ...
            yyyy, modis_dir, numel(match))
      end
      % Read the daily series at the build location, then interpolate this
      % year's days onto the matching slice of the target axis.
      [albedo, Tdaily] = icemodel.forcing.readGeusModis( ...
         string(fullfile(match.folder, match.name)), location, method, ...
         remap=remap);
      inyear = year(Time) == yyyy;
      modis(inyear) = icemodel.forcing.helpers.dailyToHourly( ...
         albedo, Tdaily, Time(inyear));
   end
end
