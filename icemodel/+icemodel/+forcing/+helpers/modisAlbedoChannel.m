function [modis, metadata] = modisAlbedoChannel(modis_dir, years, location, method, ...
      remap, Time)
   %MODISALBEDOCHANNEL GEUS MODIS daily albedo on a time axis, per location.
   %
   %  [modis, metadata] = ...
   %     icemodel.forcing.helpers.modisAlbedoChannel(modis_dir, years, ...
   %     location, method, remap, Time)
   %
   % Reads the GEUS MODIS daily albedo for each requested year and interpolates
   % it onto the hourly TIME axis. It returns one column aligned to TIME. Every
   % gridded-source builder (buildMarData / buildMerraData / buildRacmoData)
   % calls this function for its optional MODIS channel, so all three resolve
   % MODIS the same way.
   %
   % LOCATION is the ORIGINAL build request: a [lat lon] point or an EPSG:3413
   % polyshape. readGeusModis maps it onto the GEUS 5 km grid with the same
   % point (nearest/natural) or polygon (conservative/equal) selection as the
   % builder's other gridded channels. A catchment build therefore gets the
   % area-weighted ROI mean, not the single nearest cell. A missing year stays
   % NaN, because MODIS is an optional diagnostic channel. A duplicate match is
   % an error, because it makes the source layout ambiguous.
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
   %  modis   - daily MODIS albedo interpolated to TIME [-] (NaN where no year
   %            covers a sample)
   %  metadata - product/status/exact-coverage provenance. A caller
   %             writing an artifact omits MODIS when coverage is empty.
   %
   % See also: icemodel.forcing.readGeusModis,
   %  icemodel.forcing.helpers.dailyToHourly, icemodel.forcing.buildMarData

   % Initialize with NaN, so any sample without source coverage stays missing.
   % The target axis defines the requested artifact years in metadata, not the
   % caller's vector orientation or duplicate years.
   modis = nan(numel(Time), 1);
   requested_years = unique(year(Time))';
   source_years = unique(reshape(years, 1, []), 'stable');
   % Each source year contributes at most one covered year, so the coverage
   % buffer is sized to the source list once and trimmed after the loop.
   coverage_years = zeros(1, numel(source_years));
   n_coverage = 0;
   for yyyy = source_years
      inyear = year(Time) == yyyy;
      if ~any(inyear)
         continue
      end

      % Missing optional MODIS years stay NaN. Ambiguous matches are still a
      % source-layout error because the requested year cannot be selected.
      match = dir(fullfile(modis_dir, sprintf('*_%d_*.nc', yyyy)));
      if isempty(match)
         continue
      end
      if numel(match) > 1
         error('icemodel:forcing:modisAlbedoChannel:ambiguousFile', ...
            'expected one MODIS file for %d in %s, found %d', ...
            yyyy, modis_dir, numel(match))
      end
      % Read the daily series at the build location, then interpolate this
      % year's days onto the matching slice of the target axis.
      [albedo, Tdaily] = icemodel.forcing.readGeusModis( ...
         string(fullfile(match.folder, match.name)), location, method, ...
         remap=remap);
      values = icemodel.forcing.helpers.dailyToHourly( ...
         albedo, Tdaily, Time(inyear));
      valid = isfinite(values) & imag(values) == 0 ...
         & real(values) >= 0 & real(values) <= 1;
      if ~any(valid)
         error('icemodel:forcing:modisAlbedoChannel:noPhysicalValues', ...
            ['MODIS source file for %d produced no finite physical values ' ...
            'on the target axis'], yyyy)
      end
      modis(inyear) = values;
      n_coverage = n_coverage + 1;
      coverage_years(n_coverage) = yyyy;
   end
   % Drop the unused tail so no-coverage builds keep the 1-by-0 shape the
   % provenance helper expects.
   coverage_years = coverage_years(1:n_coverage);

   % Reuse the exact source matches above: provenance adds no directory scan or
   % NetCDF read beyond the physical channel construction.
   metadata = icemodel.forcing.helpers.geusModisCoverageMetadata( ...
      requested_years, coverage_years);
end
