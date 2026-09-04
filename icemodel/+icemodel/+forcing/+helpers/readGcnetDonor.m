function donor = readGcnetDonor(filename)
   %READGCNETDONOR Load one GC-Net surface file as an observed-only donor.
   %
   %  donor = icemodel.forcing.helpers.readGcnetDonor(filename)
   %
   % Each mapped channel keeps its own origin mask. A sample without a
   % per-sample origin flag is not eligible, because nothing in the file shows
   % that the value is native.
   %
   % See also: icemodel.forcing.reconstruct.fillPromiceStation,
   %  icemodel.forcing.helpers.gcnetHourlyAxis

   arguments
      filename (1, 1) string
   end

   map = [ ...
      "Ta_2m", "tair"
      "RH_2m", "rh"
      "WS_10m", "wspd"
      "SRin", "swd"
      "LRin", "lwd"];
   try
      t = ncread(filename, 'time');
      % The raw fractional-day coordinate drifts and lands off the hour
      % for most rows. The shared row-index convention gives the exact
      % hourly axis that both the builder and the donor use.
      times = icemodel.forcing.helpers.gcnetHourlyAxis( ...
         icemodel.forcing.helpers.gcnetTime(t, ...
         ncreadatt(filename, 'time', 'units')));
      series = timetable(times);
      observed = timetable(times);
      for m = 1:size(map, 1)
         x = double(ncread(filename, map(m, 1)));
         series.(map(m, 2)) = x(:);
         origin = nan(numel(x), 1);
         try
            origin = double(ncread(filename, map(m, 1) + "_origin"));
         catch
            % No origin flag means no proof that a value is native.
         end
         % A reconstructed sample in one channel does not invalidate a
         % simultaneous observed sample in another channel.
         observed.(map(m, 2)) = origin(:) == 0;
      end
      % A location declared in the file takes precedence. The real
      % Vandecrux surface NetCDFs carry no location attributes, so the
      % station catalog is the fallback for the donor geometry. That
      % catalog comes from the Dataverse metadata of the dataset.
      [~, base] = fileparts(filename);
      station = string(extractBefore(base + "_", "_surface_"));
      try
         location = struct( ...
            'lat_wgs84', double(ncreadatt(filename, '/', 'latitude')), ...
            'lon_wgs84', double(ncreadatt(filename, '/', 'longitude')), ...
            'elev_m', double(ncreadatt(filename, '/', 'elevation')));
      catch
         station_info = icemodel.forcing.helpers ...
            .gcnetVandecruxStationMetadata(station);
         location = station_info.site_location;
      end
      if ~isfinite(location.lat_wgs84) || ~isfinite(location.lon_wgs84)
         error('icemodel:forcing:readGcnetDonor:unknownStation', ...
            'no coordinates for GC-Net station %s', station);
      end
      donor = struct('series', series, ...
         'station', station, ...
         'family', "gcnet", ...
         'location', location, ...
         'observed_mask', observed);
   catch
      % A file without the required channels or coordinates cannot meet
      % the donor contract. Skip it and do not infer the missing metadata.
      donor = [];
   end
end
