function donor = readGcnetDonor(filename)
   %READGCNETDONOR Load one GC-Net surface file as an observed-only donor.
   %
   %  donor = icemodel.forcing.helpers.readGcnetDonor(filename)
   %
   % Each canonical channel keeps its own origin mask. Samples without a
   % per-sample origin flag are ineligible because their native provenance
   % cannot be proven.

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
      % for most rows; the shared row-index convention gives the exact
      % hourly axis every consumer (builder and donor alike) agrees on.
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
      % A self-describing artifact's declared location wins; the REAL
      % Vandecrux surface NetCDFs carry no location attributes at all, so
      % the station catalog (fed by the dataset's own Dataverse metadata)
      % is the authoritative fallback for donor geometry.
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
      % Files without required channels or coordinates cannot pass the
      % donor contract; skip them rather than infer missing metadata.
      donor = [];
   end
end
