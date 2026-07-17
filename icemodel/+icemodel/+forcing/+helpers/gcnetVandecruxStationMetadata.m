function info = gcnetVandecruxStationMetadata(station)
   %GCNETVANDECRUXSTATIONMETADATA Return canonical station aliases and location.
   %
   %  info = icemodel.forcing.helpers.gcnetVandecruxStationMetadata("dye2")
   %
   % Fetch, inventory, and builders share this metadata so aliases and station
   % coordinates cannot drift between discovery and staged artifacts.
   arguments
      station (1, :) string
   end

   station = icemodel.forcing.helpers.gcnetVandecruxStation(station);
   proto = struct('station', "", 'aliases', strings(1, 0), ...
      'site_location', struct('lat_wgs84', NaN, 'lon_wgs84', NaN, ...
      'elev_m', NaN));
   info = repmat(proto, 1, numel(station));
   for k = 1:numel(station)
      info(k).station = station(k);
      switch station(k)
         case "DYE_2"
            info(k).aliases = ["DYE_2", "Dye_2", "DY2", ...
               "dye2_long", "dye2"];
            info(k).site_location = struct('lat_wgs84', 66.48, ...
               'lon_wgs84', -46.28, 'elev_m', 2165);
         case "Summit"
            info(k).aliases = ["Summit", "SUM", "sum", "summit"];
            info(k).site_location = struct('lat_wgs84', 72.58, ...
               'lon_wgs84', -38.50, 'elev_m', 3254);
      otherwise
         info(k).aliases = station(k);
      end
      info(k).site_location = icemodel.forcing.helpers.projectLocation( ...
         info(k).site_location);
   end
end
