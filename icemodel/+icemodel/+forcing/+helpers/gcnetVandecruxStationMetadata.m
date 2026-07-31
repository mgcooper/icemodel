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
   % Identity and coordinates come from the single-source catalog (fed by
   % the dataset's own Dataverse metadata); unknown stations keep the NaN
   % prototype so callers can detect them.
   catalog = icemodel.forcing.helpers.gcnetVandecruxCatalog();
   names = string({catalog.station});
   for k = 1:numel(station)
      info(k).station = station(k);
      info(k).aliases = station(k);
      match = find(names == station(k), 1);
      if ~isempty(match)
         info(k).aliases = catalog(match).aliases;
         info(k).site_location = catalog(match).site_location;
      end
      info(k).site_location = icemodel.forcing.helpers.projectLocation( ...
         info(k).site_location);
   end
end
