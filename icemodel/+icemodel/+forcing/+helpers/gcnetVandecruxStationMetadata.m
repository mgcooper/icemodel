function info = gcnetVandecruxStationMetadata(station)
   %GCNETVANDECRUXSTATIONMETADATA Return station aliases and location.
   %
   %  info = icemodel.forcing.helpers.gcnetVandecruxStationMetadata("dye2")
   %
   % Fetch, inventory, and the builders share this metadata. They therefore use
   % the same aliases and station coordinates for discovery and for the staged
   % artifacts.
   %
   % See also: icemodel.forcing.buildGcnetVandecruxData,
   %  icemodel.forcing.helpers.readGcnetDonor
   arguments
      station (1, :) string
   end

   station = icemodel.forcing.helpers.gcnetVandecruxStation(station);
   proto = struct('station', "", 'aliases', strings(1, 0), ...
      'site_location', struct('lat_wgs84', NaN, 'lon_wgs84', NaN, ...
      'elev_m', NaN));
   info = repmat(proto, 1, numel(station));
   % The catalog supplies the identity and the coordinates. The catalog comes
   % from the dataset's own Dataverse metadata. An unknown station keeps the
   % NaN prototype, so callers can detect it.
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
