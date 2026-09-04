function station = gcnetVandecruxStation(station)
   %GCNETVANDECRUXSTATION Normalize Vandecrux/GC-Net station aliases.
   %
   %  station = icemodel.forcing.helpers.gcnetVandecruxStation("dye2")
   %
   % Fetch validation, inventory discovery, builders, and the donor
   % self-exclusion all match the same Vandecrux station names. The
   % self-exclusion applies POLICY A8: the target is never its own donor,
   % through every alias spelling. This function reads the alias table from
   % the catalog, so every listed alias resolves to the station's official
   % name. An unknown name passes through unchanged.
   %
   % See also: icemodel.forcing.buildGcnetVandecruxData,
   %  icemodel.forcing.reconstruct.fillPromiceStation

   arguments
      station (1, :) string
   end

   catalog = icemodel.forcing.helpers.gcnetVandecruxCatalog();
   station = reshape(station, 1, []);
   for k = 1:numel(station)
      token = icemodel.forcing.helpers.normalizedFileToken(station(k));
      for c = 1:numel(catalog)
         if ismember(token, ...
               icemodel.forcing.helpers.normalizedFileToken( ...
               catalog(c).aliases))
            station(k) = catalog(c).station;
            break
         end
      end
   end
end
