function station = gcnetVandecruxStation(station)
   %GCNETVANDECRUXSTATION Normalize Vandecrux/GC-Net station aliases.
   %
   %  station = icemodel.forcing.helpers.gcnetVandecruxStation("dye2")
   %
   % Fetch validation, inventory discovery, builders, and the donor
   % self-exclusion (POLICY A8: the target is never its own donor,
   % through every alias spelling) all match the same Vandecrux station
   % names. Canonicalization derives from the single-source catalog's
   % alias table so no spelling can escape it; unknown names pass
   % through unchanged.

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
