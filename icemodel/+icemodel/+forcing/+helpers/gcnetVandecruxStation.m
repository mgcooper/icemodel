function station = gcnetVandecruxStation(station)
   %GCNETVANDECRUXSTATION Normalize Vandecrux/GC-Net station aliases.
   %
   %  station = icemodel.forcing.helpers.gcnetVandecruxStation("dye2")
   %
   % Fetch validation, inventory discovery, and builders all match the same
   % Vandecrux cache filenames. Keep alias normalization here so every path agrees
   % on the canonical station tokens.

   arguments
      station (1, :) string
   end

   station = reshape(station, 1, []);
   for k = 1:numel(station)
      switch icemodel.forcing.helpers.normalizedFileToken(station(k))
         case {"dye2", "dy2", "dye2long"}
            station(k) = "DYE_2";
         case {"sum", "summit"}
            station(k) = "Summit";
         otherwise
            station(k) = string(station(k));
      end
   end
end
