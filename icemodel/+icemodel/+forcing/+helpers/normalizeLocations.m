function [locations, batch] = normalizeLocations(location)
   %NORMALIZELOCATIONS Normalize one location or a point list to a row cell.
   %
   %  [locations, batch] = ...
   %     icemodel.forcing.helpers.normalizeLocations(location)
   %
   % A polyshape or one [lat lon] row remains one scalar location. An N-by-2
   % numeric point list with N greater than one becomes a 1-by-N cell array and
   % sets BATCH true so gridded builders can preserve their public return shape.

   % Only a multirow numeric point list activates the batch contract; every
   % other supported location object remains one opaque location value.
   if isnumeric(location) && size(location, 2) == 2 && size(location, 1) > 1
      locations = num2cell(location, 2)';
      batch = true;
   else
      locations = {location};
      batch = false;
   end
end
