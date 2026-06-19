function v = slabMean(field, start, count, inslab)
   %SLABMEAN Mean of a static grid field over hyperslab target cells.
   %
   %  v = icemodel.forcing.helpers.slabMean(field, start, count, inslab)
   %
   % Companion to icemodel.forcing.helpers.gridLocation: averages a
   % static 2-D grid field (latitude, elevation, ...) over the
   % in-target cells of the hyperslab gridLocation selected.
   %
   % See also: icemodel.forcing.helpers.gridLocation

   arguments
      field (:, :) double
      start (1, 2) double
      count (1, 2) double
      inslab (:, 1) double
   end

   sub = field(start(1):start(1) + count(1) - 1, ...
      start(2):start(2) + count(2) - 1);
   v = mean(sub(inslab));
end
