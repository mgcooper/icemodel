function time = gcnetHourlyAxis(time)
   %GCNETHOURLYAXIS Use the documented hourly row-index time convention.
   %
   %  time = icemodel.forcing.helpers.gcnetHourlyAxis(time)
   %
   % Role
   %  The real Vandecrux surface files encode time as fractional days, but
   %  the coordinate drifts inside leap years and jumps at several year
   %  boundaries while the row count and endpoint still describe a
   %  continuous hourly series. The product is hourly, so every consumer
   %  (met building and donor reading alike) uses the unambiguous
   %  row-index axis anchored at the first timestamp. Shared here so the
   %  builder and the donor reader can never disagree about the axis.
   %
   % Returns
   %  time : column datetime, strictly hourly from time(1).
   %
   % See also: icemodel.forcing.buildGcnetVandecruxData,
   %  icemodel.forcing.helpers.readGcnetDonor

   time = time(:);
   if numel(time) < 2
      return
   end
   idx = 0:numel(time) - 1;
   time = time(1) + hours(idx(:));
end
