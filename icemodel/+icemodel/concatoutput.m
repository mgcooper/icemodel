function [ice1, ice2] = concatoutput(ice1, ice2, thisice1, thisice2)
   %CONCATOUTPUT Concatenate yearly icemodel output structures.
   %
   %  [ice1, ice2] = icemodel.concatoutput(ice1, ice2, thisice1, thisice2)
   %
   % Concatenate one year's output onto an existing output accumulator. ICE1
   % may be either the raw struct returned by the core model or a postprocessed
   % timetable loaded from disk. ICE2 is a struct with time varying fields
   % stored column-wise and static depth fields stored once.

   if isempty(ice1)
      ice1 = thisice1;
      ice2 = thisice2;
      return
   end

   if istimetable(ice1)
      ice1 = [ice1; thisice1];
   else
      fields = fieldnames(ice1);
      for n = 1:numel(fields)
         thisfield = fields{n};
         ice1.(thisfield) = cat(1, ice1.(thisfield), thisice1.(thisfield));
      end
   end

   fields = fieldnames(ice2);
   for n = 1:numel(fields)
      thisfield = fields{n};
      if strcmp(thisfield, 'Z')
         continue
      elseif strcmp(thisfield, 'Time')
         ice2.Time = [ice2.Time; thisice2.Time];
      else
         ice2.(thisfield) = cat(2, ice2.(thisfield), thisice2.(thisfield));
      end
   end
end
