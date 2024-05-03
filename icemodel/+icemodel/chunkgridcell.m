function [S, E] = chunkgridcell(firstcell, finalcell, numjob, whichjob)

   % Decided not to complete this function, but this is how it would have been
   % called from rungrid_1:
   %
   % firstcell = 1;
   % finalcell = numel(gridnums);
   % numjob = 2;
   % whichjob = 1;
   % [~, E] = icemodel.chunkgridcell(firstcell, finalcell, numjob, whichjob);
   
   numcell = (finalcell - firstcell) + 1;

   % This replicates how I divided icemodel/mar and icemodel/modis into two
   % equal chunks each:
   switch whichjob
      case 1
         S = 1;
         E = floor(numcell / numjob);
      case 2
         S = ceil(numcell / 2);
         E = numcell;
   end


   % % Simplest case:
   % S = 1;
   % E = numel(gridnums);

   % % Divide into N ~equal chunks of size ni and startpoint si:
   % N = 1;
   % si = 1;
   % ni = 154;
   %
   % S = si + ni * (N-1);
   % E = si + ni * (N) - 1;

   % % Divide into two equal chunks for smbmodel/mar
   % S = 1;
   % E = floor(numel(gridnums) / 2);
end
