function fields = cumulativeoutputs()
   %CUMULATIVEOUTPUTS Return cumulative column diagnostic channels.
   %
   % These channels store a running total, so fixed-step aggregation takes the
   % final native sample rather than a mean. Note runoff is a balance rather
   % than a strict cumulative sum: it can fall when evaporation exceeds the
   % melt and overflow arriving in the same step.

   fields = {'melt', 'runoff', 'freeze', 'dlayer'};
end
