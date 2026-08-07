function fields = cumulative_output_fields()
   %CUMULATIVE_OUTPUT_FIELDS Return cumulative column diagnostic channels.
   %
   % These channels store endpoint state and therefore use the final native
   % sample, rather than a mean, when fixed-step output is aggregated.

   fields = {'melt', 'runoff', 'freeze', 'dlayer'};
end
