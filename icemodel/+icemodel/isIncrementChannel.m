function tf = isIncrementChannel(names)
   %ISINCREMENTCHANNEL True for per-step increment channels.
   %
   % One place that defines the df_ prefix rule. A df_ channel holds
   % one forcing step's change, so retiming sums it over a bin; averaging an
   % increment would divide it by the samples per bin.
   %
   % startsWith is the rule, and it lives here rather than in each consumer so
   % retiming and postprocessing cannot drift apart on what counts.
   %
   % Accepts a char row, a string, or an array of either, and returns a
   % logical of the same shape so both scalar and vectorized callers can use
   % it. Uses startsWith rather than cell set operations to stay
   % codegen-compatible.
   %
   % Inputs
   %  names - channel name or array of channel names
   %
   % Outputs
   %  tf - true where the name is a per-step increment channel
   %
   %#codegen

   tf = startsWith(names, 'df_');
end
