function tf = isIncrementChannel(names)
   %ISINCREMENTCHANNEL True for per-step increment channels.
   %
   % A channel is an increment channel when its name starts with df_. The
   % per-forcing-step coupling recovery count is also additive. Retiming sums
   % these channels over a bin; averaging would divide their totals by the
   % samples per bin.
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

   recovery_count_field = icemodel.namelists.surfaceoutputs( ...
      'icemodel_diagnostic_suffix');
   tf = startsWith(names, 'df_') | strcmp(names, recovery_count_field{1});
end
