function tf = isIncrementChannel(names)
   %ISINCREMENTCHANNEL True for per-step increment channels.
   %
   % Use this to determine which variables should be retimed using sums over
   % time bins rather than averaging, which would divide their totals by the
   % samples per bin. A channel is an increment channel when its name starts
   % with df_ or when the surface-output namelist declares it as additive.
   %
   % Accepts a char row, a string, or an array of either, and returns a logical
   % of the same shape so both scalar and vectorized callers can use it.
   % retimeHourlyFixedStep includes a generated-code option that calls this
   % helper, so the lookup uses strcmp instead of cell set comparisons which are
   % not supported by codegen.
   %
   % Inputs
   %  names - channel name or array of channel names
   %
   % Outputs
   %  tf - true where the name is a per-step increment channel
   %
   % See also: icemodel.postprocess, icemodel.retimeHourlyFixedStep
   %
   %#codegen

   % Get the 'additivity' property from the surface-output namelist.
   additive = icemodel.namelists.surfaceoutputs('additive_diagnostic');
   tf = startsWith(names, 'df_');
   for k = 1:numel(additive)
      tf = tf | strcmp(names, additive{k});
   end
end
