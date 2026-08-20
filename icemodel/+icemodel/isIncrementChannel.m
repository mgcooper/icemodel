function tf = isIncrementChannel(names)
   %ISINCREMENTCHANNEL True for per-step increment channels.
   %
   % A channel is an increment channel when its name starts with df_ or
   % when the surface-output namelist declares it additive. Retiming sums
   % these channels over a bin; averaging would divide their totals by the
   % samples per bin.
   %
   % Accepts a char row, a string, or an array of either, and returns a
   % logical of the same shape so both scalar and vectorized callers can use
   % it. retimeHourlyFixedStep documents a generated-code path through
   % this helper, so the lookup loops with strcmp instead of using
   % codegen-unsupported cell set operations.
   %
   % Inputs
   %  names - channel name or array of channel names
   %
   % Outputs
   %  tf - true where the name is a per-step increment channel
   %
   %#codegen

   % Additivity is a property of the channel, owned by the surface-output
   % namelist: 'additive_diagnostic' names the additive solver
   % diagnostics in one place.
   additive = icemodel.namelists.surfaceoutputs('additive_diagnostic');
   tf = startsWith(names, 'df_');
   for k = 1:numel(additive)
      tf = tf | strcmp(names, additive{k});
   end
end
