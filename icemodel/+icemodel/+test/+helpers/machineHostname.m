function hostname = machineHostname(probe)
   %MACHINEHOSTNAME Return the current machine's trimmed hostname.
   %
   %  hostname = icemodel.test.helpers.machineHostname()
   %  hostname = icemodel.test.helpers.machineHostname(probe)
   %
   % PROBE is an optional function handle with the same two outputs as
   % system('hostname').
   %
   % See also: icemodel.test.helpers.perfBaselineCompatibility,
   %  run_perf_suite, build_perf_baseline

   if nargin < 1
      probe = @() system('hostname');
   end

   % computer() identifies the platform architecture, such as MACA64.
   % hostname identifies the machine whose timings are saved.
   [status, value] = probe();
   hostname = string(strtrim(value));
   if status ~= 0 || isblanktext(hostname)
      error('icemodel:test:perf:machineIdentityUnavailable', ...
         'Cannot determine the machine hostname.')
   end
end
