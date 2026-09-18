function identity = machineHostname(probe, is_mac)
   %MACHINEHOSTNAME Return a normalized, network-independent machine identity.
   %
   %  identity = icemodel.test.helpers.machineHostname()
   %  identity = icemodel.test.helpers.machineHostname(probe)
   %  identity = icemodel.test.helpers.machineHostname(probe, is_mac)
   %
   % On macOS the `hostname` command returns a DHCP or reverse-DNS name that
   % changes with the network (for example "Mac.lan" on one network and
   % "MacBook-Air-2.local" on another), so a saved machine identity must not
   % depend on it (DesignSpec decisions 18 and 19). On macOS this function
   % probes `scutil --get LocalHostName` first; that name and the ".local"
   % Bonjour name stay fixed across networks. A nonzero status or a blank
   % result falls back to `hostname`. Other platforms use `hostname` only.
   % Either raw value is normalized by
   % icemodel.test.helpers.normalizeMachineIdentity, so a baseline that saved
   % a hostname value (for example "MacBook-Air-2.local") compares equal to
   % the LocalHostName-derived identity (for example
   % "macbook-air-2").
   %
   % PROBE is an optional function handle called as
   % [status, value] = probe(command), where COMMAND is the shell command
   % text ("scutil --get LocalHostName" or "hostname") and the two outputs
   % match system(command). A test injects each platform branch through
   % it: macOS with scutil success, macOS with scutil failure and the
   % hostname fallback, non-macOS with hostname only, and a blank or failed
   % result from every probe. IS_MAC is an optional logical that selects
   % the macOS branch.
   % It defaults to ismac() and lets a test simulate either platform
   % without depending on the host operating system.
   %
   % When no probe returns a usable name, the function raises
   % icemodel:test:perf:machineIdentityUnavailable.
   %
   % See also: icemodel.test.helpers.normalizeMachineIdentity,
   %  icemodel.test.helpers.perfBaselineCompatibility, run_perf_suite,
   %  build_perf_baseline, run_aa_acceptance

   if nargin < 1 || isempty(probe)
      probe = @(command) system(command);
   end
   if nargin < 2 || isempty(is_mac)
      is_mac = ismac();
   end

   % On macOS, prefer the Bonjour/mDNS LocalHostName: unlike `hostname`,
   % it does not change when the machine joins a different network.
   raw = "";
   if is_mac
      [status, value] = probe('scutil --get LocalHostName');
      candidate = icemodel.test.helpers.normalizeMachineIdentity(value);
      if status == 0 && ~isblanktext(candidate)
         raw = string(strtrim(value));
      end
   end

   % Fall back to `hostname` when scutil is unavailable, fails, returns
   % blank, or the platform is not macOS.
   if isblanktext(raw)
      [status, value] = probe('hostname');
      value = string(strtrim(value));
      if status ~= 0 || isblanktext(value)
         error('icemodel:test:perf:machineIdentityUnavailable', ...
            'Cannot determine the machine hostname.')
      end
      raw = value;
   end

   identity = icemodel.test.helpers.normalizeMachineIdentity(raw);
end
