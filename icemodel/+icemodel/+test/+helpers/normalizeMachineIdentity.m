function identity = normalizeMachineIdentity(raw)
   %NORMALIZEMACHINEIDENTITY Fold a machine name to a stable identity.
   %
   %  identity = icemodel.test.helpers.normalizeMachineIdentity(raw)
   %
   % RAW is machine-name text (char or string), such as a hostname probe
   % result or a baseline's saved meta.hostname value. IDENTITY lowercases
   % the text, strips leading and trailing whitespace, and strips one
   % trailing ".local" suffix.
   %
   % On macOS the `hostname` command returns a DHCP or reverse-DNS name
   % that changes with the network (for example "Mac.lan"), while
   % `scutil --get LocalHostName` and the ".local" Bonjour name stay fixed
   % (DesignSpec decisions 18 and 19). With this normalization a baseline
   % that saved a hostname value (for example "MacBook-Air-2.local")
   % compares equal to the LocalHostName-derived identity (for example
   % "macbook-air-2"), so a rolling baseline built with the hostname probe
   % stays comparable with the LocalHostName probe.
   %
   % A blank RAW returns unchanged (still blank), so a caller's
   % isblanktext check applies to the result.
   %
   % Examples
   %  icemodel.test.helpers.normalizeMachineIdentity("MacBook-Air-2.local")
   %     returns "macbook-air-2"
   %  icemodel.test.helpers.normalizeMachineIdentity("X.LOCAL")
   %     returns "x"
   %  icemodel.test.helpers.normalizeMachineIdentity("a.local.local")
   %     returns "a.local"
   %
   % See also: icemodel.test.helpers.machineHostname,
   %  icemodel.test.helpers.perfBaselineCompatibility, run_aa_acceptance

   % A blank value names no machine. Return it blank rather than raise, so
   % the caller's isblanktext check decides.
   text = strtrim(string(raw));
   if isblanktext(text)
      identity = text;
      return
   end

   % Case does not distinguish machines, so fold it away before the
   % suffix strip.
   text = lower(text);

   % Strip exactly one trailing ".local", the Bonjour/mDNS domain macOS
   % appends to the LocalHostName. An interior ".local" stays:
   % "a.local.local" gives "a.local".
   identity = regexprep(text, '\.local$', '');
end
