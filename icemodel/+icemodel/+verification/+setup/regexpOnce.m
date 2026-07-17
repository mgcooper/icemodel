function value = regexpOnce(text, pattern)
   %REGEXPONCE Return one stripped regexp token or an empty string.
   %
   %  value = icemodel.verification.setup.regexpOnce(text, pattern)
   %
   % Setup XML/text metadata readers share the forcing helper so source parsers
   % do not depend on setup orchestration internals.

   value = icemodel.forcing.helpers.regexpOnce(text, pattern);
end
