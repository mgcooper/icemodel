function quoted = shellQuote(value)
   %SHELLQUOTE Quote one value as a literal argument for the host shell.
   %
   %  quoted = icemodel.shellQuote(value)
   %
   % Anything that builds a command string for system() needs this, so it sits
   % at the runtime level rather than under one consumer. Report rendering,
   % fixture packing, and release preflight all call it.
   %
   % On POSIX the value is wrapped in single quotes and any embedded single
   % quote is closed, escaped, and reopened, so no character is interpreted.
   % On Windows double quotes are the only option, and cmd still expands %VAR%
   % inside them, so a percent sign is rejected rather than silently expanded.
   %
   % Inputs
   %  value - path or argument to quote
   %
   % Outputs
   %  quoted - the quoted argument, as a string

   value = char(string(value));
   if ispc
      if contains(value, '%')
         error('icemodel:shellQuote:unsafeWindowsPath', ...
            'Shell command arguments cannot contain %% on Windows')
      end
      quoted = string([char(34), value, char(34)]);
      return
   end
   quote = char(39);
   escaped = strrep(value, quote, [quote, '"', quote, '"', quote]);
   quoted = string([quote, escaped, quote]);
end
