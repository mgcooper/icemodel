function quoted = shellQuote(value)
   %SHELLQUOTE Quote one value as a literal POSIX shell argument.
   quote = char(39);
   escaped = strrep(char(value), quote, [quote '"' quote '"' quote]);
   quoted = string([quote escaped quote]);
end
