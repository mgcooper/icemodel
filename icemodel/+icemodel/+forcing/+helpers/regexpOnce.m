function value = regexpOnce(text, pattern)
   %REGEXPONCE Return one stripped regexp token or an empty string.
   %
   %  value = icemodel.forcing.helpers.regexpOnce(text, pattern)
   %
   % Source parsers use this for lightweight text/XML metadata extraction.

   token = regexp(char(text), pattern, 'tokens', 'once');
   if isempty(token)
      value = "";
   else
      value = regexprep(string(token{1}), '<[^>]+>', ' ');
      value = strip(regexprep(value, '\s+', ' '));
   end
end
