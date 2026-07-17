function token = normalizedFileToken(value)
   %NORMALIZEDFILETOKEN Compare filenames case-insensitively across separators.
   token = lower(string(value));
   token = regexprep(token, '[\s_\-]', '');
end
