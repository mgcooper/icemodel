function tag = sanitizeTag(tag)
   %SANITIZETAG Replace punctuation and whitespace for filename-safe tags.

   % Normalize once so the replacement rules work for char or string input.
   % Keep the transform narrow so saved baseline names stay readable.
   tag = replace(string(tag), ".", "_");
   tag = replace(tag, "-", "_");
   tag = replace(tag, " ", "_");
end
