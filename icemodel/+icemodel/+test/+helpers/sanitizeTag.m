function tag = sanitizeTag(tag)
%SANITIZETAG Replace punctuation and whitespace for filename-safe tags.

   tag = replace(string(tag), ".", "_");
   tag = replace(tag, "-", "_");
   tag = replace(tag, " ", "_");
end
