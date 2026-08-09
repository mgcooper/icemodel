function text = sanitizeText(value)
   %SANITIZETEXT Collapse controls and neutralize raw HTML delimiters.

   text = icemodel.verification.report.safeLabel(value);
   text = replace(text, "&", "&amp;");
   text = replace(text, "<", "&lt;");
   text = replace(text, ">", "&gt;");
end
