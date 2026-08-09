function label = safeLabel(value)
   %SAFELABEL Collapse control characters for MATLAB graphics text.
   %
   % A newline or tab in saved metadata would otherwise break an axis label or
   % the markdown line it sits on. This is the one place the control-character
   % class is defined; sanitizeText, escapeMarkdownText, and markdownCode all
   % start here and then add their own encoding.

   label = strtrim(regexprep(string(value), '[\x00-\x1F\x7F]+', ' '));
end
