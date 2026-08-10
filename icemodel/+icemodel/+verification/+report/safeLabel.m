function label = safeLabel(value)
   %SAFELABEL Collapse control characters for MATLAB graphics text.
   %
   % Replaces each run of control characters (0x00-0x1F and 0x7F) with a
   % single space and trims the result. A newline or tab in saved metadata
   % breaks an axis label or the markdown line it sits on. sanitizeText,
   % escapeMarkdownText, and markdownCode call this first, then apply their
   % own encoding.

   label = strtrim(regexprep(string(value), '[\x00-\x1F\x7F]+', ' '));
end
