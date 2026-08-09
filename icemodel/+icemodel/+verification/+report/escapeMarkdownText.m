function text = escapeMarkdownText(value)
   %ESCAPEMARKDOWNTEXT Show saved text literally, without markup or raw HTML.
   %
   % Two layers. sanitizeText encodes the HTML delimiters first, so no rendered
   % document ever contains a literal <script sequence even before Markdown
   % runs. Then every remaining ASCII punctuation character is backslash
   % escaped so the value cannot introduce emphasis, a link, or a table cell.
   %
   % Ampersand and semicolon are left unescaped, or the entities sanitizeText
   % just produced would be broken back apart and shown as \&amp\; instead of
   % rendering as the original character.
   %
   % Backslash is escaped first, or the escapes added afterwards would
   % themselves be escaped.

   text = icemodel.verification.report.sanitizeText(value);
   punctuation = setdiff([92, 33:47, 58:64, 91, 93:96, 123:126], ...
      [38, 59], 'stable');
   escape = string(char(92));
   for k = 1:numel(punctuation)
      token = string(char(punctuation(k)));
      text = replace(text, token, escape + token);
   end
end
