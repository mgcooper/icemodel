function text = escapeMarkdownText(value)
   %ESCAPEMARKDOWNTEXT Show saved text literally, without markup or raw HTML.
   %
   % Two layers. sanitizeText encodes the HTML delimiters first, so no rendered
   % document ever contains a literal <script sequence even before Markdown
   % runs. Then every remaining ASCII punctuation character is backslash
   % escaped so the value cannot introduce emphasis, a link, or a table cell.
   %
   % This function leaves the ampersand and the semicolon unescaped.
   % Otherwise it would break apart the entities that sanitizeText produced,
   % and the output would show \&amp\; instead of the original character.
   %
   % This function escapes the backslash first. Otherwise it would escape the
   % escapes that it adds afterwards.

   text = icemodel.verification.report.sanitizeText(value);
   punctuation = setdiff([92, 33:47, 58:64, 91, 93:96, 123:126], ...
      [38, 59], 'stable');
   escape = string(char(92));
   for k = 1:numel(punctuation)
      token = string(char(punctuation(k)));
      text = replace(text, token, escape + token);
   end
end
