function text = markdownCode(value)
   %MARKDOWNCODE Wrap saved metadata in a code span that renders literally.
   %
   % A code span renders its contents literally, so the angle brackets cannot
   % execute. This function entity-encodes them anyway, so no rendered report
   % contains a raw <script> sequence. The fence also grows past any backticks
   % in the value, so the span cannot close early.
   %
   % sanitizeText encodes the ampersand too, so a value that contains one shows
   % as &amp; rather than &. That costs display accuracy and adds no safety
   % inside a code span. It is accepted so that one sanitizer serves every
   % report writer. A path or tag with a literal ampersand shows this effect.

   text = icemodel.verification.report.sanitizeText(value);
   runs = regexp(char(text), '`+', 'match');
   fence_length = 1;
   if ~isempty(runs)
      fence_length = max(cellfun(@numel, runs)) + 1;
   end
   fence = string(repmat('`', 1, fence_length));
   if startsWith(text, "`") || endsWith(text, "`")
      text = fence + " " + text + " " + fence;
   else
      text = fence + text + fence;
   end
end
