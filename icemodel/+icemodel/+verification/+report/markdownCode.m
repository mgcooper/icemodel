function text = markdownCode(value)
   %MARKDOWNCODE Wrap saved metadata in a code span that renders literally.
   %
   % A code span renders its contents literally, so the angle brackets are
   % already inert. They are entity-encoded anyway so that no rendered report
   % ever contains a raw <script> sequence. The fence also grows past any
   % backticks in the value so the span cannot be closed early.
   %
   % sanitizeText encodes the ampersand too, so a value containing one shows
   % as &amp; rather than &. That costs display fidelity and buys no safety
   % inside a code span; it is accepted so one sanitizer serves every report
   % writer. A path or tag with a literal ampersand is the
   % case where that shows.

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
