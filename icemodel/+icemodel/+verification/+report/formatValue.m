function text = formatValue(value)
   %FORMATVALUE Format one scalar table value for Markdown.

   if iscell(value)
      value = value{1};
   end
   if isdatetime(value)
      if isnat(value)
         text = "NA";
      else
         text = string(value, "yyyy-MM-dd HH:mm:ss z");
      end
   elseif isnumeric(value)
      % markdownTable slices a whole row of a variable, so a multi-column
      % variable arrives here as a vector. Join it rather than testing a
      % vector in an if, which would error.
      if ~isscalar(value) && ~isempty(value)
         parts = strings(numel(value), 1);
         for k = 1:numel(value)
            if isfinite(value(k))
               parts(k) = string(sprintf('%.5g', value(k)));
            else
               parts(k) = "NA";
            end
         end
         text = strjoin(parts, ", ");
      elseif isempty(value) || ~isfinite(value)
         text = "NA";
      else
         text = string(sprintf('%.5g', value));
      end
   elseif islogical(value)
      text = string(value);
   else
      text = join(string(value), ", ");
   end
   text = icemodel.verification.report.escapeMarkdownText(text);
end
