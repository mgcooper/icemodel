function lines = markdownTable(values)
   %MARKDOWNTABLE Convert a compact table to inert Markdown.

   vars = string(values.Properties.VariableNames);
   if isempty(vars)
      lines = "No tabular fields were saved.";
      return
   end
   header = "| " + join(replace(vars, "_", " "), " | ") + " |";
   divider = "| " + join(repmat("---", size(vars)), " | ") + " |";
   lines = strings(height(values) + 2, 1);
   lines(1:2) = [header; divider];
   for row = 1:height(values)
      cells = strings(size(vars));
      for col = 1:numel(vars)
         column = values.(vars(col));
         cells(col) = icemodel.verification.report.formatValue(column(row, :));
      end
      lines(row + 2) = "| " + join(cells, " | ") + " |";
   end
end
