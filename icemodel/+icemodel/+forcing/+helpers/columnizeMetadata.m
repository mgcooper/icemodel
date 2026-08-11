function value = columnizeMetadata(value)
   %COLUMNIZEMETADATA Store metadata vectors as columns for inspection.
   %
   %  value = icemodel.forcing.helpers.columnizeMetadata(value)
   %
   % Users inspect metadata saved into MAT files interactively. MATLAB displays
   % Nx1 string, cell, and numeric vectors more readably than 1xN rows. This
   % helper therefore turns each single-row multi-element metadata vector into a
   % column. It works recursively and leaves char text and true matrices
   % unchanged.

   if istimetable(value)
      return
   end

   if istable(value)
      if height(value) == 1 && width(value) > 1
         value = rowTableToColumnTable(value);
      end
      return
   end

   if isstruct(value)
      for n = 1:numel(value)
         names = fieldnames(value(n));
         for k = 1:numel(names)
            field = names{k};
            value(n).(field) = icemodel.forcing.helpers.columnizeMetadata( ...
               value(n).(field));
         end
      end
      if isrow(value) && ~isscalar(value)
         value = value(:);
      end
      return
   end

   if iscell(value)
      for k = 1:numel(value)
         value{k} = icemodel.forcing.helpers.columnizeMetadata(value{k});
      end
   end

   % Do not transpose character text; a char row is one scalar string in the
   % older MATLAB representation, not a metadata vector.
   if ischar(value)
      return
   end

   if isrow(value) && numel(value) > 1
      value = value(:);
   end
end

function value = rowTableToColumnTable(value)
   %ROWTABLETOCOLUMNTABLE Make one-row metadata tables inspectable.
   names = string(value.Properties.VariableNames(:));
   cells = table2cell(value);
   cells = cells(:);
   if all(cellfun(@(x) isnumeric(x) && isscalar(x), cells))
      values = cellfun(@double, cells);
   elseif all(cellfun(@(x) islogical(x) && isscalar(x), cells))
      values = cellfun(@logical, cells);
   else
      values = cells;
   end
   value = table(names, values, 'VariableNames', {'variable', 'value'});
end
