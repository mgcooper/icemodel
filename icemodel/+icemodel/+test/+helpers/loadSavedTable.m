function T = loadSavedTable(pathname, varnames)
   %LOADSAVEDTABLE Load a saved table-like object from a MAT file.
   %
   %  T = icemodel.test.helpers.loadSavedTable(pathname, varnames)

   % Prefer the requested variable names before falling back to the first
   % saved payload in the MAT file.
   S = load(pathname);
   raw = [];
   for i = 1:numel(varnames)
      if isfield(S, varnames(i))
         raw = S.(varnames(i));
         break
      end
   end
   if isempty(raw)
      fn = fieldnames(S);
      raw = S.(fn{1});
   end

   % Normalize both table and struct payloads into a table for callers.
   if istable(raw)
      T = raw;
   elseif isstruct(raw)
      T = structToTable(raw);
   else
      error('unsupported saved table format: %s', class(raw))
   end
end

function T = structToTable(s)
   %STRUCTTOTABLE Robust struct-to-table conversion for saved test artifacts.

   % Expand scalar structs with vector-valued fields into a row-oriented table.
   if isscalar(s)
      fn = fieldnames(s);
      scalar_fields = true;
      n_rows = [];
      for i = 1:numel(fn)
         value = s.(fn{i});
         if iscell(value) || isnumeric(value) || isstring(value) ...
               || isdatetime(value)
            if isscalar(value)
               continue
            end
            if isempty(n_rows)
               n_rows = numel(value);
            end
            if numel(value) ~= n_rows
               scalar_fields = false;
               break
            end
         else
            scalar_fields = false;
            break
         end
      end
      if scalar_fields && ~isempty(n_rows)
         % Broadcast scalar fields so every output row retains the metadata.
         out = struct();
         for i = 1:numel(fn)
            value = s.(fn{i});
            if isscalar(value)
               out.(fn{i}) = broadcastScalar(value, n_rows);
            else
               out.(fn{i}) = toColumn(value);
            end
         end
         T = struct2table(out);
         return
      end
   end

   % Fall back to MATLAB's default conversion for ordinary struct arrays.
   T = struct2table(s);
end

function value = broadcastScalar(value, n_rows)
   %BROADCASTSCALAR Repeat one scalar value across N_ROWS table rows.

   if isstring(value) || ischar(value)
      value = repmat(string(value), n_rows, 1);
   elseif isnumeric(value) || isdatetime(value)
      value = repmat(value, n_rows, 1);
   else
      value = repmat({value}, n_rows, 1);
   end
end

function value = toColumn(value)
   %TOCOLUMN Normalize one payload into the column shape struct2table expects.

   if ischar(value)
      value = string({value});
   elseif isstring(value) || isnumeric(value) || isdatetime(value) ...
         || iscell(value)
      value = value(:);
   end
end
