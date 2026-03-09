function T = loadSavedTable(pathname, varnames)
%LOADSAVEDTABLE Load a saved table/struct from a MAT file.

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

   if istable(raw)
      T = raw;
   elseif isstruct(raw)
      T = test.helpers.structToTable(raw);
   else
      error('unsupported saved table format: %s', class(raw))
   end
end
