function colocation = mergeColocation(colocation, add)
   %MERGECOLOCATION Copy every field from ADD onto a colocation struct.
   %
   % The dataset-family importers call this helper after each staged source,
   % so that one function performs the manifest-leg merge.

   fields = fieldnames(add);
   for k = 1:numel(fields)
      colocation.(fields{k}) = add.(fields{k});
   end
end
