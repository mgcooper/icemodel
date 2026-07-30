function colocation = mergeColocation(colocation, add)
   %MERGECOLOCATION Copy every field from ADD onto a colocation struct.
   %
   % This tiny helper centralizes the manifest-leg merge used by dataset-family
   % importers after each staged source.

   fields = fieldnames(add);
   for k = 1:numel(fields)
      colocation.(fields{k}) = add.(fields{k});
   end
end
