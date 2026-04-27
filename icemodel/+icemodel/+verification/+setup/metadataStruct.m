function metadata = metadataStruct(field_pairs)
   %METADATASTRUCT Build a metadata struct from a 2-column cell array.
   %
   %  metadata = icemodel.verification.setup.metadataStruct(field_pairs)
   %
   % Inputs
   %  field_pairs   N-by-2 cell array whose first column is field names and
   %                second column is field values.
   %
   % Outputs
   %  metadata      Struct with one field for each name/value pair.
   %
   % Role
   %  Setup helper shared by importers when assembling artifact metadata and
   %  manifest notes.

   names = string(field_pairs(:, 1));
   values = field_pairs(:, 2);
   metadata = cell2struct(values, names, 1);
end
