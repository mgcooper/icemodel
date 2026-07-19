function inventory = retmipOutputInventory(filename)
   %RETMIPOUTPUTINVENTORY Return variables in a RetMIP model-output NetCDF.
   %
   %  inventory = icemodel.verification.setup.retmipOutputInventory(filename)
   %
   % Role
   %  Lightweight model-output indexer. RetMIP model outputs are comparison
   %  products, not forcing, so importers need variable inventory without loading
   %  entire NetCDF arrays.

   arguments
      filename (1, 1) string
   end

   % ncinfo reads headers only, which is enough to validate model-output product
   % contents in tests and manifests.
   info = ncinfo(filename);
   vars = string({info.Variables.Name});
   inventory = struct( ...
      'filename', filename, ...
      'variables', vars(:));
end
