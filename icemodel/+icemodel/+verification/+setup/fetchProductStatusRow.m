function row = fetchProductStatusRow(cache_dir, product, doi, patterns, kwargs)
   %FETCHPRODUCTSTATUSROW Build the standard fetch product status record.
   %
   % Fetchers use this for product rows whose presence is determined by one or
   % more local file patterns. Dataset-specific fetchers may add fields after
   % construction when they need richer provenance.

   arguments
      cache_dir (1, 1) string
      product (1, 1) string
      doi (1, 1) string
      patterns (1, :) string
      kwargs.cache_subdir (1, 1) string = product
      kwargs.exclude_folders (1, :) string = strings(1, 0)
      kwargs.exclude_names (1, :) string = strings(1, 0)
   end

   files = icemodel.verification.setup.fetchProductFiles( ...
      cache_dir, patterns, exclude_folders=kwargs.exclude_folders, ...
      exclude_names=kwargs.exclude_names);
   row = struct( ...
      'product', product, ...
      'doi', doi, ...
      'landing_url', "https://doi.org/" + doi, ...
      'present', ~isempty(files), ...
      'cache_dir', string(fullfile(cache_dir, kwargs.cache_subdir)), ...
      'found_files', files);
end
