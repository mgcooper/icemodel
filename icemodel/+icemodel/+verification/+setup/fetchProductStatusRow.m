function row = fetchProductStatusRow(cache_dir, product, doi, patterns, kwargs)
   %FETCHPRODUCTSTATUSROW Build the standard fetch product status record.
   %
   % Fetchers use this function for product rows. One or more local file
   % patterns decide whether a product is present. A dataset-specific fetcher
   % can add fields to the row when it needs more provenance.

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
