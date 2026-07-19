function status = buildFetchProductStatus(cache_dir, products, registry, ...
      row_builder, kwargs)
   %BUILDFETCHPRODUCTSTATUS Build ordered status rows from a product registry.
   %
   %  status = icemodel.verification.setup.buildFetchProductStatus( ...
   %     cache_dir, products, registry, row_builder)
   %
   % The shared loop owns registry lookup and preallocation. The supplied row
   % builder retains only family-specific file matching and provenance fields.

   arguments
      cache_dir (1, 1) string
      products (1, :) string
      registry (:, 2) string
      row_builder (1, 1) function_handle
      kwargs.prototype (1, 1) struct = ...
         icemodel.verification.setup.emptyFetchProductStatusRow()
   end

   status = repmat(kwargs.prototype, 1, numel(products));
   for k = 1:numel(products)
      product = products(k);
      doi = registry(registry(:, 1) == product, 2);
      status(k) = row_builder(cache_dir, product, doi);
   end
end
