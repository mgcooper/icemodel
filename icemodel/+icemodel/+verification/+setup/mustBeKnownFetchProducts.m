function mustBeKnownFetchProducts(products, registry, error_id, error_label)
   %MUSTBEKNOWNFETCHPRODUCTS Reject selectors absent from a fetch registry.
   %
   % Family fetchers supply their stable error id and display label so the
   % shared validation preserves each public error contract.

   arguments
      products (1, :) string
      registry (:, 2) string
      error_id (1, 1) string
      error_label (1, 1) string
   end

   names = icemodel.verification.setup.fetchProductNames(registry);
   bad = setdiff(products, names, 'stable');
   if ~isempty(bad)
      error(error_id, 'unknown %s product(s): %s', ...
         error_label, strjoin(bad, ', '));
   end
end
