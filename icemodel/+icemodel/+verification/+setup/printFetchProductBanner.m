function printFetchProductBanner(cache_dir, status, missing, source_label)
   %PRINTFETCHPRODUCTBANNER Print shared manual-cache retrieval instructions.
   %
   % Family fetchers supply only their display label; status rows already carry
   % the product names, DOI provenance, and landing URLs needed by the banner.

   arguments
      cache_dir (1, 1) string
      status (1, :) struct
      missing (1, :) string
      source_label (1, 1) string
   end

   product_width = max([strlength(string({status.product})), 1]);
   fprintf('\n');
   fprintf('=== %s source cache incomplete ===\n', source_label);
   fprintf('Cache directory: %s\n', cache_dir);
   fprintf('Missing products: %s\n', strjoin(missing, ', '));
   fprintf('\nRetrieval:\n');
   for k = 1:numel(status)
      fprintf('  %-*s DOI: %s\n', product_width, ...
         status(k).product, status(k).doi);
      fprintf('  %*s URL: %s\n', product_width, "", status(k).landing_url);
   end
   fprintf('\nPlace downloaded files under cache_dir/<product>/ or use files\n');
   fprintf('whose names contain the product token.\n\n');
end
