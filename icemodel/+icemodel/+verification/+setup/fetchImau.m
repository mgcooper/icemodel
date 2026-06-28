function [source_dir, status] = fetchImau(kwargs)
   %FETCHIMAU Locate or verify local IMAU PANGAEA source caches.
   %
   %  source_dir = icemodel.verification.setup.fetchImau()
   %  [source_dir, status] = icemodel.verification.setup.fetchImau( ...
   %     cache_dir="/path", strict=false)
   %
   % Role
   %  Cache validator for the IMAU hourly S21/S22/S23 network plus the daily
   %  19-site SEB product used for QA/provenance. It does not auto-download in
   %  tests; it reports the DOI landing URLs and validates local files.

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.products (1, :) string = ["hourly", "daily"]
      kwargs.strict (1, 1) logical = true
      kwargs.silent (1, 1) logical = false
   end

   cache_dir = kwargs.cache_dir;
   icemodel.helpers.ensureDirExists(cache_dir);

   % Build status first so callers can use this as a dry-run retrieval manifest.
   status = productStatus(cache_dir, kwargs.products);
   missing = string({status(~[status.present]).product});

   if isempty(missing)
      source_dir = cache_dir;
      return
   end

   % The banner keeps the two PANGAEA products distinct: hourly sites are the
   % case inventory, daily stations are QA/provenance only.
   if ~kwargs.silent
      fprintf('\n');
      fprintf('=== IMAU source cache incomplete ===\n');
      fprintf('Cache directory: %s\n', cache_dir);
      fprintf('Missing products: %s\n', strjoin(missing, ', '));
      fprintf('\nRetrieval:\n');
      for k = 1:numel(status)
         fprintf('  %-6s DOI: %s\n', status(k).product, status(k).doi);
         fprintf('         URL: %s\n', status(k).landing_url);
      end
      fprintf('\nPlace downloaded files under cache_dir/<product>/ or use files\n');
      fprintf('whose names contain the product token.\n\n');
   end

   if kwargs.strict
      error('icemodel:verification:fetchImau:missingSources', ...
         'IMAU source cache incomplete in %s. Missing: %s', ...
         cache_dir, strjoin(missing, ', '));
   end

   source_dir = cache_dir;
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical IMAU source-cache directory.
   pathname = string(fullfile(icemodel.getpath('data'), ...
      'verification', 'imau'));
end

function status = productStatus(cache_dir, products)
   %PRODUCTSTATUS Return one validation/URL row per requested IMAU product.
   status = repmat(productRow("", "", cache_dir), 1, numel(products));
   for k = 1:numel(products)
      status(k) = productRow(products(k), productDoi(products(k)), cache_dir);
   end
end

function row = productRow(product, doi, cache_dir)
   %PRODUCTROW Build one stable status row.
   product_dir = fullfile(cache_dir, product);
   files = [dir(fullfile(product_dir, '*')); ...
      dir(fullfile(cache_dir, "*" + product + "*"))];
   files = files(~[files.isdir]);
   row = struct( ...
      'product', product, ...
      'doi', doi, ...
      'landing_url', "https://doi.org/" + doi, ...
      'present', ~isempty(files), ...
      'cache_dir', string(product_dir));
end

function doi = productDoi(product)
   %PRODUCTDOI Map IMAU products to PANGAEA DOI ids.
   switch string(product)
      case "hourly"
         doi = "10.1594/PANGAEA.971647";
      case "daily"
         doi = "10.1594/PANGAEA.970127";
      otherwise
         error('icemodel:verification:fetchImau:unknownProduct', ...
            'unknown IMAU product: %s', product);
   end
end
