function [source_dir, status] = fetchRetmip(kwargs)
   %FETCHRETMIP Locate or verify local RetMIP source caches.
   %
   %  source_dir = icemodel.verification.setup.fetchRetmip()
   %  [source_dir, status] = icemodel.verification.setup.fetchRetmip( ...
   %     cache_dir="/path", strict=false)
   %
   % Role
   %  Cache validator and retrieval-instruction helper. It validates local
   %  RetMIP forcing/protocol, model-output, and script caches without pulling
   %  multi-GB products during tests.

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.products (1, :) string = ["forcing", "outputs", "scripts"]
      kwargs.strict (1, 1) logical = true
      kwargs.silent (1, 1) logical = false
   end

   cache_dir = kwargs.cache_dir;
   icemodel.helpers.ensureDirExists(cache_dir);

   % Build status first so callers/tests can inspect DOI/URL construction even
   % when the local cache is intentionally incomplete.
   status = productStatus(cache_dir, kwargs.products);
   missing = string({status(~[status.present]).product});

   if isempty(missing)
      source_dir = cache_dir;
      return
   end

   % Print actionable instructions instead of hiding the missing cache behind a
   % generic file-not-found error.
   if ~kwargs.silent
      fprintf('\n');
      fprintf('=== RetMIP source cache incomplete ===\n');
      fprintf('Cache directory: %s\n', cache_dir);
      fprintf('Missing products: %s\n', strjoin(missing, ', '));
      fprintf('\nRetrieval:\n');
      for k = 1:numel(status)
         fprintf('  %-8s DOI: %s\n', status(k).product, status(k).doi);
         fprintf('           URL: %s\n', status(k).landing_url);
      end
      fprintf('\nPlace downloaded files under cache_dir/<product>/ or use files\n');
      fprintf('whose names contain the product token.\n\n');
   end

   if kwargs.strict
      error('icemodel:verification:fetchRetmip:missingSources', ...
         'RetMIP source cache incomplete in %s. Missing: %s', ...
         cache_dir, strjoin(missing, ', '));
   end

   source_dir = cache_dir;
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical RetMIP source-cache directory.
   pathname = string(fullfile(icemodel.getpath('data'), ...
      'verification', 'retmip'));
end

function status = productStatus(cache_dir, products)
   %PRODUCTSTATUS Return one validation/URL row per requested RetMIP product.
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
      'dataverse_api_url', ...
         "https://dataverse.geus.dk/api/datasets/:persistentId/?persistentId=doi:" + doi, ...
      'present', ~isempty(files), ...
      'cache_dir', string(product_dir));
end

function doi = productDoi(product)
   %PRODUCTDOI Map RetMIP products to their GEUS Dataverse DOI.
   switch string(product)
      case {"forcing", "scripts"}
         doi = "10.22008/FK2/GZ3CSN";
      case "outputs"
         doi = "10.22008/FK2/CVPUJL";
      otherwise
         error('icemodel:verification:fetchRetmip:unknownProduct', ...
            'unknown RetMIP product: %s', product);
   end
end
