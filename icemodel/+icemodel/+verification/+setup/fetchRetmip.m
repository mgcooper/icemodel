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
   %  An empty products selection returns an empty status without creating or
   %  scanning cache_dir.

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.products (1, :) string = ...
         icemodel.verification.setup.fetchProductNames(productRegistry())
      kwargs.strict (1, 1) logical = true
      kwargs.silent (1, 1) logical = false
      kwargs.create_cache_dir (1, 1) logical = true
   end

   % Validate the registry selection before cache creation or discovery.
   icemodel.verification.setup.mustBeKnownFetchProducts( ...
      kwargs.products, productRegistry(), ...
      "icemodel:verification:fetchRetmip:unknownProduct", "RetMIP");

   % Resolve and optionally create the family cache root.
   cache_dir = icemodel.verification.setup.resolveFetchCacheDir( ...
      kwargs.cache_dir, defaultCacheDir());
   if kwargs.create_cache_dir && ~isempty(kwargs.products)
      icemodel.helpers.ensureDirExists(cache_dir);
   end

   % Build family-specific cache rows before shared strict handling.
   status = productStatus(cache_dir, kwargs.products);
   source_dir = icemodel.verification.setup.finishFetchStatus( ...
      cache_dir, status, strict=kwargs.strict, silent=kwargs.silent, ...
      error_id="icemodel:verification:fetchRetmip:missingSources", ...
      error_label="RetMIP", banner_callback=@(cache_dir, status, missing) ...
      icemodel.verification.setup.printFetchProductBanner( ...
      cache_dir, status, missing, "RetMIP"));
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical RetMIP source-cache directory.
   pathname = icemodel.forcing.helpers.verificationSourceDir("", "retmip");
end

function status = productStatus(cache_dir, products)
   %PRODUCTSTATUS Return one validation/URL row per requested RetMIP product.
   prototype = icemodel.verification.setup.emptyFetchProductStatusRow();
   prototype.dataverse_api_url = "";
   status = icemodel.verification.setup.buildFetchProductStatus( ...
      cache_dir, products, productRegistry(), @productRow, ...
      prototype=prototype);
end

function row = productRow(cache_dir, product, doi)
   %PRODUCTROW Build one stable status row.
   patterns = [fullfile(product, "*"), "*" + product + "*"];
   if string(product) == "outputs"
      patterns = [patterns, "*_*_3hourly_*.nc"];
   end
   row = icemodel.verification.setup.fetchProductStatusRow( ...
      cache_dir, product, doi, patterns);
   row.dataverse_api_url = ...
      "https://dataverse.geus.dk/api/datasets/:persistentId/?persistentId=doi:" + doi;
end

function registry = productRegistry()
   %PRODUCTREGISTRY Keep RetMIP selectors and DOI provenance in one table.
   registry = [ ...
      "forcing", "10.22008/FK2/GZ3CSN"
      "outputs", "10.22008/FK2/CVPUJL"
      "scripts", "10.22008/FK2/GZ3CSN"];
end
