function [source_dir, status] = fetchKtransect(kwargs)
   %FETCHKTRANSECT Locate or verify the local K-transect PANGAEA source cache.
   %
   %  source_dir = icemodel.verification.setup.fetchKtransect()
   %  [source_dir, status] = icemodel.verification.setup.fetchKtransect( ...
   %     cache_dir="/path", strict=false)
   %
   % Role
   %  Cache validator for the Smeets et al. (2022) PANGAEA.947483 K-transect
   %  annual AWS series and its sensor-height workbook. It does not
   %  auto-download in tests; it reports the DOI landing URL and validates
   %  local files. An empty products selection returns an empty status without
   %  creating or scanning cache_dir.
   %
   % Name-value
   %  cache_dir : string  Local family cache root.
   %  products : string vector  "datasets" and/or "heights" cache products.
   %  strict : logical  Error when a requested product is incomplete.
   %  silent : logical  Suppress retrieval guidance for incomplete products.
   %  create_cache_dir : logical  Create cache_dir before validation.
   %
   % Returns
   %  source_dir : string  Resolved cache root.
   %  status : struct array  One shared fetch-status row per product.

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
      "icemodel:verification:fetchKtransect:unknownProduct", "K-transect");

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
      error_id="icemodel:verification:fetchKtransect:missingSources", ...
      error_label="K-transect", banner_callback=@(cache_dir, status, missing) ...
      icemodel.verification.setup.printFetchProductBanner( ...
      cache_dir, status, missing, "K-transect"));
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical K-transect source-cache directory.
   pathname = icemodel.forcing.helpers.verificationSourceDir("", "ktransect");
end

function status = productStatus(cache_dir, products)
   %PRODUCTSTATUS Return one validation/URL row per requested product.
   status = icemodel.verification.setup.buildFetchProductStatus( ...
      cache_dir, products, productRegistry(), @productRow);
end

function row = productRow(cache_dir, product, doi)
   %PRODUCTROW Build one stable status row.
   switch string(product)
      case "datasets"
         % The annual children extract to datasets/; accept flat and recursive
         % layouts for manually arranged caches.
         patterns = [fullfile("datasets", "K-transect_*.tab"), ...
            "K-transect_*.tab", fullfile("**", "K-transect_*.tab")];
         row = icemodel.verification.setup.fetchProductStatusRow( ...
            cache_dir, product, doi, patterns);
      case "heights"
         % The sensor-height workbook is a series-level PANGAEA attachment.
         patterns = ["metadata_GRL_AWS_*.xlsx", ...
            fullfile("**", "metadata_GRL_AWS_*.xlsx")];
         row = icemodel.verification.setup.fetchProductStatusRow( ...
            cache_dir, product, doi, patterns);
      otherwise
         row = icemodel.verification.setup.fetchProductStatusRow( ...
            cache_dir, product, doi, "*" + product + "*");
   end
end

function registry = productRegistry()
   %PRODUCTREGISTRY Keep K-transect selectors and DOI provenance in one table.
   % Both products publish under the single PANGAEA.947483 series DOI; the
   % per-year child DOIs are pinned from each annual file's own header.
   registry = [ ...
      "datasets", "10.1594/PANGAEA.947483"
      "heights", "10.1594/PANGAEA.947483"];
end
