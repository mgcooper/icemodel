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
   %  An empty products selection returns an empty status without creating or
   %  scanning cache_dir.
   %
   % Name-value
   %  cache_dir : string  Local family cache root.
   %  products : string vector  "hourly" and/or "daily" cache products.
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
      "icemodel:verification:fetchImau:unknownProduct", "IMAU");

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
      error_id="icemodel:verification:fetchImau:missingSources", ...
      error_label="IMAU", banner_callback=@(cache_dir, status, missing) ...
      icemodel.verification.setup.printFetchProductBanner( ...
      cache_dir, status, missing, "IMAU"));
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical IMAU source-cache directory.
   pathname = icemodel.forcing.helpers.verificationSourceDir("", "imau");
end

function status = productStatus(cache_dir, products)
   %PRODUCTSTATUS Return one validation/URL row per requested IMAU product.
   status = icemodel.verification.setup.buildFetchProductStatus( ...
      cache_dir, products, productRegistry(), @productRow);
end

function row = productRow(cache_dir, product, doi)
   %PRODUCTROW Build one stable status row.
   switch string(product)
      case "hourly"
         patterns = [fullfile("hourly", "*.tab"), "*.tab", ...
            fullfile("**", "*.tab")];
         row = icemodel.verification.setup.fetchProductStatusRow( ...
            cache_dir, product, doi, patterns, ...
            exclude_folders=filesep + "daily", exclude_names="GRL_");
      case "daily"
         patterns = [fullfile("daily", "GRL_*_AWS.tab"), ...
            "GRL_*_AWS.tab", fullfile("**", "GRL_*_AWS.tab")];
         row = icemodel.verification.setup.fetchProductStatusRow( ...
            cache_dir, product, doi, patterns);
      otherwise
         row = icemodel.verification.setup.fetchProductStatusRow( ...
            cache_dir, product, doi, "*" + product + "*");
   end
end

function registry = productRegistry()
   %PRODUCTREGISTRY Keep IMAU selectors and DOI provenance in one table.
   registry = [ ...
      "hourly", "10.1594/PANGAEA.971647"
      "daily", "10.1594/PANGAEA.970127"];
end
