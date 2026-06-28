function [source_dir, status] = fetchPromice(kwargs)
   %FETCHPROMICE Locate or verify local PROMICE pypromice L3 source caches.
   %
   %  source_dir = icemodel.verification.setup.fetchPromice()
   %  [source_dir, status] = icemodel.verification.setup.fetchPromice( ...
   %     cache_dir="/path", stations="KAN_M", strict=false)
   %
   % Role
   %  Cache validator for local PROMICE station NetCDF files and required AWS
   %  metadata tables. It does not download data; it reports missing local
   %  files so staging can fail early in strict mode or skip per site.

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.stations (1, :) string = strings(1, 0)
      kwargs.products (1, :) string {mustBeMember(kwargs.products, ...
         ["hour", "day", "month"])} = "hour"
      kwargs.strict (1, 1) logical = true
      kwargs.silent (1, 1) logical = false
   end

   cache_dir = kwargs.cache_dir;
   if cache_dir == ""
      cache_dir = defaultCacheDir();
   end
   icemodel.helpers.ensureDirExists(cache_dir);

   % Build status before throwing so callers can inspect the cache manifest in
   % dry-run or skip-missing mode.
   status = cacheStatus(cache_dir, kwargs.products, kwargs.stations);
   missing = string({status(~[status.present]).product});

   if isempty(missing)
      source_dir = cache_dir;
      return
   end

   % Keep the message actionable without embedding download behavior in tests or
   % staging helpers.
   if ~kwargs.silent
      fprintf('\n');
      fprintf('=== PROMICE source cache incomplete ===\n');
      fprintf('Cache directory: %s\n', cache_dir);
      fprintf('Missing products: %s\n', strjoin(missing, ', '));
      fprintf('\nExpected pypromice L3 files under cache_dir/<product>/ and\n');
      fprintf('AWS_sites_metadata.csv, AWS_stations_metadata.csv, AWS_variables.csv\n');
      fprintf('at the cache root. Source portal: https://promice.org\n\n');
   end

   if kwargs.strict
      error('icemodel:verification:fetchPromice:missingSources', ...
         'PROMICE source cache incomplete in %s. Missing: %s', ...
         cache_dir, strjoin(missing, ', '));
   end

   source_dir = cache_dir;
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical PROMICE source-cache directory.
   pathname = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice'));
end

function status = cacheStatus(cache_dir, products, stations)
   %CACHESTATUS Return metadata and station-product validation rows.
   status = repmat(emptyRow(), 1, numel(products) + 1);
   status(1) = metadataRow(cache_dir);
   for k = 1:numel(products)
      status(k + 1) = productRow(cache_dir, products(k), stations);
   end
end

function row = emptyRow()
   %EMPTYROW Prototype status row shared by metadata and product checks.
   row = struct('product', "", 'landing_url', "https://promice.org", ...
      'present', false, 'cache_dir', "", ...
      'missing_files', strings(1, 0), ...
      'missing_stations', strings(1, 0));
end

function row = metadataRow(cache_dir)
   %METADATAROW Validate the AWS metadata tables required by PROMICE readers.
   required = ["AWS_sites_metadata.csv", "AWS_stations_metadata.csv", ...
      "AWS_variables.csv"];
   missing = required(~isfile(fullfile(cache_dir, required)));
   row = emptyRow();
   row.product = "metadata";
   row.cache_dir = cache_dir;
   row.present = isempty(missing);
   row.missing_files = missing;
end

function row = productRow(cache_dir, product, stations)
   %PRODUCTROW Validate one pypromice L3 product directory.
   product_dir = fullfile(cache_dir, product);
   row = emptyRow();
   row.product = product;
   row.cache_dir = string(product_dir);

   if isempty(stations)
      files = dir(fullfile(product_dir, "*_" + product + ".nc"));
      row.present = any(~[files.isdir]);
      return
   end

   missing = strings(1, numel(stations));
   n_missing = 0;
   for station = reshape(stations, 1, [])
      if ~hasStationFile(product_dir, product, station)
         n_missing = n_missing + 1;
         missing(n_missing) = station;
      end
   end
   missing = missing(1:n_missing);
   row.present = isempty(missing);
   row.missing_stations = missing;
end

function tf = hasStationFile(product_dir, product, station)
   %HASSTATIONFILE Match station ids with or without underscore/case variants.
   files = dir(fullfile(product_dir, "*_" + product + ".nc"));
   if isempty(files)
      tf = false;
      return
   end

   names = string({files.name});
   suffix = "_" + product + ".nc";
   found = lower(erase(erase(names, suffix), "_"));
   wanted = lower(erase(string(station), "_"));
   tf = any(found == wanted);
end
