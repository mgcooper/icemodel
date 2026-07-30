function [source_dir, status] = fetchGcnet(kwargs)
   %FETCHGCNET Locate or verify local Vandecrux/GC-Net source caches.
   %
   %  source_dir = icemodel.verification.setup.fetchGcnet()
   %  [source_dir, status] = icemodel.verification.setup.fetchGcnet( ...
   %     cache_dir="/path", stations=["DYE_2", "Summit"], strict=false)
   %
   % Role
   %  Cache validator for the Vandecrux/GC-Net Arctic Data packages used by
   %  RetMIP Dye-2-long and Summit staging. It accepts the user's current flat
   %  manual cache and DOI package subfolders, but never downloads files in
   %  tests or importer setup.
   %  An empty products selection returns an empty status without creating or
   %  scanning cache_dir.
   %
   % Name-value
   %  cache_dir : string  Local family cache root.
   %  stations : string vector  Physical GC-Net station subset.
   %  products : string vector  Vandecrux/GC-Net product subset.
   %  strict : logical  Error when a requested product is incomplete.
   %  silent : logical  Suppress retrieval guidance for incomplete products.
   %  create_cache_dir : logical  Create cache_dir before validation.
   %
   % Returns
   %  source_dir : string  Resolved cache root.
   %  status : struct array  One station-aware row per requested product.

   arguments
      kwargs.cache_dir (1, 1) string = defaultCacheDir()
      kwargs.stations (1, :) string = ["DYE_2", "Summit"]
      kwargs.products (1, :) string ...
         {icemodel.verification.validators.mustBeGcnetProductSelection( ...
         kwargs.products)} = ...
         icemodel.verification.setup.gcnetProductNames()
      kwargs.strict (1, 1) logical = true
      kwargs.silent (1, 1) logical = false
      kwargs.create_cache_dir (1, 1) logical = true
   end

   cache_dir = icemodel.verification.setup.resolveFetchCacheDir( ...
      kwargs.cache_dir, defaultCacheDir());
   if kwargs.create_cache_dir && ~isempty(kwargs.products)
      icemodel.helpers.ensureDirExists(cache_dir);
   end

   % Resolve dataset-specific request aliases before status construction.
   stations = canonicalStations(kwargs.stations);

   % Build cache status before strict handling.
   status = productStatus(cache_dir, kwargs.products, stations);
   source_dir = icemodel.verification.setup.finishFetchStatus( ...
      cache_dir, status, strict=kwargs.strict, silent=kwargs.silent, ...
      error_id="icemodel:verification:fetchGcnet:missingSources", ...
      error_label="GC-Net/Vandecrux", ...
      banner_callback=@printRetrievalBanner);
end

function printRetrievalBanner(cache_dir, status, missing)
   %PRINTRETRIEVALBANNER Print GC-Net/Vandecrux retrieval instructions.
   fprintf('\n');
   fprintf('=== GC-Net/Vandecrux source cache incomplete ===\n');
   fprintf('Cache directory: %s\n', cache_dir);
   fprintf('Missing products: %s\n', strjoin(missing, ', '));
   fprintf('\nRetrieval:\n');
   for k = 1:numel(status)
      fprintf('  %-16s DOI: %s\n', status(k).product, status(k).doi);
      fprintf('                   URL: %s\n', status(k).landing_url);
   end
   fprintf('\nMissing files/patterns:\n');
   for k = find(~[status.present])
      for j = 1:numel(status(k).missing_files)
         fprintf('  - %s\n', status(k).missing_files(j));
      end
   end
   % Preserve the established banner for ordinary missing files; add an
   % ambiguity section only when multiple candidates require user cleanup.
   has_ambiguity = arrayfun(@(row) ~isempty(row.ambiguous_files), status);
   if any(has_ambiguity)
      fprintf('\nAmbiguous file matches:\n');
      for k = find(has_ambiguity)
         for j = 1:numel(status(k).ambiguous_files)
            fprintf('  - %s\n', status(k).ambiguous_files(j));
         end
      end
   end
   fprintf('\nPlace Arctic Data downloads either directly under cache_dir\n');
   fprintf('or in DOI/package subfolders. No automatic download is attempted.\n\n');
end

function pathname = defaultCacheDir()
   %DEFAULTCACHEDIR Canonical Vandecrux/GC-Net source-cache directory.
   pathname = icemodel.forcing.helpers.verificationSourceDir("", "gcnet");
end

function stations = canonicalStations(stations)
   %CANONICALSTATIONS Map RetMIP aliases to Vandecrux file prefixes.
   if isempty(stations)
      stations = ["DYE_2", "Summit"];
   end
   stations = icemodel.forcing.helpers.gcnetVandecruxStation(stations);
   stations = unique(stations, 'stable');
end

function status = productStatus(cache_dir, products, stations)
   %PRODUCTSTATUS Return one validation/URL row per requested product.
   status = repmat(emptyRow(), 1, numel(products));
   for k = 1:numel(products)
      status(k) = productRow(products(k), cache_dir, stations);
   end
end

function row = emptyRow()
   %EMPTYROW Extend the literal shared fetch row with GC-Net-only diagnostics.
   row = icemodel.verification.setup.emptyFetchProductStatusRow();
   row.stations = strings(1, 0);
   row.resolved_files = repmat(resolvedEntry("", "", ""), 1, 0);
   row.missing_files = strings(1, 0);
   row.ambiguous_files = strings(1, 0);
   row.provenance_file = "";
   row.missing_provenance = strings(1, 0);
end

function row = productRow(product, cache_dir, stations)
   %PRODUCTROW Build one stable row from uniquely resolved station files.
   spec = icemodel.verification.setup.gcnetProductSpec(product);
   expected = expectedStationEntries(spec, stations);
   row = icemodel.verification.setup.fetchProductStatusRow( ...
      cache_dir, product, spec.doi, spec.data_patterns, cache_subdir="");

   % Require one, and only one, basename-normalized match for every registry
   % rule. This preserves accepted separator/case aliases without allowing a
   % partial basename or duplicate package copy to claim product completeness.
   [found, resolved, missing, ambiguous] = stationFileStatus( ...
      row.found_files, expected);
   row.present = isempty(missing) && isempty(ambiguous);
   row.found_files = found;
   row.stations = stations;
   row.resolved_files = resolved;
   row.missing_files = missing;
   row.ambiguous_files = ambiguous;
   [row.provenance_file, row.missing_provenance] = provenanceStatus( ...
      cache_dir, spec);
end

function [found, resolved, missing, ambiguous] = ...
      stationFileStatus(files, expected)
   %STATIONFILESTATUS Resolve every registry rule without arbitrary selection.
   [~, names, extensions] = fileparts(files);
   tokens = normalizeFileToken(names + extensions);
   found = strings(1, 0);
   resolved = repmat(resolvedEntry("", "", ""), 1, 0);
   missing = strings(1, 0);
   ambiguous = strings(1, 0);

   % A zero-match rule is missing; a multi-match rule is both missing from the
   % usable product and explicitly ambiguous for retrieval diagnostics.
   for k = 1:numel(expected)
      matches = expectedFileMatches(tokens, expected(k));
      if nnz(matches) == 1
         filename = files(matches);
         found(end + 1) = filename; %#ok<AGROW>
         resolved(end + 1) = resolvedEntry( ...
            expected(k).station, expected(k).suffix, filename); %#ok<AGROW>
      else
         missing(end + 1) = expected(k).display; %#ok<AGROW>
         if nnz(matches) > 1
            ambiguous = [ambiguous, files(matches)]; %#ok<AGROW>
         end
      end
   end
   found = unique(found, 'stable');
   ambiguous = unique(ambiguous, 'stable');
end

function matches = expectedFileMatches(tokens, expected)
   %EXPECTEDFILEMATCHES Apply the exact or numbered registry filename rule.
   if expected.mode == "exact"
      matches = tokens == expected.token;
      return
   end

   % Simulated-firn suffixes end immediately before a numeric bin id. Anchoring
   % both ends rejects partial downloads, backup suffixes, and other products.
   expression = ['^' regexptranslate('escape', char(expected.token)) ...
      '[0-9]+\.nc$'];
   matches = ~cellfun('isempty', regexp(cellstr(tokens), expression, 'once'));
   matches = reshape(matches, size(tokens));
end

function entries = expectedStationEntries(spec, stations)
   %EXPECTEDSTATIONENTRIES Build required station filename checks for one product.
   suffixes = reshape(string(spec.station_suffixes), 1, []);
   n_entries = numel(stations) * numel(suffixes);
   entries = repmat(expectedEntry("", "", "", "", ""), 1, n_entries);

   % Station files prove the RetMIP-relevant DYE_2/Summit products are usable.
   % Package XML is tracked separately as optional provenance.
   n = 0;
   for station = reshape(stations, 1, [])
      for suffix = suffixes
         n = n + 1;
         display = station + suffix;
         if spec.station_file_mode == "numbered"
            display = display + "<index>.nc";
         end
         token = station + suffix;
         entries(n) = expectedEntry( ...
            station, suffix, display, token, spec.station_file_mode);
      end
   end
end

function [filename, missing] = provenanceStatus(cache_dir, spec)
   %PROVENANCESTATUS Report optional Arctic Data XML provenance availability.
   metadata_file = string(spec.metadata_file);
   files = icemodel.verification.setup.fetchProductFiles( ...
      cache_dir, [metadata_file, fullfile("**", metadata_file)]);
   if isempty(files)
      filename = "";
      missing = metadata_file;
   else
      filename = files(1);
      missing = strings(1, 0);
   end
end

function entry = expectedEntry(station, suffix, display, token, mode)
   %EXPECTEDENTRY Store one display string and one registry-backed match rule.
   entry = struct( ...
      'station', string(station), ...
      'suffix', string(suffix), ...
      'display', string(display), ...
      'token', normalizeFileToken(token), ...
      'mode', string(mode));
end

function entry = resolvedEntry(station, suffix, filename)
   %RESOLVEDENTRY Bind one registry identity to its unique classified file.
   entry = struct('station', string(station), 'suffix', string(suffix), ...
      'filename', string(filename));
end

function token = normalizeFileToken(value)
   %NORMALIZEFILETOKEN Compare filenames case-insensitively across separators.
   token = icemodel.forcing.helpers.normalizedFileToken(value);
end
