function files = fetchProductFiles(cache_dir, patterns, kwargs)
   %FETCHPRODUCTFILES Return files matching product-cache patterns.
   %
   %  files = icemodel.verification.setup.fetchProductFiles(cache_dir, patterns)
   %
   % Patterns are relative to cache_dir and may include recursive globs. Optional
   % exclusions keep product-specific filename rules in callers while sharing the
   % filesystem collection and de-duplication behavior.

   arguments
      cache_dir (1, 1) string
      patterns (1, :) string
      kwargs.exclude_folders (1, :) string = strings(1, 0)
      kwargs.exclude_names (1, :) string = strings(1, 0)
   end

   hits = strings(1, 0);
   for pattern = patterns
      listing = dir(fullfile(cache_dir, pattern));
      listing = listing(~[listing.isdir]);
      if isempty(listing)
         continue
      end
      paths = string(fullfile({listing.folder}, {listing.name}));
      hits = [hits, paths]; %#ok<AGROW>
   end
   hits = unique(hits, 'stable');
   if isempty(hits)
      files = hits;
      return
   end

   [folders, names, extensions] = fileparts(hits);
   basenames = names + extensions;
   keep = true(size(hits));
   for folder = kwargs.exclude_folders
      keep = keep & ~contains(folders, folder);
   end
   for name = kwargs.exclude_names
      keep = keep & ~startsWith(basenames, name);
   end
   files = hits(keep);
end
