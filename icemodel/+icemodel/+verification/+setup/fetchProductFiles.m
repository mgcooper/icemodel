function files = fetchProductFiles(cache_dir, patterns, kwargs)
   %FETCHPRODUCTFILES Return files matching product-cache patterns.
   %
   %  files = icemodel.verification.setup.fetchProductFiles(cache_dir, patterns)
   %
   % Patterns are relative to cache_dir and can include recursive globs. The
   % optional exclusions leave the product-specific filename rules with the
   % callers, while this function collects and de-duplicates the files.

   arguments
      cache_dir (1, 1) string
      patterns (1, :) string
      kwargs.exclude_folders (1, :) string = strings(1, 0)
      kwargs.exclude_names (1, :) string = strings(1, 0)
   end

   % Each pattern contributes one match block. Collect the blocks in a buffer
   % sized to the pattern list, then concatenate once so the scan does not
   % reallocate the hit list per pattern.
   hit_blocks = repmat({strings(1, 0)}, 1, numel(patterns));
   for k = 1:numel(patterns)
      listing = dir(fullfile(cache_dir, patterns(k)));
      listing = listing(~[listing.isdir]);
      if isempty(listing)
         continue
      end
      hit_blocks{k} = string(fullfile({listing.folder}, {listing.name}));
   end
   hits = unique([strings(1, 0), hit_blocks{:}], 'stable');
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
