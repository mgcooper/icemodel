function filename = locateNormalizedFile(root, expected, kwargs)
   %LOCATENORMALIZEDFILE Find a file by normalized path/name token.
   %
   % The lookup matches fetch validators that compare case-insensitively and
   % ignore spaces, underscores, and hyphens. The recursive scan accepts both
   % flat manual caches and DOI/package subfolders.

   arguments
      root (1, 1) string
      expected (1, 1) string
      kwargs.error_id (1, 1) string = "icemodel:forcing:fileNotFound"
      kwargs.error_label (1, 1) string = expected
   end

   files = icemodel.forcing.helpers.listSourceFiles(root);
   tokens = icemodel.forcing.helpers.normalizedFileToken(files);
   wanted = icemodel.forcing.helpers.normalizedFileToken(expected);
   idx = find(contains(tokens, wanted), 1, 'first');
   if isempty(idx)
      error(kwargs.error_id, '%s not found: %s under %s', ...
         kwargs.error_label, expected, root)
   end
   filename = files(idx);
end
