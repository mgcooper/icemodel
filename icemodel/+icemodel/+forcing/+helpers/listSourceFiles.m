function files = listSourceFiles(root)
   %LISTSOURCEFILES Return flat and nested files under a source/cache root.
   %
   %  files = icemodel.forcing.helpers.listSourceFiles(root)
   %
   % The scan accepts flat manual caches and unpacked DOI/package subfolders, then
   % removes duplicate paths produced by MATLAB recursive globs.

   arguments
      root (1, 1) string
   end

   top = dir(fullfile(root, '*'));
   nested = dir(fullfile(root, '**', '*'));
   all_files = [top(:); nested(:)];
   all_files = all_files(~[all_files.isdir]);
   if isempty(all_files)
      files = strings(1, 0);
   else
      files = unique(string(fullfile({all_files.folder}, {all_files.name})), ...
         'stable');
   end
end
