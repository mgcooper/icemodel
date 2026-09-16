function makecontents(option, opts)
   %MAKECONTENTS Write a Contents.m file in each namespace folder.
   %
   %  icemodel.internal.makecontents()
   %  icemodel.internal.makecontents("-backup")
   %  icemodel.internal.makecontents(_, Folder=PATH)
   %
   %  Run this after you add, rename, or remove a file in a namespace folder,
   %  and after you edit an H1 line.
   %
   %  makecontents() writes a Contents.m file in each namespace folder (a
   %  folder whose name starts with "+") under the icemodel code folder,
   %  icemodel.internal.fullpath("icemodel"). Each listing also covers the
   %  namespace folders below it, so a namespace that holds only namespaces,
   %  such as icemodel.test, gets a listing too. `help icemodel.column` then
   %  prints the listing of that namespace.
   %
   %  makecontents("-backup") moves each earlier Contents.m to a dated file
   %  in tempdir and prints its path. Backups are off by default because git
   %  tracks every Contents.m. Note: a backup beside the original would be a
   %  callable member of that namespace.
   %
   %  makecontents(_, Folder=PATH) works on the folder PATH instead, so a
   %  test can run it on a scratch copy.
   %
   % See also: updatecontents

   arguments
      option (1, 1) string ...
         {mustBeMember(option, ["-backup", "-nobackup"])} = "-nobackup"
      opts.Folder (1, 1) string = icemodel.internal.fullpath("icemodel")
   end

   % A missing folder is a caller error, so stop before any file moves.
   if ~isfolder(opts.Folder)
      error('icemodel:internal:makecontents:folderNotFound', ...
         'folder not found: %s', opts.Folder)
   end

   % Every namespace folder at any depth gets its own listing.
   listing = dir(fullfile(opts.Folder, '**', '+*'));
   listing = listing([listing.isdir]);
   folders = sort(string(fullfile({listing.folder}, {listing.name})));
   for folder = reshape(folders, 1, [])
      writeOneContents(folder, option == "-backup");
   end
end

function writeOneContents(folder, dobackup)
   %WRITEONECONTENTS Write the Contents.m file of one namespace folder.
   %
   % The earlier Contents.m moves to a temporary file first, so the new
   % listing does not name it and a failed write can restore it.

   contentsfile = fullfile(folder, 'Contents.m');
   previous = "";
   if isfile(contentsfile)
      previous = string(tempname);
      movefile(contentsfile, previous);
   end

   % Restore the earlier listing if the write fails, then report the error.
   try
      updatecontents(folder);
   catch err
      if previous ~= ""
         movefile(previous, contentsfile);
      end
      rethrow(err)
   end

   % Keep or discard the earlier listing after a successful write.
   if previous == ""
      return
   end
   if dobackup
      [~, name] = fileparts(folder);
      backup = fullfile(tempdir, extractAfter(name, 1) + "_Contents_" ...
         + string(datetime('now', 'Format', 'yyyyMMdd''T''HHmmssSSS')) + ".m");
      movefile(previous, backup);
      fprintf('Backed up %s to %s\n', contentsfile, backup);
   else
      delete(previous)
   end
end
