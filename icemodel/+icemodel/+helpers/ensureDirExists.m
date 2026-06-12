function pathname = ensureDirExists(pathname)
   %ENSUREDIREXISTS Create one directory when it does not already exist.
   %
   %  pathname = icemodel.helpers.ensureDirExists(pathname)
   %
   % Return the input pathname so callers can compose path creation inline.

   arguments
      pathname (1, :) string
   end

   if ~isfolder(pathname)
      mkdir(pathname);
   end
end
