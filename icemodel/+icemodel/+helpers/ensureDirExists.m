function pathname = ensureDirExists(pathname)
   %ENSUREDIREXISTS Create one directory when it does not already exist.
   %
   %  pathname = icemodel.helpers.ensureDirExists(pathname)
   %
   % Return the input pathname so callers can compose path creation inline.

   arguments
      pathname (1, :) string
   end

   if exist(pathname, 'dir') ~= 7
      mkdir(pathname);
   end
end
