function relpath = fixtureRelativePosix(root, pathname)
   %FIXTURERELATIVEPOSIX Return one root-relative path with POSIX separators.
   %
   %  relpath = icemodel.verification.setup.fixtureRelativePosix(root, pathname)
   %
   % Both slash forms are normalized explicitly so Windows-style paths can be
   % validated on every host rather than only when filesep is a backslash.

   arguments
      root (1, 1) string
      pathname (1, 1) string
   end

   % Remove trailing root separators before slicing the known root prefix.
   backslash = string(char(92));
   separators = ["/", backslash];
   while strlength(root) > 0 && any(endsWith(root, separators))
      root = extractBefore(root, strlength(root));
   end
   relpath = extractAfter(pathname, strlength(root));

   % Strip the path-boundary separator, then normalize every remaining level.
   while strlength(relpath) > 0 && any(startsWith(relpath, separators))
      relpath = extractAfter(relpath, 1);
   end
   relpath = replace(relpath, backslash, "/");
   relpath = replace(relpath, string(filesep), "/");
end
