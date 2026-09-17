function pathname = canonicalPath(pathname)
   %CANONICALPATH Return the canonical absolute form of a file path.
   %
   %  pathname = icemodel.helpers.canonicalPath(pathname)
   %
   % Resolves symbolic links in the existing components of PATHNAME and
   % removes "." and ".." segments. The target does not have to exist. A
   % relative PATHNAME resolves against the current MATLAB folder (pwd), not
   % against the Java user.dir property (see icemodel.helpers.absolutePath).
   %
   % Returns
   %  pathname : string scalar, the canonical absolute path.
   %
   % See also: icemodel.helpers.absolutePath, icemodel.isPathInside

   % Java canonicalization resolves existing symbolic links and dot segments.
   pathname = icemodel.helpers.absolutePath(pathname);
   pathname = string(java.io.File(char(pathname)).getCanonicalPath());
end
