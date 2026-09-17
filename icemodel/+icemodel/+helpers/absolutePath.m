function pathname = absolutePath(pathname)
   %ABSOLUTEPATH Return a path anchored at the current MATLAB folder.
   %
   %  pathname = icemodel.helpers.absolutePath(pathname)
   %
   % An absolute PATHNAME returns unchanged. A relative PATHNAME returns
   % fullfile(pwd, PATHNAME). The function does not resolve symbolic links or
   % dot segments.
   %
   % Call it before a Java path API sees a relative path. MATLAB cd does not
   % keep the Java user.dir property current, so java.io.File and
   % java.nio.file.Paths resolve a relative path against an earlier folder.
   %
   % Returns
   %  pathname : string scalar.
   %
   % See also: icemodel.helpers.canonicalPath

   % java.io.File classifies absolute paths the same way on every platform.
   pathname = string(pathname);
   if ~java.io.File(char(pathname)).isAbsolute()
      pathname = string(fullfile(pwd, pathname));
   end
end
