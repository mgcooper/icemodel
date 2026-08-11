function writeJson(pathname, value)
   %WRITEJSON Write pretty-printed JSON as UTF-8 with one trailing newline.
   %
   %  icemodel.verification.setup.writeJson(pathname, value)
   %
   % This function writes the readiness ledgers, the preview evidence, and the
   % QA JSON. Each file is UTF-8 with a trailing newline. writeManifest is an
   % exception: it compares bytes before it rewrites, so it keeps its own write
   % with no trailing newline and changes nothing when the bytes match.

   fid = fopen(pathname, 'w', 'n', 'UTF-8');
   if fid < 0
      error('icemodel:verification:setup:writeJson:openFailed', ...
         'cannot open output file: %s', pathname)
   end
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(value, PrettyPrint=true));
   delete(cleanup)
end
