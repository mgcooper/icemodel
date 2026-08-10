function writeJson(pathname, value)
   %WRITEJSON Write pretty-printed JSON as UTF-8 with one trailing newline.
   %
   %  icemodel.verification.setup.writeJson(pathname, value)
   %
   % Writes the readiness ledgers, preview evidence, and QA JSON. All are
   % write UTF-8 with a trailing newline. writeManifest is an exception: it
   % compares bytes before rewriting, so it keeps its own no-newline write and
   % stays a no-op when nothing changed.

   fid = fopen(pathname, 'w', 'n', 'UTF-8');
   if fid < 0
      error('icemodel:verification:setup:writeJson:openFailed', ...
         'cannot open output file: %s', pathname)
   end
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(value, PrettyPrint=true));
   delete(cleanup)
end
