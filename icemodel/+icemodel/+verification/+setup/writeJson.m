function writeJson(pathname, value)
   %WRITEJSON Write pretty-printed JSON as UTF-8 with one trailing newline.
   %
   %  icemodel.verification.setup.writeJson(pathname, value)
   %
   % One writer for the readiness ledgers, preview evidence, and QA JSON, so
   % they agree on UTF-8 and on the trailing newline instead of each caller
   % choosing. writeManifest is a deliberate exception: it compares bytes
   % before rewriting, so it keeps its own no-newline write to stay a no-op
   % when nothing changed.

   fid = fopen(pathname, 'w', 'n', 'UTF-8');
   if fid < 0
      error('icemodel:verification:setup:writeJson:openFailed', ...
         'cannot open output file: %s', pathname)
   end
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(value, PrettyPrint=true));
   delete(cleanup)
end
