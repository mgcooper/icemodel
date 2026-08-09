function exportAndClose(fig, filename)
   %EXPORTANDCLOSE Export one report figure and release graphics state.

   cleanup = onCleanup(@() close(fig));
   exportgraphics(fig, filename, Resolution=160)
   clear cleanup
end
