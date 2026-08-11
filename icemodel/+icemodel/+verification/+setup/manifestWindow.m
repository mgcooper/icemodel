function window = manifestWindow(window_start, window_end)
   %MANIFESTWINDOW Serialize one start/end pair for a JSON manifest.
   %
   %  window = icemodel.verification.setup.manifestWindow( ...
   %     window_start, window_end)
   %
   % This helper builds the two-field record, so importers and RCM staging
   % handle UTC, midnight, and open bounds the same way.

   window = struct( ...
      'start', icemodel.verification.setup.formatManifestTime(window_start), ...
      'end', icemodel.verification.setup.formatManifestTime(window_end));
end
