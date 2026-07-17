function window = manifestWindow(window_start, window_end)
   %MANIFESTWINDOW Serialize one start/end pair for a JSON manifest.
   %
   %  window = icemodel.verification.setup.manifestWindow( ...
   %     window_start, window_end)
   %
   % Keeping this two-field record in one helper prevents importers and RCM
   % staging from drifting in their handling of UTC, midnight, and open bounds.

   window = struct( ...
      'start', icemodel.verification.setup.formatManifestTime(window_start), ...
      'end', icemodel.verification.setup.formatManifestTime(window_end));
end
