function writeManifest(pathname, manifest)
   %WRITEMANIFEST Write one verification family manifest JSON file.
   %
   %  icemodel.verification.setup.writeManifest(pathname, manifest)
   %
   % Inputs
   %  pathname   Destination manifest JSON path.
   %  manifest   Family manifest struct.
   %
   % Role
   %  Setup helper used by import/update tooling. Normal verification workflow
   %  functions read manifests but do not write them.

   % Write pretty JSON so committed manifest diffs stay reviewable.
   fid = fopen(pathname, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, jsonencode(manifest, PrettyPrint=true), 'char');
end
