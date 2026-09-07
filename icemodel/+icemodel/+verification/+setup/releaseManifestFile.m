function pathname = releaseManifestFile(version, kwargs)
   %RELEASEMANIFESTFILE Return the release-data manifest path for a version.
   %
   %  pathname = icemodel.verification.setup.releaseManifestFile()
   %  pathname = icemodel.verification.setup.releaseManifestFile("v1.1")
   %  filename = icemodel.verification.setup.releaseManifestFile("v1.1", ...
   %     directory="")
   %
   % VERSION defaults to the version in CITATION.cff and carries the leading
   % "v".
   %
   % Name-value
   %  directory  Folder that holds the manifest. The default is the tracked
   %             test/assets folder. Pass a staging folder to name a manifest
   %             that packFixtures writes. A blank directory returns the bare
   %             filename, which resolves against the current folder.
   %
   % See also: icemodel.verification.setup.fixtureDataRoot,
   %  icemodel.verification.setup.fixtureFileList,
   %  icemodel.verification.setup.packFixtures

   arguments
      version (1, 1) string = ...
         "v" + string(icemodel.internal.version())
      kwargs.directory (1, 1) string = ...
         string(icemodel.internal.fullpath('test', 'assets'))
   end

   % One template owns the manifest filename, so every producer and consumer
   % resolves the same name for a release.
   filename = "icemodel-" + version + "-data-manifest.json";

   % fullfile drops a blank part, so a blank directory returns the bare
   % filename without a separate branch.
   pathname = string(fullfile(kwargs.directory, filename));
end
