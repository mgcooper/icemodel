function pathname = fixtureDataRoot(version)
   %FIXTUREDATAROOT Return the data root registered for a release.
   %
   %  pathname = icemodel.verification.setup.fixtureDataRoot()
   %  pathname = icemodel.verification.setup.fixtureDataRoot("v1.1")
   %
   % VERSION defaults to the version in CITATION.cff. Every release capability
   % is staged under test/data.
   %
   % See also: icemodel.test.helpers.formalBaselinePolicy,
   %  icemodel.verification.setup.fixtureFileList

   arguments
      version (1, 1) string = ...
         "v" + string(icemodel.internal.version())
   end

   % Validate the release registration, then return the asset staging root.
   icemodel.test.helpers.formalBaselinePolicy(version);
   pathname = string(icemodel.internal.fullpath('test', 'data'));
end
