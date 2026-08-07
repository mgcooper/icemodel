function assertArtifactSha256(pathname, expected_sha256)
   %ASSERTARTIFACTSHA256 Require current bytes to match a pinned SHA-256.
   %
   %  icemodel.verification.helpers.assertArtifactSha256( ...
   %     pathname, expected_sha256)

   arguments
      pathname (1, 1) string
      expected_sha256 (1, 1) string
   end

   % A readiness identity is usable only while both the artifact and its
   % complete hash remain available at the execution boundary.
   if ~isfile(pathname)
      error('icemodel:verification:artifactIdentity:missing', ...
         'Pinned artifact is unavailable: %s', pathname)
   end
   if isempty(regexp(expected_sha256, '^[0-9a-fA-F]{64}$', 'once'))
      error('icemodel:verification:artifactIdentity:sha256', ...
         'Pinned artifact has no valid SHA-256 identity: %s', pathname)
   end

   % Rehash the current bytes so a readiness ledger cannot silently outlive
   % the artifact version that it admitted.
   actual_sha256 = ...
      icemodel.verification.setup.fileSha256(pathname);
   if ~strcmpi(actual_sha256, expected_sha256)
      error('icemodel:verification:artifactIdentity:mismatch', ...
         'Pinned artifact bytes changed after readiness: %s', pathname)
   end
end
