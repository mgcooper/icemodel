function sha256 = policySha256()
   %POLICYSHA256 Return the SHA-256 fingerprint of reconstruction POLICY.md.
   %
   %  sha256 = icemodel.forcing.reconstruct.policySha256()
   %
   % The producer stamps this digest into every promice_filled artifact.
   % A consumer compares against it by default, so a staged artifact is valid
   % only under the policy text the running code ships. A frozen release
   % overrides that comparison with the digest it registered:
   % formalBaselinePolicy pins promice_filled_policy_sha256, and
   % assertPromiceFilledArtifact prefers that pin. Editing POLICY.md therefore
   % changes the rolling comparison, not a registered release comparison.

   % Resolve POLICY.md beside this helper so callers cannot substitute a
   % workspace-relative or current-directory-dependent policy file.
   policy_file = fullfile(fileparts(mfilename('fullpath')), 'POLICY.md');
   sha256 = icemodel.verification.setup.fileSha256(policy_file);
end
