function sha256 = policySha256()
   %POLICYSHA256 Return the SHA-256 fingerprint of reconstruction POLICY.md.
   %
   %  sha256 = icemodel.forcing.reconstruct.policySha256()
   %
   % Producers and consumers call this single source so a staged
   % promice_filled artifact is accepted only under the policy text that the
   % running code ships.

   % Resolve POLICY.md beside this helper so callers cannot substitute a
   % workspace-relative or current-directory-dependent policy file.
   policy_file = fullfile(fileparts(mfilename('fullpath')), 'POLICY.md');
   sha256 = icemodel.verification.setup.fileSha256(policy_file);
end
