function pathname = assertRootRelativeArtifactSha256( ...
      root, relative_path, expected_sha256)
   %ASSERTROOTRELATIVEARTIFACTSHA256 Verify one root-scoped artifact identity.
   %
   %  pathname = ...
   %     icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
   %     root, relative_path, expected_sha256)

   arguments
      root (1, 1) string
      relative_path (1, 1) string
      expected_sha256 (1, 1) string
   end

   % Producer manifests use selected-data-root-relative POSIX paths. Reject
   % absolute and parent-relative spellings before fullfile can reinterpret them.
   relative_path = replace(relative_path, "\", "/");
   segments = split(relative_path, "/");
   unsafe = relative_path == "" || startsWith(relative_path, "/") ...
      || ~isempty(regexp(relative_path, '^[A-Za-z]:', 'once')) ...
      || any(segments == "..");
   pathname = string(fullfile(root, relative_path));
   if unsafe || ~icemodel.isPathInside(pathname, root)
      error('icemodel:verification:artifactIdentity:relativePath', ...
         'Pinned artifact path must remain relative to its selected root: %s', ...
         relative_path)
   end

   % Reuse the shared hash check once the path is known to be ours.
   icemodel.verification.helpers.assertArtifactSha256( ...
      pathname, expected_sha256);
end
