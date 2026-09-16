function tf = isPathInside(pathname, root)
   %ISPATHINSIDE True when a canonical path resolves within a selected root.
   %
   %  tf = icemodel.isPathInside(pathname, root)
   %
   % Role
   %  Checks whether a canonical path is inside a selected root. Root-scoped
   %  artifact checks use it (reconstruction driver, runtime readiness gate,
   %  report builder). Canonicalization through icemodel.helpers.canonicalPath
   %  resolves symlinks and relative segments, so a symlink or a ".." segment
   %  cannot move a path outside the root without detection. A relative path
   %  resolves against the current MATLAB folder.
   %
   % Returns
   %  tf : logical scalar, true when pathname is root or inside it.
   %
   % See also: icemodel.forcing.reconstruct.fillPromiceStation,
   %  icemodel.forcing.reconstruct.verifyPromiceFilledReadiness,
   %  icemodel.verification.report.buildGapFillReport

   % Compare resolved paths so aliases and dot segments cannot hide an escape.
   pathname = icemodel.helpers.canonicalPath(pathname);
   root = icemodel.helpers.canonicalPath(root);
   tf = pathname == root || startsWith(pathname, root + filesep);
end
