function tf = isPathInside(pathname, root)
   %ISPATHINSIDE True when a canonical path resolves within a selected root.
   %
   %  tf = icemodel.isPathInside(pathname, root)
   %
   % Role
   %  Canonical-path containment predicate for root-scoped artifact checks
   %  (reconstruction driver, runtime readiness gate, report builder).
   %  Canonicalization through java.io.File resolves symlinks and relative
   %  segments, so a path that escapes the root cannot hide behind
   %  aliasing.
   %
   % Returns
   %  tf : logical scalar, true when pathname is root or inside it.
   %
   % See also: icemodel.forcing.reconstruct.fillPromiceStation,
   %  icemodel.forcing.reconstruct.verifyPromiceFilledReadiness,
   %  icemodel.verification.report.buildGapFillReport

   pathname = string(java.io.File(char(pathname)).getCanonicalPath());
   root = string(java.io.File(char(root)).getCanonicalPath());
   tf = pathname == root || startsWith(pathname, root + filesep);
end
