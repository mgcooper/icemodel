function tf = isPathInside(pathname, root)
   %ISPATHINSIDE True when a canonical path resolves within a selected root.
   %
   %  tf = icemodel.internal.isPathInside(pathname, root)
   %
   % Role
   %  Single canonical-path containment predicate for every root-scoped
   %  artifact check (reconstruction driver, runtime readiness gate,
   %  report builder). Canonicalization through java.io.File resolves
   %  symlinks and relative segments so escapes cannot hide behind
   %  aliasing.
   %
   % Returns
   %  tf : logical scalar, true when pathname is root or inside it.
   %
   % See also: icemodel.forcing.reconstruct.fillPromiceStation,
   %  icemodel.verifyPromiceFilledReadiness,
   %  icemodel.verification.report.buildGapFillReport

   pathname = string(java.io.File(char(pathname)).getCanonicalPath());
   root = string(java.io.File(char(root)).getCanonicalPath());
   tf = pathname == root || startsWith(pathname, root + filesep);
end
