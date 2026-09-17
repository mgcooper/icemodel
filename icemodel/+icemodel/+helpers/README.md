# icemodel.helpers

Purpose: shared runtime utilities used by two or more namespaces.

Contents:

- `copyFields`
  - Copies one scalar struct's fields onto another. By default, it adds all
    source fields. Set `only_matching = true` to copy only fields already in
    the target. `icemodel.couplers.update_solver_diag` uses this option to
    promote the accepted `diag.substep` record.
- `ensureDirExists`
  - Creates a directory when needed.
- `absolutePath`
  - Anchors a relative path at the current MATLAB folder. Call it before a
    Java path API sees a relative path, because MATLAB `cd` does not keep the
    Java `user.dir` property current.
- `canonicalPath`
  - Returns the absolute path with symbolic links in existing components and
    dot segments resolved. The target does not have to exist.
  - These functions call it:
    - `icemodel.isPathInside`
    - `icemodel.internal.releaseMetadata`
    - `icemodel.forcing.reconstruct.verifyPromiceFilledReadiness`
    - `icemodel.verification.setup.packFixtures`
    - `icemodel.verification.setup.fetchFixtures`
    - `icemodel.verification.setup.fixtureCallerSymlink`
    - `icemodel.verification.setup.repairRcmArtifactMetadata`
    - `icemodel.test.helpers.resolveReleaseDataRoots`
    - `icemodel.test.helpers.worktreeRevision`
- `rmttleapinds`
  - Removes leap-day rows from a timetable.
