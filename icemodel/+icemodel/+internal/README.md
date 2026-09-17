# icemodel.internal

This namespace manages toolbox versions, releases, installation, and citation.

Runtime helpers belong one level up, in `icemodel` (see its README).

One exception: `fullpath` resolves the repository root and is called from
`icemodel.config` and `icemodel.getpath` while running a case. It stays here
because finding the installed toolbox root is a toolbox-management job.

Version handling:

- `CITATION.cff` is the software-version source.
- `version` caches that value and provides a process-local override and reset.
- `readCffVersion` keeps the version-loading path compatible with the
  model's minimum MATLAB version compatibility.

Release maintenance:

- `releaseMetadata("prepare", ...)` checks release state and validates staged
  CFF metadata.
- `releaseMetadata("observe", ...)` performs bounded public GitHub and Zenodo
  checks.
- `releaseMetadata("finalize", ...)` verifies Zenodo DOI lineage before
  updating and validating the CFF identifiers.

Namespace listings:

- `makecontents` writes a generated `Contents.m` in every namespace folder
  under `icemodel/`, so `help icemodel.column` prints the functions of that
  namespace with their H1 lines. Run `icemodel.internal.makecontents()` after
  you add, rename, or remove a file in a namespace folder, and after you edit
  an H1 line. Do not edit a `Contents.m` by hand.
