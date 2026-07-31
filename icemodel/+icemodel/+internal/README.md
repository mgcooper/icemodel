# icemodel.internal

This namespace contains repository support code, not public IceModel APIs.

Version handling:

- `CITATION.cff` is the persisted software-version source.
- `version` caches that value and retains its process-local override and reset
  behavior.
- `readCffVersion` keeps the version-loading path compatible with the core
  model's documented MATLAB floor.

Release maintenance:

- `releaseMetadata("prepare", ...)` checks release state and validates staged
  CFF metadata.
- `releaseMetadata("observe", ...)` performs bounded public GitHub and Zenodo
  checks.
- `releaseMetadata("finalize", ...)` verifies Zenodo DOI lineage before
  updating and validating the CFF identifiers.

`releaseMetadata` is maintainer tooling and may use modern MATLAB features. It
does not merge, tag, push, publish a release, create a pull request, or invoke
the shared release skill. The maintainer release section in the repository
README gives the complete sequence.
