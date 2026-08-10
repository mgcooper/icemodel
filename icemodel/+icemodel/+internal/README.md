# icemodel.internal

This namespace holds code that manages the toolbox itself: version handling,
release metadata, install support, and citation. It is not model code and not a
public API. Nothing here should be needed to configure or run the model.

Runtime helpers belong one level up, in `icemodel` (see its README). If a
function is used while running a case, it does not belong here even if it feels
like a utility.

One exception: `fullpath` resolves the repository root and is called from
`icemodel.config` and `icemodel.getpath` while running a case. It stays here
because finding the installed toolbox root is a toolbox-management job, and
over 100 call sites depend on the current name.

## A note on the name

MathWorks uses `+internal` as a convention rather than a language feature. The
only special treatment is that tooling hides it: it does not appear in tab
completion, `help` listings, or documentation search. The signal is "not public
API, no compatibility promise."

Their convention scopes it by *ownership*: a `+internal` sits beside the public
code it serves and holds that code's implementation details. This repository
scopes it by *purpose* instead, using it for toolbox management across the whole
project. The alternative would scatter version and install code through the
tree.

Version handling:

- `CITATION.cff` is the persisted software-version source.
- `version` caches that value and provides a process-local override and reset.
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
