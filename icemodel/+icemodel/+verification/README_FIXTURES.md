# Test/Demo Fixture Data: Release-With-Assets Strategy

The committed binary fixture DATA used by the unit and verification suites is
heavy: ~119 MB of `.mat` files under `demo/data/`. To keep the repo lean this
data can be moved out of version control and into a **versioned GitHub release
asset** that is fetched / verified on demand, while the repo keeps only the lean
**manifests** and a **fetcher**.

This document describes the tooling, the workflow, and the one-time migration.

> Status: the **tooling and migration path** are delivered. The actual flip
> (deleting the committed fixtures and uploading the asset) is a deliberate
> user step because it requires creating a GitHub release. It is tracked as
> **`icemodel-1ps.17`**. Until the flip, the committed fixtures remain in the
> repo and offline CI is unchanged.

## What moves vs what stays

| Path | After the flip |
|------|----------------|
| `demo/data/eval/**/*.mat` (obs/reference bundles, ~14 MB) | moves to the asset |
| `demo/data/input/met/*.mat` (forcing, ~91 MB) | moves to the asset |
| `demo/data/input/userdata/*.mat` (per-year userdata, ~15 MB) | moves to the asset |
| `demo/data/eval/<family>/manifest.json` | **stays committed** (lean, reviewable) |
| `.gitkeep` placeholders | **stays committed** (preserve tree structure) |
| `data/verification/**` raw source caches | already gitignored (unchanged) |

The fixture file set is enumerated in exactly one place,
`icemodel.verification.setup.fixtureFileList`, so the packed set and the verified
set can never drift.

## Tooling

### `packFixtures` — producer

```matlab
result = icemodel.verification.setup.packFixtures("v0.1.0");
```

Bundles the fixture data into `release-staging/icemodel-fixtures-v0.1.0.tar.gz`
plus a sidecar `release-staging/icemodel-fixtures-v0.1.0.MANIFEST.json`. The
manifest records, for every fixture file, its repo-relative path, SHA-256, and
byte size, plus the version and a UTC creation timestamp, so the bundle is
verifiable. `release-staging/` is gitignored; the bundle is never committed.

`result` reports the committed footprint, the compressed asset size, and the
saving. Archives compress the `.mat` payload substantially, so the asset is much
smaller than the committed footprint.

### `fetchFixtures` — consumer / verifier (mirrors `fetchSumup`)

```matlab
% Verify the on-disk fixtures against a staged/asset manifest:
icemodel.verification.setup.fetchFixtures("v0.1.0");

% Provision from a downloaded archive (extract + verify):
icemodel.verification.setup.fetchFixtures("v0.1.0", ...
   archive="icemodel-fixtures-v0.1.0.tar.gz", extract=true);
```

`fetchFixtures` follows the same verify-cache contract as `fetchSumup` /
`fetchEsmSnowmip`: it verifies a local cache and, when fixtures are missing,
prints actionable instructions to download the release asset — it **never
auto-downloads**. Modes, in priority order:

1. **archive given** → extract into `demo/data`, then SHA-256-verify every file.
2. **manifest resolvable** (explicit, archive sibling, or version-named staged) →
   SHA-256-verify the on-disk fixtures.
3. **no manifest** → presence-only check. This is the non-breaking path used
   today: the committed fixtures are present and no bundle manifest exists, so
   the call is a clean no-op.

### Bootstrap wiring (non-breaking)

`icemodel.test.helpers.bootstrapTestEnvironment` calls `fetchFixtures` in lenient
mode (`strict=false, silent=true`). With the committed fixtures present it is a
no-op; after the flip it is the on-demand provisioning hook. It never aborts
bootstrap on its own — individual tests still fail clearly if a fixture they
genuinely need is absent.

## One-time migration (the flip — `icemodel-1ps.17`)

Run these once, when you are ready to move the fixtures to a release asset:

```bash
# 1. Pack the current committed fixtures into a versioned bundle.
/Applications/MATLAB_R2024b.app/bin/matlab -batch \
  "icemodel.verification.setup.packFixtures('v0.1.0')"

# 2. Create the GitHub release and upload the archive + manifest as assets.
gh release create fixtures-v0.1.0 \
  release-staging/icemodel-fixtures-v0.1.0.tar.gz \
  release-staging/icemodel-fixtures-v0.1.0.MANIFEST.json \
  --repo mgcooper/icemodel \
  --title "Test/demo fixtures v0.1.0" \
  --notes "Demo/test fixture data (.mat) for the icemodel unit + verification suites. Provision with icemodel.verification.setup.fetchFixtures('v0.1.0', archive=..., extract=true)."

# 3. Remove the heavy fixture .mat from version control (keep the working-tree
#    copies; --cached deletes only the tracked entry). Manifests/.gitkeep stay.
git rm --cached $(git ls-files 'demo/data/eval/**/*.mat')
git rm --cached $(git ls-files 'demo/data/input/met/*.mat')
git rm --cached $(git ls-files 'demo/data/input/userdata/*.mat')

# 4. Activate the .gitignore flip block (uncomment the "release-asset flip"
#    block in .gitignore so the fixture .mat are ignored while manifest.json /
#    .gitkeep stay committed).

# 5. Commit the removal + .gitignore change.
git add .gitignore
git commit -m "ops(fixtures): move demo fixture .mat to release asset fixtures-v0.1.0"
```

After the flip, a fresh clone has no fixture `.mat`. To provision:

```bash
# Download the asset, then extract + verify into demo/data:
/Applications/MATLAB_R2024b.app/bin/matlab -batch \
  "icemodel.verification.setup.fetchFixtures('v0.1.0', archive='icemodel-fixtures-v0.1.0.tar.gz', extract=true)"
```

In CI, run the equivalent before the suite (e.g. `gh release download
fixtures-v0.1.0 --pattern '*.tar.gz'` then `fetchFixtures(..., archive=..., extract=true)`).
The bootstrap hook then becomes the verifying provisioner instead of a no-op.

> The flip removes the fixtures from the working tree and future commits but not
> from past history; a separate history rewrite (`git filter-repo`) is optional
> and out of scope for the flip itself.
