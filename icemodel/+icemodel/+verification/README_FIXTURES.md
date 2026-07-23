# v1.1 release-data provisioning

The v1.1 data boundary separates tracked demo data, provisioned formal/public
verification data, optional source-integration data, and the full interactive
scientific archive.

| Capability | Install root | Required |
|---|---|---|
| `formal-core` | `test/data` | yes |
| `verification-showcase` | `test/data` | yes |
| `forcing-integration` | `test/data/forcing` | no |

The complete interactive verification archive remains under top-level `data/`
and is not a release capability. `demo/data` remains tracked with only the two
15-minute demo forcing files and four spectral tables.

## Manifest

`test/assets/icemodel-v1.1-data-manifest.json` is the authoritative source. It
declares each archive and file's capability, required status, relative install
path, byte size, and SHA-256. The three release archive names are:

- `icemodel-v1.1-formal-core.tar.gz`
- `icemodel-v1.1-verification-showcase.tar.gz`
- `icemodel-v1.1-forcing-integration.tar.gz`

The manifest contains final local archive metadata and 141 optional
forcing-integration file rows. Publishing those artifacts remains a separate,
explicit approval gate.

## Producer

`packFixtures` makes one archive per selected capability and writes the filtered
release manifest beside them:

```matlab
result = icemodel.verification.setup.packFixtures("v1.1", ...
   capabilities=["formal-core", "verification-showcase", ...
      "forcing-integration"], ...
   root="/path/to/staged/test/data");
```

Packing refuses missing or hash-drifted source files. Output goes to the
gitignored `release-staging/` directory by default. On macOS, packing uses the
native USTAR writer with metadata copying disabled so undeclared AppleDouble
members cannot enter an archive.

## Consumer

Calling the provisioning API without overrides installs the two mandatory v1.1
capabilities and downloads missing release archives:

```matlab
result = icemodel.verification.setup.fetchFixtures("v1.1");
```

Pass `download=false` for network-free verification. Missing mandatory data then
reports this explicit provisioning command:

```matlab
result = icemodel.verification.setup.fetchFixtures("v1.1", download=false);

icemodel.verification.setup.fetchFixtures("v1.1", ...
   capabilities=["formal-core", "verification-showcase"], download=true)
```

Local archives and manifests support offline or pre-publication provisioning:

```matlab
result = icemodel.verification.setup.fetchFixtures("v1.1", ...
   capabilities="formal-core", ...
   manifest="/path/icemodel-v1.1-data-manifest.json", ...
   archive="/path/icemodel-v1.1-formal-core.tar.gz");
```

A scalar local archive requires one selected capability. Multiple selected
capabilities require one archive per capability in the same order.

Before canonical data are touched, `fetchFixtures` verifies the archive size and
SHA-256, checks raw tar headers, rejects unsafe paths/types and undeclared or
missing members, extracts to same-filesystem temporary storage, and verifies
every file. Promotion backs up declared prior paths and restores them on any
failure. Unrelated files are preserved, and an already-valid capability returns
success without rewriting data.
