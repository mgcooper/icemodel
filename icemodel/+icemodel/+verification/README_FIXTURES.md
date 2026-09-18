# Release-data provisioning

Release data live in four separate trees: tracked demo data, provisioned
formal and public verification data, optional source-integration data, and the
full scientific archive.

| Capability | Install root | Required |
|---|---|---|
| `formal-core` | `test/data` | yes |
| `verification-showcase` | `test/data` | yes |
| `forcing-integration` | `test/data/forcing` | no |

The full interactive verification archive remains under top-level `data/`
and is not a release capability. `demo/data` remains tracked with only the two
15-minute demo forcing files and four spectral tables.

## Manifest

The tracked `test/assets/icemodel-<version>-data-manifest.json` governs release
data. It declares each archive and file's capability, required status, relative
install path, byte size, and SHA-256. The three v1.3 archive names are:

- `icemodel-v1.3-formal-core.tar.gz`
- `icemodel-v1.3-verification-showcase.tar.gz`
- `icemodel-v1.3-forcing-integration.tar.gz`

The manifest contains final local archive metadata and 141 optional
forcing-integration file rows. Publishing those artifacts remains a separate,
explicit approval gate. Published releases retain their versioned manifests
and archives.

## Producer

`packFixtures` makes one archive per selected capability and writes the filtered
release manifest beside them:

```matlab
result = icemodel.verification.setup.packFixtures( ...
   capabilities=["formal-core", "verification-showcase", ...
      "forcing-integration"], ...
   root="/path/to/staged/test/data");
```

Packing rejects a missing source file or a file whose hash does not match the
manifest. By default, it writes output to the gitignored `release-staging/`
directory. On macOS, packing uses the native USTAR writer with metadata copying
disabled so undeclared AppleDouble members cannot enter an archive.

## Consumer

Calling the provisioning API without a version installs the two mandatory
capabilities for the version in `CITATION.cff` and downloads missing release
archives:

```matlab
result = icemodel.verification.setup.fetchFixtures();
```

Pass `download=false` for network-free verification. Missing mandatory data then
reports this explicit provisioning command:

```matlab
result = icemodel.verification.setup.fetchFixtures(download=false);

icemodel.verification.setup.fetchFixtures( ...
   capabilities=["formal-core", "verification-showcase"], download=true)
```

Pass an explicit version to reproduce an earlier release or to use local
archives and manifests before publication:

```matlab
result = icemodel.verification.setup.fetchFixtures("v1.1", ...
   capabilities="formal-core", ...
   manifest="/path/icemodel-v1.1-data-manifest.json", ...
   archive="/path/icemodel-v1.1-formal-core.tar.gz");
```

A scalar local archive requires one selected capability. Multiple selected
capabilities require one archive per capability in the same order.

`fetchFixtures` finishes every check before it changes installed data. For
each selected archive that needs installation it checks the size, the SHA-256,
and the raw tar headers. It rejects unsafe paths and types, undeclared
members, and missing members. Extraction goes to temporary storage on the same
filesystem, where every file is verified before promotion. Promotion backs up
the declared existing paths and restores them after any failure, and it leaves
unrelated files alone. A capability that is already valid returns success
without rewriting data.
