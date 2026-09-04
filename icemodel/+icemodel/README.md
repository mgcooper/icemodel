# icemodel

Use this namespace for everything needed to configure and run a model simulation.

- **Run configuration.** `setopts` builds defaults, `resetopts` applies caller
  overrides, `configureRun` resolves derived fields, and `getopts` reads them
  back. `config` and `setpath` resolve the workspace. `setcase` and
  `dependencies` wire a case and its external repositories.
- **Forcing.** `loadmet` loads met files, `processmet` derives runtime forcing
  fields, `interpmet` interpolates forcing values, and
  `createMetFileNames` builds the configured met-file list.
- **Output.** `postprocess` converts runtime values to model outputs.
  `buildOutputPayload`, `prepareRunOutput`, `updateoutput`, `writeoutput`,
  `concatoutput`, `retimeHourlyFixedStep`, and `outputYears` build, update,
  write, combine, retime, and select output records. `isIncrementChannel`
  identifies channels that aggregate by summation.
- **Restart.** `restartfile`, `saveRestartState`, `loadRestartState`.
- **Small runtime utilities** used across those steps: `pairedWindow`,
  `isPathInside`, `time2iter`, `chunkgridcell`, `parameterLookup`,
  `physicalConstant`.
- **Cross-namespace utility.** `shellQuote` quotes arguments for the git and
  cffconvert commands in `+internal` and for fixture and report commands in
  `+verification`.
