# icemodel

Use this namespace for everything needed to configure and run a model simulation.
This file lists the top-level functions. Each subnamespace has its own README
that lists its members.

- **Run configuration.**
  - `setopts` builds default options.
  - `resetopts` applies caller overrides.
  - `configureRun` resolves derived fields.
  - `getopts` reads options back by name.
  - `saveRunOpts` saves the resolved `opts` for a run.
  - `config` and `getpath` resolve the workspace paths.
  - `mkfolders` creates the yearly output folders.
  - `setcase` and `dependencies` wire a case and its external repositories.
- **Forcing.**
  - `loadmet` loads met files.
  - `processmet` derives runtime forcing fields.
  - `interpmet` interpolates forcing values.
  - `createMetFileNames` builds the configured met-file list.
  - `resolvePrecipPhase` selects the runtime rain and snow split.
- **Output.**
  - `postprocess` converts runtime values to model outputs.
  - `buildOutputPayload` returns the output cell arrays in `opts` order.
  - `prepareRunOutput` prepares the output folders.
  - `updateoutput` stores one timestep of output.
  - `writeoutput` saves the output files.
  - `concatoutput` combines yearly output.
  - `retimeHourlyFixedStep` aggregates fixed-step data to hourly values.
  - `outputYears` selects the simulation years that saved output keeps.
  - `isIncrementChannel` identifies channels that aggregate by summation.
  - `loadresults` loads saved yearly output files and concatenates them.
  - `extractice2` extracts a time and depth subset of `ice2` data.
- **Restart.**
  - `restartfile` returns the restart file path for a run and year.
  - `saveRestartState` saves the year-boundary state.
  - `loadRestartState` loads it.
- **Small runtime utilities** used across those steps:
  - `pairedWindow` normalizes an optional start and end pair to UTC.
  - `isPathInside` checks that a path resolves inside a root.
  - `time2iter` converts a time to a met timestep index.
  - `chunkgridcell` returns the first and last cell of job 1 or job 2 of a
    two-job split. Its own comment marks it incomplete.
  - `parameterLookup` returns a model parameter value.
  - `physicalConstant` returns a physical constant.
  - `cvconvert` converts constituent amounts between control-volume forms,
    for example from mass to volume fraction.
- **Cross-namespace utilities.**
  - `shellQuote` quotes arguments for the git and cffconvert commands in
    `+internal` and for fixture and report commands in `+verification`.
  - `completions` returns the argument completion lists that the
    `icemodel.namelists` functions define.
