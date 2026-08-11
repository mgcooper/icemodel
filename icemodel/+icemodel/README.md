# icemodel

The top-level namespace is the runtime layer. It holds everything needed to
configure a run, feed it, run it, and get results back out. The sub-namespaces
(`+column`, `+surface`, `+couplers`, `+vapor`, `+radiation`, `+numerics`) hold
the physics. This level orchestrates them.

## What belongs here

A function belongs at this level when it is part of running the model rather
than part of solving the physics:

- **Run configuration.** `setopts` builds defaults, `resetopts` applies caller
  overrides, `configureRun` resolves derived fields, and `getopts` reads them
  back. `config` and `setpath` resolve the workspace. `setcase` and
  `dependencies` wire a case and its external repositories.
- **Forcing.** `loadmet` is the canonical met-loading layer, with `processmet`,
  `interpmet`, and `createMetFileNames` around it.
- **Output.** `postprocess` is the canonical model-output layer.
  `buildOutputPayload`, `prepareRunOutput`, `updateoutput`, `writeoutput`,
  `concatoutput`, `retimeHourlyFixedStep`, and `outputYears` shape and persist
  results. `isIncrementChannel` decides how a channel aggregates.
- **Restart.** `restartfile`, `saveRestartState`, `loadRestartState`.
- **Small runtime utilities** used across those steps: `pairedWindow`,
  `isPathInside`, `time2iter`, `chunkgridcell`, `parameterLookup`,
  `physicalConstant`.
- **Utilities shared by namespaces that must not depend on each other.**
  `shellQuote` is the one such case: `+internal` needs it to build git and
  cffconvert commands, and `+verification` needs it for fixture packing and
  the report builders. Neither namespace may depend on the other, so their
  shared helper sits here.

## What does not belong here

- Physics kernels. Those go in the namespace that owns the process.
- Fixed lists of supported names, ids, or channels. Those are namelists and go
  in `+namelists`.
- Toolbox management, such as version handling and install support. That is
  `+internal`.
- Verification, staging, and reporting. Those are `+verification`.
