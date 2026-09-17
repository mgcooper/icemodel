# icemodel.validators

Purpose: Argument-block validators. Each validator checks its input against
the matching `icemodel.namelists` function, so the allowed names have one
definition.

Contents:

- Run inputs:
  - `mustBeSmbmodelName`
  - `mustBeSiteName`
  - `mustBeForcingName`
  - `mustBeUserdataName`
  - `mustBeUservarName`
- Formal test-suite selectors:
  - `mustBeFormalSmbmodelName`
  - `mustBeTestSmbmodelSelector`
  - `mustBeTestTierName`
  - `mustBeSolverFilter`
  - `mustBeRollingBaselineName`
  - `mustBeTestVerbosityName`
  - `mustBeBenchmarkSamplingProfileName`
