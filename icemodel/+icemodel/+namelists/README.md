# icemodel.namelists

Purpose: One function per list of supported names. Validators, argument
defaults, and output code read these lists, so each list has one definition.

Contents:

- Output channels:
  - `surfaceoutputs` returns the `ice1` channels.
  - `budgetoutputs` returns the mass-budget channels.
  - `cumulativeoutputs` returns the cumulative column diagnostics.
- Run selectors:
  - `smbmodel`
  - `solver`
  - `forcings`
  - `sitename`
  - `userdata`
  - `uservars`
  - `physicalconstant`
- Function selectors, each for the top-level function with the same name:
  - `completions`
  - `config`
  - `cvconvert`
  - `getpath`
- Test-suite selectors:
  - `testtier`
  - `testsmbmodel`
  - `testverbosity`
  - `unittest`
  - `benchmark`
  - `benchmarksamplingprofile`
  - `rollingbaseline`

`icemodel.validators` wraps these lists in argument validators.
