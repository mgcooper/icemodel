# icemodel.run

Purpose: Convenience entry points that set options, run one model, and
return postprocessed output.

Contents:

- `point` runs `icemodel` or `skinmodel` for one site.
  - Its name-value inputs:
    - `sitename`, the site
    - `forcings`, the forcing source
    - `userdata` and `uservars`, the userdata source and variables
    - `smbmodel`, `icemodel` or `skinmodel`
    - `simyears` and `n_spinup_years`, the simulation and spinup years
    - `gridcell`, an optional sector grid cell
    - `testname`, `saveflag`, and `backupflag`, the output options
  - It returns:
    - `ice1`
    - `ice2`
    - `met`
    - the resolved `opts`
  - With `saveflag=true` it loads the saved output with
    `icemodel.loadresults`. Otherwise it calls `icemodel.postprocess`.
