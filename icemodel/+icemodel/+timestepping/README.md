# icemodel.timestepping

Purpose: Full-step and substep control, retries, resets, and timestep changes.

Contents:

- `initialize_timesteps`
- `newtimestep`
- `checksubstep`
  - Accepts, retries, or force-advances one substep. Reads the attempt flags
    from `diag.substep`. If both solves succeed, the outer loop fails, and
    settings are not in recovery mode, `checksubstep` retries once. The retry
    uses `cpl_recovery_alpha` and disables acceleration. If settings already
    match recovery mode, `checksubstep` shortens dt. It also shortens dt after
    a failed retry. After `maxsubstep` failures, it forces an advance. It counts
    `diag.n_failed_substeps` and `diag.n_forced_advances`.
- `resetsubstep`
- `acceptsubstep`
  - Accepts one substep: checkpoints the state, credits the substep time, adds
    the accepted `diag.substep` to the timestep record with
    `icemodel.couplers.update_solver_diag`, and restores the primary
    `settings0`. A forced advance calls it with the restored checkpoint and
    accepts elapsed time only.
- `nexttimestep`
- `getforcings` and `getsubstepforcings`
  - Legacy forcing-access functions.
