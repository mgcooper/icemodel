# icemodel.timestepping

Purpose: full-step/substep control, retry/reset flow, and timestep adaptation.

Current contents:
- `initialize_timesteps`
- `newtimestep`
- `checksubstep`
  - accept, retry, or force-advance one substep. Owns the cross-step
    forced-advance streak guard (one full forcing step of consecutive
    forced time is the limit)
- `resetsubstep`
- `acceptsubstep`
  - accept one substep: checkpoint the state pass-through and credit the
    substep time. A forced advance calls it with the restored checkpoint,
    accepting elapsed time only
- `nexttimestep`
- `getforcings` and `getsubstepforcings`
  - legacy met-struct scalarizers with no production caller

Rules:
- keep runtime control logic here, even when state payloads include column
  variables
- do not move generic math into this namespace

This namespace holds the active runtime control helpers.
