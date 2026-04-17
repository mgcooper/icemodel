# icemodel.timestepping

Purpose: full-step/substep control, retry/reset flow, and timestep adaptation.

Current contents:
- `check_substep`
- `reset_substep`
- `update_substep`
- `init_timesteps`
- `new_timestep`
- `next_step`
- `update_force_advance_guard`

Rules:
- keep runtime control logic here, even when state payloads include column variables
- do not move generic math into this namespace

Migration status: active runtime control helpers have been migrated here.
