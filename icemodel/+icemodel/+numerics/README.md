# icemodel.numerics

Purpose: generic math algorithms shared across model domains.

Current contents:
- `trisolve`
- `aitkenscalar`
- `fsearchzero`
- `complexstep`
- `complex_step_derivative`
- `sign_or_one`

Rules:
- only truly domain-generic methods belong here
- surface or column solve orchestration stays in the owning physics namespace

Migration status: the active shared numerics helpers have been migrated here.
