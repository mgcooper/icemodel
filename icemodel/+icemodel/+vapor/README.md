# icemodel.vapor

Purpose: vapor thermodynamics and atmosphere-surface vapor conversions shared across domains.

Planned contents:
- saturation vapor pressure relations
- RH/dew-point/wet-bulb conversions
- vapor density, diffusivity, and conductivity helpers

Rules:
- keep general thermodynamic transforms here
- keep column microstructure transport, such as `VAPORTRANSFER`, in `icemodel.column`

Migration status: planned namespace scaffold only; file migration still open.
