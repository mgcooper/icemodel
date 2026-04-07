# icemodel.surface.turbulence

Purpose: internal organization for surface turbulent-flux schemes.

Namespaces:
- `bulk_richardson`
- `monin_obukhov`

Public surface code should call:
- `icemodel.surface.diagnose_turbulent_heat_fluxes`

Not this layer directly, unless a test or a tightly scoped internal workflow needs a scheme-specific function.

Rules:
- keep conceptually parallel helpers visibly parallel across schemes
- use semantic names, not `br_*` prefixes
- keep comments/docstrings when moving or splitting helpers

Migration status: active and complete for the currently supported THF schemes.
