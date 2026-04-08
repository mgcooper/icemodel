# icemodel.radiation

Purpose: optical properties and radiative-transfer machinery.

Public entry points:
- `icemodel.radiation.initialize_spectral_model`
- `icemodel.radiation.get_scattering_coefficients`
- `icemodel.radiation.get_solar_spectrum`
- `icemodel.radiation.update_extinction_coefficients`
- `icemodel.radiation.spectral_extinction_coefficients`
- `icemodel.radiation.rescale_spectral_extinction_coefficients`
- `icemodel.radiation.bulk_extinction_coefficients`
- `icemodel.radiation.bulk_extinction_coefficients_lookup`
- `icemodel.radiation.solvetwostream`
- `icemodel.radiation.smooth_twostream_fluxes`

Contents:
- extinction and scattering property assembly
- solar forcing helpers
- two-stream solver machinery

Rules:
- keep general radiation algorithms here
- keep column-owned absorbed shortwave assembly in `icemodel.column`

Migration status: active optical and two-stream helpers now live here.
