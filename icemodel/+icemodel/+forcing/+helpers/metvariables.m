function [required, optional, pptunit] = metvariables()
   %METVARIABLES Canonical met-file variable names for the forcing builders.
   %
   %  [required, optional] = icemodel.forcing.helpers.metvariables()
   %  [required, optional, pptunit] = icemodel.forcing.helpers.metvariables()
   %
   % Returns the variable names of the icemodel met-file contract as
   % string arrays. The required set is what icemodel.loadmet and the
   % model need to run (units in brackets):
   %
   %    tair   [K]       air temperature
   %    swd    [W m-2]   downwelling shortwave radiation
   %    lwd    [W m-2]   downwelling longwave radiation
   %    albedo [-]       surface albedo
   %    wspd   [m s-1]   wind speed
   %    rh     [%]       relative humidity
   %    psfc   [Pa]      surface air pressure
   %    ppt    [m s-1]   total precipitation as a water-equivalent RATE
   %
   % Canonical precipitation unit (PPTUNIT, third output): the precipitation
   % channels (ppt, and the optional rainf/snowf split) are a water-equivalent
   % RATE in metres per second, the unit icemodel.surface.advective_heat_flux
   % consumes directly. This is timestep-independent, so it is consistent
   % across every source and sampling interval. ESM-SnowMIP already produces
   % m s-1 (mass flux / ro_liq); the gridded sources (MAR/MERRA/RACMO) emit
   % their precipitation channels directly in m s-1 from the Data builders
   % (the source mWE/h posting is converted inside the builder), so the unit
   % is uniform across every source. The shared canonical unit map
   % icemodel.forcing.helpers.variableUnits records m s-1 for them.
   %
   % The optional set covers source-specific diagnostics the builders
   % pass through when available (rain/snow precipitation split, melt,
   % runoff, surface mass balance, turbulent fluxes, surface temperature,
   % cloud fraction, snow depth, PROMICE boom height, wind direction, and
   % optional MODIS albedo).
   %
   % See also: icemodel.forcing.helpers.validatemet, icemodel.forcing.data2met,
   %  icemodel.loadmet

   required = ["tair", "swd", "lwd", "albedo", "wspd", "rh", "psfc", "ppt"];

   optional = ["rainf", "snowf", "snow_depth", "melt", "runoff", "smb", ...
      "shf", "lhf", "tsfc", "cfrac", "snowd", "boom_height", "wdir", ...
      "modis"];

   % Canonical precipitation-rate unit (water-equivalent metres per second).
   pptunit = "m s-1";
end
