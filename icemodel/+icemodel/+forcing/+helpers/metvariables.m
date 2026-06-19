function [required, optional] = metvariables()
   %METVARIABLES Canonical met-file variable names for the forcing builders.
   %
   %  [required, optional] = icemodel.forcing.helpers.metvariables()
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
   %    ppt    [m]       total precipitation (water equivalent)
   %
   % The optional set covers source-specific diagnostics the builders
   % pass through when available (rain/snow precipitation split, melt,
   % runoff, turbulent fluxes, surface temperature, cloud fraction, snow
   % depth, wind direction).
   %
   % See also: icemodel.forcing.helpers.validatemet, icemodel.loadmet

   required = ["tair", "swd", "lwd", "albedo", "wspd", "rh", "psfc", "ppt"];

   optional = ["rainf", "snowf", "snow_depth", "melt", "runoff", ...
      "shf", "lhf", "tsfc", "cfrac", "snowd", "wdir"];
end
