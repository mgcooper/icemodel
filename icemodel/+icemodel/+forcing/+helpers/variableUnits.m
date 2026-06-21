function units = variableUnits(names)
   %VARIABLEUNITS Canonical unit string for each forcing-builder channel.
   %
   %  units = icemodel.forcing.helpers.variableUnits(names)
   %
   % Maps every channel emitted by any icemodel.forcing builder (MAR,
   % MERRA-2, RACMO, PROMICE) to one canonical unit, so the Data/met
   % outputs of every source agree on units. NAMES is a string array of
   % channel names; UNITS is a cellstr of the matching unit strings, ready
   % to assign to a timetable's Properties.VariableUnits.
   %
   % Canonical unit choices (one per channel, source-independent):
   %  tair, tsfc                   K        air / surface temperature
   %  swd, swu, lwd, lwu, swn, lwn,
   %    netr, shf, lhf, thf        W m-2    radiative and turbulent fluxes
   %  albedo, modis, cfrac         -        fractions (0-1)
   %  rh                           %        relative humidity (percent)
   %  wspd                         m s-1    wind speed
   %  wdir                         degrees  wind direction
   %  psfc                         Pa       surface air pressure
   %  ppt, rainf, snowf, rain,
   %    snow, precip               m s-1    precipitation as a water-
   %                                        equivalent RATE (canonical met
   %                                        precip unit; see metvariables)
   %  melt, runoff, evap, smb,
   %    refreeze, subl, sndiv,
   %    meltin                     mWE/h    diagnostic mass-flux rates
   %                                        (metres water equiv. per hour)
   %  swe                          kg m-2   snow water-equivalent store
   %  snowd, snow_depth, ablation,
   %    boom_height, stake_height,
   %    transducer_depth, elev     m        heights / depths
   %  tice1..ticeN                 K        ice-temperature string
   %
   % An emitted channel missing from the map is an error: the builders must
   % never ship an unlabeled column. Add the channel here when a builder
   % gains a new output.
   %
   % See also: icemodel.forcing.helpers.metvariables,
   %  icemodel.forcing.data2met

   arguments
      names (1, :) string
   end

   unitmap = canonicalUnitMap();

   units = strings(1, numel(names));
   for k = 1:numel(names)
      name = names(k);
      if isKey(unitmap, name)
         units(k) = unitmap(name);
      elseif ~isempty(regexp(name, '^tice\d+$', 'once'))
         units(k) = "K";   % ice-temperature string of arbitrary length
      else
         error('icemodel:forcing:variableUnits:unmappedChannel', ...
            ['no canonical unit for channel "%s"; add it to ' ...
            'icemodel.forcing.helpers.variableUnits'], name);
      end
   end
   units = cellstr(units);
end

%% Local functions
function m = canonicalUnitMap()
   %CANONICALUNITMAP The channel -> canonical unit dictionary.
   keys = [ ...
      "tair", "tsfc", ...
      "swd", "swu", "lwd", "lwu", "swn", "lwn", "netr", ...
      "shf", "lhf", "thf", ...
      "albedo", "modis", "cfrac", ...
      "rh", "wspd", "wdir", "psfc", ...
      "ppt", "rainf", "snowf", "rain", "snow", "precip", ...
      "melt", "runoff", "evap", "smb", "refreeze", "subl", ...
      "sndiv", "meltin", "swe", ...
      "snowd", "snow_depth", "ablation", "boom_height", ...
      "stake_height", "transducer_depth", "elev"];
   vals = [ ...
      "K", "K", ...
      "W m-2", "W m-2", "W m-2", "W m-2", "W m-2", "W m-2", "W m-2", ...
      "W m-2", "W m-2", "W m-2", ...
      "-", "-", "-", ...
      "%", "m s-1", "degrees", "Pa", ...
      "m s-1", "m s-1", "m s-1", "m s-1", "m s-1", "m s-1", ...
      "mWE/h", "mWE/h", "mWE/h", "mWE/h", "mWE/h", "mWE/h", ...
      "mWE/h", "mWE/h", "kg m-2", ...
      "m", "m", "m", "m", ...
      "m", "m", "m"];
   m = dictionary(keys, vals);
end
