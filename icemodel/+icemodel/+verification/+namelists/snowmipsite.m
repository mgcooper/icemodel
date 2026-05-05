function sitenames = snowmipsite()
   %SNOWMIPSITE Canonical ESM-SnowMIP site-name namelist.
   %
   %  sitenames = icemodel.verification.namelists.snowmipsite()
   %
   % Outputs
   %  sitenames   String column of supported ESM-SnowMIP 3-letter site codes.
   %
   % Role
   %  Canonical site-name list for the "esm_snowmip" dataset family. This is
   %  the namelist consumed by argument validation and by case-id resolution.
   %  For the richer site catalog (display name, location, time-window
   %  metadata), use icemodel.verification.namelists.snowmipcatalog.

   % Order is alphabetical so namelists / case-id lists are stable across
   % releases.
   sitenames = [ ...
      "cdp"; "oas"; "obs"; "ojp"; "rme"; ...
      "sap"; "snb"; "sod"; "swa"; "wfj"];
end
