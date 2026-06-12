function mustBeSnowmipSite(sitename)
   %MUSTBESNOWMIPSITE Validate sitename against the canonical ESM-SnowMIP namelist.
   %
   %  Used in arguments blocks where the namelist function cannot be
   %  called directly (MATLAB validation funcs allow only literals or
   %  previously declared args).

   valid = icemodel.verification.namelists.snowmipsite();
   if ~ismember(sitename, valid)
      error('icemodel:verification:validators:mustBeSnowmipSite', ...
         'unknown ESM-SnowMIP sitename %s. Valid: %s', ...
         sitename, strjoin(valid, ', '));
   end
end
