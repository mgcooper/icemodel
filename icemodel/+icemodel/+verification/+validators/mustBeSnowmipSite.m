function mustBeSnowmipSite(sitename)
   %MUSTBESNOWMIPSITE Validate sitename against the canonical ESM-SnowMIP
   % namelist.
   %
   %  Use this validator in arguments blocks that cannot call the namelist
   %  function directly. A MATLAB validation function accepts only literals
   %  or arguments declared earlier in the same block.

   valid = icemodel.verification.namelists.snowmipsite();
   if ~ismember(sitename, valid)
      error('icemodel:verification:validators:mustBeSnowmipSite', ...
         'unknown ESM-SnowMIP sitename %s. Valid: %s', ...
         sitename, strjoin(valid, ', '));
   end
end
