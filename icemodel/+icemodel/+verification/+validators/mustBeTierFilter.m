function mustBeTierFilter(arg)
   %MUSTBETIERFILTER Validate one optional snow-verification tier selector.
   %
   %  icemodel.verification.validators.mustBeTierFilter(arg)
   %
   % Inputs
   %  arg   Empty value or one tier selector.
   %
   % Role
   %  Argument-block validator for optional tier filters.

   % Blank means "all tiers" in callers, so it is valid here.
   if isempty(arg) || isblanktext(arg)
      return
   end

   % Compare against the canonical namelist so runner and case filters stay
   % synchronized.
   valid = icemodel.verification.namelists.tier();
   if ~ismember(arg, valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
