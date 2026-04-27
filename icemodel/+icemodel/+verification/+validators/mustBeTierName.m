function mustBeTierName(arg)
   %MUSTBETIERNAME Validate one snow-verification tier selector.
   %
   %  icemodel.verification.validators.mustBeTierName(arg)
   %
   % Inputs
   %  arg   One required tier selector.
   %
   % Role
   %  Argument-block validator for required runner tier options.

   % Required selectors cannot be blank; they must match the canonical namelist.
   valid = icemodel.verification.namelists.tier();
   if ~ismember(arg, valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
