function mustBeTestTierName(arg)
   %MUSTBETESTTIERNAME Validate one formal suite tier selector.

   valid = icemodel.namelists.testtier();
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
