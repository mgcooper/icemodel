function mustBeTestSmbmodelSelector(arg)
   %MUSTBETESTSMBMODELSELECTOR Validate one formal suite smbmodel selector.

   valid = icemodel.namelists.testsmbmodel();
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
