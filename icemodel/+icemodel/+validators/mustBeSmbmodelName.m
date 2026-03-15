function mustBeSmbmodelName(arg)
%MUSTBESMBMODELNAME Validate that input is a valid core smbmodel name.

   valid = icemodel.namelists.smbmodel();
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
