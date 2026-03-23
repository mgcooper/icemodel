function mustBeFormalSmbmodelName(arg)
   %MUSTBEFORMALSMBMODELNAME Validate one concrete formal-suite smbmodel.

   valid = icemodel.namelists.smbmodel("test");
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
