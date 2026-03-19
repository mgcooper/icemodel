function mustBeForcingName(arg)
   %MUSTBEFORCINGNAME Validate that input is a valid forcing-source name.

   valid = icemodel.namelists.forcings();
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
