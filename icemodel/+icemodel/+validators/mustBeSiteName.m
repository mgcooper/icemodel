function mustBeSiteName(arg)
   %MUSTBESITENAME Validate that input is a valid core point-run site name.

   valid = icemodel.namelists.sitename();
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
