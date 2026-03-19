function mustBeTestVerbosityName(arg)
   %MUSTBETESTVERBOSITYNAME Validate one unit-runner verbosity selector.

   valid = icemodel.namelists.testverbosity();
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
