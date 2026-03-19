function mustBeRollingBaselineName(arg)
   %MUSTBEROLLINGBASELINENAME Validate the mutable build-baseline selector.

   valid = icemodel.namelists.rollingbaseline();
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
