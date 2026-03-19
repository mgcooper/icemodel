function mustBeBenchmarkSamplingProfileName(arg)
   %MUSTBEBENCHMARKSAMPLINGPROFILENAME Validate one benchmark profile name.

   valid = icemodel.namelists.benchmarksamplingprofile();
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
