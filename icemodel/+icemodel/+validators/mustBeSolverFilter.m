function mustBeSolverFilter(arg)
   %MUSTBESOLVERFILTER Validate an optional solver-filter vector.
   %
   %  Valid values are the subset of `icemodel.namelists.solver()`. Empty is
   %  accepted so callers can express "no filter".

   if isempty(arg)
      return
   end

   mustBeNumeric(arg)
   mustBeInteger(arg)

   valid = icemodel.namelists.solver();
   if any(~ismember(arg, valid))
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('solver must be empty or a subset of [%s]', ...
         strjoin(string(valid), ' '));
      throwAsCaller(MException(eid, msg));
   end
end
