function mustBeCaseIdSubset(arg)
   %MUSTBECASEIDSUBSET Validate a subset of canonical snow-verification case ids.
   %
   %  icemodel.verification.validators.mustBeCaseIdSubset(arg)
   %
   % Inputs
   %  arg   Empty value or string array of requested case ids.
   %
   % Role
   %  Argument-block validator for user-facing verification workflow functions.

   % Empty means "all cases" in callers, so it is valid here.
   if isempty(arg)
      return
   end

   % Compare against the canonical namelist so runner and setup selectors stay
   % synchronized.
   valid = icemodel.verification.namelists.caseid();
   invalid = setdiff(reshape(arg, [], 1), valid);
   if ~isempty(invalid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Values must be chosen from:\n%s', ...
         strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
