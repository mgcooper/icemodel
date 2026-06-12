function mustBeDatasetFamilyFilter(arg)
   %MUSTBEDATASETFAMILYFILTER Validate one optional dataset-family selector.
   %
   %  icemodel.verification.validators.mustBeDatasetFamilyFilter(arg)
   %
   % Inputs
   %  arg   Empty value or one dataset-family selector.
   %
   % Role
   %  Argument-block validator for listcases dataset-family filtering.

   % Blank means "all families" in callers, so it is valid here.
   if isempty(arg) || isblanktext(arg)
      return
   end

   % Compare against the canonical namelist so setup and workflow selectors stay
   % synchronized.
   valid = icemodel.verification.namelists.datasetfamily();
   if ~ismember(arg, valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
