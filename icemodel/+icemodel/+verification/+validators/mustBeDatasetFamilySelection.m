function mustBeDatasetFamilySelection(arg)
   %MUSTBEDATASETFAMILYSELECTION Validate dataset-family selectors plus "all".
   %
   %  icemodel.verification.validators.mustBeDatasetFamilySelection(arg)
   %
   % Plotting and preview helpers accept vectors of canonical family ids and the
   % convenience selector "all".

   values = reshape(string(arg), 1, []);
   if isempty(values) || all(isblanktext(values))
      return
   end

   valid = ["all"; icemodel.verification.namelists.datasetfamily()];
   bad = setdiff(values, valid, 'stable');
   if ~isempty(bad)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
