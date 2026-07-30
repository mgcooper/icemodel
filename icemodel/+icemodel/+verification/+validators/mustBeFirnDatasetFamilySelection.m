function mustBeFirnDatasetFamilySelection(arg)
   %MUSTBEFIRNDATASETFAMILYSELECTION Validate firn-family selectors plus "all".
   values = reshape(string(arg), 1, []);
   if isempty(values) || all(isblanktext(values))
      return
   end

   valid = ["all", icemodel.verification.namelists.firndatasetfamily()];
   bad = setdiff(values, valid, 'stable');
   if ~isempty(bad)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
